# Hosting SequencePlatform on Tencent Cloud Lighthouse

How the platform (`SequencePlatform.*`) is hosted publicly: one Tencent Cloud Lighthouse instance in Frankfurt runs the Compose stack
that the AppHost generates, with Caddy in front for HTTPS. The files next to this README are the server's configuration, with
addresses replaced by placeholders. Set up 2026-11 from commit `86c7736`.

| File | Goes to |
|---|---|
| [`Caddyfile`](Caddyfile) | `/etc/caddy/Caddyfile` |
| [`cat-platform.service`](cat-platform.service) | `/etc/systemd/system/cat-platform.service` |
| [`sshd-01-hardening.conf`](sshd-01-hardening.conf) | `/etc/ssh/sshd_config.d/01-hardening.conf` |
| [`fail2ban-sshd.local`](fail2ban-sshd.local) | `/etc/fail2ban/jail.d/sshd.local` |
| [`seed-database.md`](seed-database.md) | Instructions: copy a local catalogue to the server |

## Why this setup

The stack is PostgreSQL 17, a one-off `migrate`, the API with the web UI, and exactly one search worker. It's CPU-only, uses about
450 MB of RAM in total, and its database holds profiles, never DNA (about 0.5 GB for the 2026-10 catalogue; 11 GB since
2026-11-06, with the gene layer). With that catalogue the API peaks at about 0.8 GB, during localization (measured locally).

| Choice | Why |
|---|---|
| **Lighthouse**, not CVM | Same CPU and RAM for much less, and the price includes the public IP, an SSD and 2.5 TB/month of traffic. Nothing here needs a VPC, a load balancer or autoscaling. |
| **Not serverless** | The worker is long-running, listens for PostgreSQL notifications (`search_queued`) and holds in-memory indexes, and there must be only one. Serverless also needs a managed PostgreSQL that costs more than the whole server. |
| **2 vCPU / 8 GB / 80 GB** (`bundle_starter_nmc_lin_med8_01`) | 4 GB serves fine, but leaves little room for image builds, the Ingest's maintenance jobs (`build-null`, `build-neighbours`) and the optional localization cache (+823 MB). |
| **Frankfurt** (`eu-frankfurt`) | Close to the users, inside the EU, and no ICP filing (outside mainland China). |
| **Image `lhbp-2cacsycc`** | Ubuntu Server 24.04 LTS with Docker preinstalled. |
| **Monthly, manual renewal** | The 12-month price was exactly 12 × the monthly one, so prepaying saved nothing. |

Price at purchase: **$7.50/month** (75% off the $10.00 list price; assume the list price at renewal). On this machine (Xeon Gold 6133,
AVX2): idle memory api ~56 MB, worker ~38 MB, postgres ~51 MB, and the P4 end-to-end scenario passes in 1.7 s over HTTPS.

## Network

| Port | Open to | Why |
|---|---|---|
| 22/tcp | everyone | SSH, key-only, with fail2ban. Narrow it to admin IPs in the Lighthouse firewall if they're stable. |
| 80/tcp | everyone | Let's Encrypt HTTP challenge, then a redirect to HTTPS |
| 443/tcp | everyone | The platform |
| ICMP | everyone | ping |

The API (5080) and PostgreSQL (5432) are never published beyond the machine: the generated compose file binds the API to
`127.0.0.1:5080`, and Postgres only exposes 5432 on the Compose network.

## HTTPS without a hostname

The platform is served on the bare IP. Caddy (2.10 or later, from Caddy's own apt repository; Ubuntu's 2.6 is too old) gets a
**Let's Encrypt certificate for the IP address** through the ACME `shortlived` profile, valid about 6 days and renewed automatically.

Clients send no SNI when they connect to an IP, so without `default_sni` Caddy can't choose the certificate and the TLS handshake
fails with an internal-error alert (`openssl s_client -servername IP` still works, which hides the problem).

This is the S6.4 step the API's README waits for: bearer keys now travel only over TLS, and Kestrel stays on localhost.

**The API's database role (S6.4, since 2026-11-10).** The `api` and `worker` containers connect as `sequence_app`, not `postgres`. Migration
0021 gives that role read access to every table and write access only to what uploads and search jobs write. It can't change the schema,
write the computed layers or the users, or run programs on the server (`COPY ... PROGRAM`). The deploy generates its password as
`SEQUENCE_APP_PASSWORD` in `.env.Production`, beside `POSTGRES_PASSWORD`, and `migrate` (still `postgres`) gives the role its login on
every start. Admin commands through the `migrate` image keep running as `postgres`. If the API logs `permission denied for table ...`, a
migration added a table the API writes without granting it: grant it in a new migration, not by hand.

**Compression is Caddy's job** (`encode zstd gzip`, since 2026-11-06). The API's container sets `ASPNETCORE_FORWARDEDHEADERS_ENABLED=true`,
so a request Caddy forwards with `X-Forwarded-Proto: https` counts as HTTPS, and ASP.NET Core's response compression is off for HTTPS by
default: until then every API response left the server uncompressed (E. coli K-12's gene list: 760 KB; with `encode`, 157 KB zstd or
162 KB gzip). Static files still arrive Brotli-compressed by the API, and Caddy leaves an already-encoded response alone. BREACH, the
reason compression over TLS needs thought, takes a secret and attacker-controlled text in one response; the API never echoes the key.

## Setting up a new server

1. **Instance.** Create it from the console or `tccli lighthouse CreateInstances` with the bundle and image above and an imported SSH
   key pair (`ImportKeyPair`), so no password is ever set. Run it with `--DryRun true` first.
2. **Firewall.** Set the Lighthouse firewall as in [Network](#network). The default opens 22 and 80 but not 443.
3. **OS.**
   - `apt-get full-upgrade`, enable `unattended-upgrades`, then reboot if `/var/run/reboot-required` exists.
   - Set the clock to UTC (`timedatectl set-timezone UTC`). Lighthouse images default to Asia/Shanghai.
   - The image already has a 1.9 GB `/swap.img`. Set `vm.swappiness=10`.
4. **SSH and fail2ban.** Install [`sshd-01-hardening.conf`](sshd-01-hardening.conf), check it with `sshd -t`, then
   `systemctl reload ssh`. Test a new login before closing the current session. Then `apt-get install fail2ban` and install
   [`fail2ban-sshd.local`](fail2ban-sshd.local).
5. **Build tools.**
   - `apt-get install dotnet-sdk-10.0 git`: Ubuntu's own .NET 10 package, kept current by `unattended-upgrades`.
   - Add the admin users to the `docker` group.
6. **Platform.**

   ```bash
   sudo install -d -o ubuntu -g ubuntu /opt/cat && cd /opt/cat
   git clone --branch feature/localPlatfromBuild https://github.com/HristoS/CATSequenceAnalysis.git && cd CATSequenceAnalysis
   dotnet run --project SequencePlatform.AppHost --no-launch-profile -- --operation publish --step deploy --output-path deploy/compose
   ```

   - About 2 minutes on 2 vCPU.
   - `--output-path` is resolved against the **AppHost project**, so the files land in
     `SequencePlatform.AppHost/deploy/compose/`.
   - The live env file is **`.env.Production`**: image tags plus the Postgres password. Run
     `chmod 600 .env* ~/.aspire/deployments/*/production.json*`. Aspire also keeps the password in that deployment state.
   - `docker compose ls` shows the project name the service file needs.
7. **Start at boot.** The generated services have no restart policy. Install
   [`cat-platform.service`](cat-platform.service) with the project name filled in, then `systemctl enable --now cat-platform`.
8. **HTTPS.** Install Caddy from Caddy's apt repository (<https://caddyserver.com/docs/install#debian-ubuntu-raspbian>), then
   install the [`Caddyfile`](Caddyfile) with the IP filled in and run `systemctl reload caddy`.
9. **Check.**
   - `curl https://SERVER_IP/health` should return `Healthy`.
   - Run `python3 -I scripts/platform-checks/e2e_scenario.py prepare`, then `run https://SERVER_IP KEY_FILE`, against an
     **empty** database with a throwaway user.
   - Disable that user afterwards.

## Operating it

Admin commands go through the `migrate` image, from `SequencePlatform.AppHost/deploy/compose/`:

```bash
docker compose -p COMPOSE_PROJECT -f docker-compose.yaml --env-file .env.Production run --rm --no-deps -T migrate <command>
```

- `<command>` is any Ingest command: `add-user NAME`, `disable-user NAME`, `stats`, `build-null` and so on.
- `add-user` prints the key once. Capture it straight into a `chmod 600` file rather than a terminal or chat.
- The `libgssapi_krb5.so.2` error these print is harmless: it's Npgsql's optional Kerberos support, which the chiseled image doesn't ship.

**Upgrading:**
1. Take a Lighthouse snapshot.
2. `git pull`.
3. Re-run the deploy command from step 6. `migrate` applies the new migrations, which are append-only. The first deploy with migration
   0021 adds `SEQUENCE_APP_PASSWORD` to `.env.Production`: run the `chmod 600` from step 6 again.
4. Run `systemctl restart cat-platform`.

**Derived tables after a deploy that adds them:** migration 0019 (protein p and E) needs one run of
`<command>` = `build-protein-null` (6 min on a 4-core workstation over the 2026-11 catalogue; expect several times that on 2 vCPU). It
fetches about 780 proteins from BV-BRC and stores only score counts. Until it has run, protein hits show no E.

**Seeding or replacing the catalogue:** see [`seed-database.md`](seed-database.md).

**Adding a layer computed locally (P6 S6.5).** Don't dump and restore for this: a restore replaces the server's own users, uploads and
search log. Since the 2026-11-06 restore, the server's records have the ids the local ones have, so a layer that only adds rows to
existing records can be rebuilt on the server by the same command from the same commit, and then compared. For the eukaryote gene layer
(S7.T):

1. Deploy the commit as in **Upgrading**; `migrate` applies migration 0020.
2. Run the pass through the `migrate` image (about 1.3 GB from NCBI and an hour for the 131 annotated assemblies, extrapolated):
   `... run --rm --no-deps -T migrate backfill-genes-ncbi --contact EMAIL`.
3. Run this on both databases; the lines must match. If a `cds_sha256` differs, NCBI re-released that annotation between the two runs:
   rerun the local pass for it (delete its `protein_kmer` rows first) rather than accept the difference.

   ```sql
   select s.accession, r.annotations #>> '{ncbi,gene_layer,cds_sha256}' as cds_sha256, count(g.*) as rows, sum(g.byte_size) as bytes,
          md5(string_agg(encode(g.sha256, 'hex'), '' order by g.bucket)) as rows_md5
   from sequence_record r join sequence_source s on s.id = r.source_id left join gene_profile g on g.sequence_id = r.id
   where r.kind = 'eukaryote' group by 1, 2 order by 1;
   ```
4. Rewrite the gene table in bucket order, which protein search reads faster (`results/GeneProfileOrder_Check.csv`):
   `... run --rm --no-deps -T migrate cluster-genes`. Check `df -h` first: it needs free disk about the table's size (6.8 GB locally).
   It holds an exclusive lock for the rewrite (293 s locally, longer on 2 vCPU), so protein searches wait meanwhile.

New records created locally are a different matter: their ids can collide with uploads made on the server, so they go to the server by
ingesting there, not by copying rows.

**The restricted deploy key (since 2026-11-11).** Claude Code deploys through a separate key that can only run
[`cat-deploy.sh`](cat-deploy.sh): `status`, `deploy <commit>`, `rollback`, `logs`, read-only `admin` queries, the maintenance
`job`s, and `user add|rotate|disable` (API keys go to a mode-600 file and are never printed). Arguments are split on whitespace,
matched against fixed patterns and never evaluated. Every call is logged to `journalctl -t cat-deploy`. Locally, the
`deploy-platform` skill drives it, and a Claude Code hook stops sessions from reading SSH keys, `docs/tencent` or transcripts.

On your machine (PowerShell):

```powershell
ssh-keygen -t ed25519 -f $HOME\.ssh\cat_deploy_ed25519 -C claude-deploy        # set a passphrase
Get-Service ssh-agent | Set-Service -StartupType Automatic; Start-Service ssh-agent   # once, as administrator
ssh-add $HOME\.ssh\cat_deploy_ed25519
```

and in `~/.ssh/config`:

```
Host cat-deploy
    HostName SERVER_IP
    User ubuntu
    IdentityFile ~/.ssh/cat_deploy_ed25519
    IdentitiesOnly yes
```

On the server, with your own login, from the checkout once it contains these files:

```bash
sudo install -o root -g root -m 755 hosting/tencent/cat-deploy.sh /usr/local/bin/cat-deploy
sudo install -o root -g root -m 644 hosting/tencent/cat-deploy.conf.example /etc/cat-deploy.conf
sudo nano /etc/cat-deploy.conf          # COMPOSE_PROJECT from `docker compose ls`; CONTACT optional
sudo -n true || { sudo install -o root -g root -m 440 hosting/tencent/sudoers-cat-deploy /etc/sudoers.d/cat-deploy && sudo visudo -c; }
echo 'command="/usr/local/bin/cat-deploy",restrict ssh-ed25519 AAAA...PUBLIC-KEY... claude-deploy' >> ~/.ssh/authorized_keys
```

The script is installed by copy, so a `git pull` can't change what the key may do: reinstall it by hand after reviewing a change.
To revoke the key, delete its `authorized_keys` line. Check it from your machine with `ssh cat-deploy status`. `ssh cat-deploy bash` must answer `refused`.

**Adding an admin's SSH key.** The Lighthouse key pair bound to the instance (`cat_platform`) is the only one `ubuntu` accepts at
first. To let another admin in without a restart, log in through the console (instance → **Log in**, the browser terminal) and
append their public key:

```bash
echo 'ssh-ed25519 AAAA... their-comment' >> /home/ubuntu/.ssh/authorized_keys
```

Importing the key in the console and binding it to the instance works too, but may restart the instance. fail2ban bans an address for
an hour after 5 failed logins in 10 minutes, so check the key's fingerprint (`ssh-keygen -l -f`) before trying it. A private key saved
on Windows may need LF line endings and a final newline before OpenSSH accepts it (`Load key ...: invalid format`).

**Still to do:**
- nightly `pg_dump -Fc` with copies kept off the server
- a restore test
- a Lighthouse disk and traffic alert
- an external uptime check on `/health`
