---
name: deploy-platform
description: Deploy, check, or administer the hosted SequencePlatform server (Tencent Lighthouse) through the restricted `ssh cat-deploy` key. Use when the user asks to deploy or redeploy the platform, roll back, check the server's status, health or logs, run a maintenance job there (build-null, build-protein-null, cluster-genes, ...), or add, rotate or disable an API user on the server.
---

# Deploy the platform

The server is reached only as `ssh cat-deploy <action>`. That key runs `/usr/local/bin/cat-deploy`
([`hosting/tencent/cat-deploy.sh`](../../../hosting/tencent/cat-deploy.sh)) and nothing else: no shell, no scp, no port forwarding.
Run `ssh cat-deploy help` for the list of actions.

## Hard rules
- **Never look for access details.** Don't read `~/.ssh`, `docs/tencent`, `authorized_keys`, `.pem`/key files or old transcripts, and
  don't search for IP addresses or usernames. A hook blocks this, and security flagged it once already. Everything needed is the
  `cat-deploy` alias.
- **If `ssh cat-deploy` fails to connect** (`Permission denied`, `agent has no identities`, timeout), stop and tell the user.
  They need to run `ssh-add` for the deploy key, or check the alias. Don't retry in a loop: fail2ban bans the address after 5 failed
  logins in 10 minutes.
- **Anything outside the action list** (seeding or restoring the database, Caddy, sshd, firewall, OS updates, editing
  `/etc/cat-deploy.conf`) is for the user, with their own login. Point them to
  [`hosting/tencent/README.md`](../../../hosting/tencent/README.md) or [`seed-database.md`](../../../hosting/tencent/seed-database.md).
- **Pushing is the user's call.** Ask before every `git push`.

## Deploy
1. `ssh cat-deploy status`: note the running commit and whether the stack is healthy.
2. Pick the target: local `HEAD` of `feature/localPlatfromBuild`, unless the user names a commit. Check that the working tree is
   clean and the commit is on `origin/feature/localPlatfromBuild` (`git fetch`, then `git branch -r --contains <sha>`). If it isn't
   pushed, ask to push.
3. Show what changes: `git log --oneline <running>..<target>` and
   `git diff --name-only <running>..<target> -- SequencePlatform.Core/Migrations/`. For each new migration, check the README's
   "Derived tables after a deploy" and "Adding a layer" sections for a job it needs (for example 0019 → `build-protein-null`).
4. Tell the user what will happen, then run `ssh cat-deploy deploy <sha>` with a 10-minute timeout. It refuses unpushed commits and a
   dirty server checkout, publishes, restarts and waits for `/health`. There are no snapshots, and **migrations are not undone by a
   rollback**. A deploy that adds a migration can only be rolled back as code.
5. `ssh cat-deploy status` again, and `ssh cat-deploy admin stats` when a migration touched data. Report: old → new commit, health,
   and the jobs still to run.

If the deploy fails: `ssh cat-deploy logs api 200` (or `migrate`, `worker`). If the service is down and the fix isn't obvious, offer
`ssh cat-deploy rollback`, which redeploys the previous commit.

## Maintenance jobs
`ssh cat-deploy job start <name>`, where name is `build-null`, `build-neighbours`, `build-origins`, `build-protein-null`,
`load-taxonomy`, `verify-genes` or `cluster-genes`. Each runs detached in its own container. Follow it with `job status` and
`job log <name> 100`, and clear it with `job clean <name>` once it has exited.
- These rewrite tables, and nothing can undo that. Say so and get a yes first.
- `cluster-genes` holds an exclusive lock on the gene table for minutes, so protein searches wait. The script checks free disk first.
- Order after new records or origins: `build-origins` → `build-null` → `build-neighbours`.
- Poll gently: the 2 vCPU server takes several times the local run times quoted in the README.

## Users
- `ssh cat-deploy user add <name>` / `user rotate <name>`: the API key is written to `~/api-keys/<name>.key` (mode 600) on the server,
  and only the path is printed. Never try to read or print the key. The user fetches it with their own login.
- `ssh cat-deploy user disable <name>`.
- `ssh cat-deploy admin users` lists names, creation dates and status. Keep `hristo` and `strahil` active.

## Read-only checks (no confirmation needed)
`status`, `logs <api|worker|migrate|postgres> [lines]`, `admin stats`, `admin users`,
`admin list [--kind K] [--rank R] [--limit N] [--after ID]`, `job status`, `job log <name> [lines]`.
