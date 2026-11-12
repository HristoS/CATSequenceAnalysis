#!/usr/bin/env bash
# /usr/local/bin/cat-deploy (root-owned, 755): the only program the restricted deploy key may run.
# authorized_keys line:  command="/usr/local/bin/cat-deploy",restrict ssh-ed25519 AAAA... claude-deploy
# sshd puts the client's command in SSH_ORIGINAL_COMMAND; this script splits it on whitespace, never evaluates it, checks every
# argument against a pattern, and runs one of the fixed actions below. Anything else is refused and logged (journalctl -t cat-deploy).
# Settings: /etc/cat-deploy.conf (see cat-deploy.conf.example). Installed by hand, so a git pull can't change what this key can do.
set -euo pipefail
umask 077

CONF=/etc/cat-deploy.conf
[[ -r $CONF ]] || { echo "missing $CONF" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONF"
: "${REPO:?}" "${BRANCH:?}" "${COMPOSE_PROJECT:?}"
KEY_DIR=${KEY_DIR:-$HOME/api-keys}
STATE_DIR=${STATE_DIR:-$HOME/.cat-deploy}
CONTACT=${CONTACT:-}
COMPOSE_DIR=$REPO/SequencePlatform.AppHost/deploy/compose
HEALTH_URL=http://127.0.0.1:5080/health

mkdir -p "$STATE_DIR" "$KEY_DIR"
read -r -a ARGS <<< "${SSH_ORIGINAL_COMMAND:-help}"
ACTION=${ARGS[0]:-help}
CLIENT=${SSH_CLIENT:-local}
log() { logger -t cat-deploy -- "from=${CLIENT%% *} $*"; }
refuse() { log "refused: ${SSH_ORIGINAL_COMMAND:-}"; echo "refused: $*" >&2; exit 2; }

for a in "${ARGS[@]}"; do [[ $a =~ ^[A-Za-z0-9._=-]{1,64}$ ]] || refuse "argument '$a' has characters outside [A-Za-z0-9._=-]"; done
(( ${#ARGS[@]} <= 8 )) || refuse "too many arguments"
log "action: ${ARGS[*]}"

compose() { (cd "$COMPOSE_DIR" && docker compose -p "$COMPOSE_PROJECT" -f docker-compose.yaml --env-file .env.Production "$@"); }
ingest() { compose run --rm --no-deps -T migrate "$@"; }
psql_q() { compose exec -T postgres sh -c 'PGPASSWORD="$POSTGRES_PASSWORD" psql -U postgres -d sequence_platform -v ON_ERROR_STOP=1 -At -c "$0"' "$1"; }

lock() { exec 9> "$STATE_DIR/lock"; flock -n 9 || refuse "another deploy or user action is running"; }

wait_healthy() {
  for _ in $(seq 1 60); do
    if curl -fsS --max-time 3 "$HEALTH_URL" 2>/dev/null; then echo; return 0; fi
    sleep 2
  done
  echo "NOT HEALTHY after 120 s; see: logs api" >&2; return 1
}

usage() {
  cat <<'EOF'
cat-deploy actions:
  status                          commit, containers, /health, disk and memory
  deploy <commit>                 check out a pushed commit of the branch, publish the Compose stack, restart, check /health
  rollback                        deploy the commit that ran before the last deploy (migrations are not undone)
  logs <api|worker|migrate|postgres> [lines]
  admin stats | users | list [--kind K] [--rank R] [--limit N] [--after ID]
  job start <build-null|build-neighbours|build-origins|build-protein-null|load-taxonomy|verify-genes|cluster-genes>
  job status | job log <name> [lines] | job clean <name>
  user add <name> | user rotate <name> | user disable <name>
EOF
}

# --- actions -----------------------------------------------------------------------------------------------------------------------

do_status() {
  echo "commit:   $(git -C "$REPO" log -1 --format='%h %cs %s')"
  echo "branch:   $(git -C "$REPO" rev-parse --abbrev-ref HEAD)"
  echo "previous: $(cat "$STATE_DIR/previous" 2>/dev/null || echo none)"
  echo "service:  $(systemctl is-active cat-platform || true)"
  echo -n "health:   "; curl -fsS --max-time 5 "$HEALTH_URL" || echo "unreachable"; echo
  compose ps --format 'table {{.Service}}\t{{.State}}\t{{.Status}}'
  df -h / | tail -1 | awk '{print "disk:     " $4 " free of " $2 " (" $5 " used)"}'
  free -m | awk '/^Mem/ {print "memory:   " $7 " MB available of " $2 " MB"}'
}

do_deploy() {
  local target=$1
  [[ $target =~ ^[0-9a-f]{7,40}$ ]] || refuse "deploy needs a commit hash (7-40 hex characters)"
  lock
  cd "$REPO"
  [[ -z $(git status --porcelain --untracked-files=no) ]] || refuse "the server checkout has local changes; fix by hand"
  git fetch -q origin "$BRANCH"
  local sha; sha=$(git rev-parse --verify -q "$target^{commit}") || refuse "commit $target not found after fetch"
  git merge-base --is-ancestor "$sha" "origin/$BRANCH" || refuse "commit $sha is not on origin/$BRANCH"
  local current; current=$(git rev-parse HEAD)
  [[ $sha != "$current" ]] || echo "already at $sha; publishing again"
  [[ $sha == "$current" ]] || echo "$current" > "$STATE_DIR/previous"
  git checkout -q -B "$BRANCH" "$sha"
  echo "== checked out $(git log -1 --format='%h %s')"
  echo "== publishing (about 2 minutes on 2 vCPU)"
  dotnet run --project SequencePlatform.AppHost --no-launch-profile -- --operation publish --step deploy --output-path deploy/compose
  chmod 600 "$COMPOSE_DIR"/.env* 2>/dev/null || true
  chmod 600 "$HOME"/.aspire/deployments/*/production.json* 2>/dev/null || true
  echo "== restarting cat-platform"
  sudo -n systemctl restart cat-platform
  wait_healthy
  log "deployed $sha (previous $current)"
  echo "== deployed $(git log -1 --format='%h %cs %s')"
}

do_logs() {
  local service=${1:-} lines=${2:-200}
  [[ $service =~ ^(api|worker|migrate|postgres)$ ]] || refuse "logs: service must be api, worker, migrate or postgres"
  [[ $lines =~ ^[0-9]{1,4}$ ]] && (( lines <= 2000 )) || refuse "logs: lines must be 1-2000"
  compose logs --no-color --tail "$lines" "$service"
}

do_admin() {
  case ${1:-} in
    stats) ingest stats ;;
    users) psql_q "select name, created_at::date, coalesce(disabled_at::date::text, 'active') from app_user order by id" ;;
    list)
      shift
      local opts=()
      while (( $# )); do
        [[ $1 =~ ^--(kind|rank|limit|after|min-length|max-length|provider)$ && $# -ge 2 ]] || refuse "list: allowed options are --kind --rank --limit --after --min-length --max-length --provider"
        opts+=("$1" "$2"); shift 2
      done
      ingest list "${opts[@]}" ;;
    *) refuse "admin: stats, users or list" ;;
  esac
}

job_name() { [[ ${1:-} =~ ^(build-null|build-neighbours|build-origins|build-protein-null|load-taxonomy|verify-genes|cluster-genes)$ ]] || refuse "job: unknown job '${1:-}'"; echo "cat-job-$1"; }

do_job() {
  case ${1:-} in
    start)
      local job=${2:-} name; name=$(job_name "$job")
      lock
      if docker ps -q --filter "name=^${name}$" | grep -q .; then refuse "$name is already running"; fi
      docker rm "$name" >/dev/null 2>&1 || true
      local extra=()
      if [[ $job == cluster-genes ]]; then
        local need free
        need=$(psql_q "select pg_total_relation_size('gene_profile')")
        free=$(( $(df --output=avail -B1 / | tail -1) ))
        echo "gene_profile: $(( need / 1024 / 1024 )) MB, free disk: $(( free / 1024 / 1024 )) MB"
        (( free > need * 12 / 10 )) || refuse "cluster-genes needs free disk of about the table's size plus 20%"
      fi
      if [[ $job =~ ^(build-protein-null|load-taxonomy)$ && -n $CONTACT ]]; then extra=(--contact "$CONTACT"); fi
      compose run -d --name "$name" --no-deps -T migrate "$job" "${extra[@]}"
      echo "started $name; follow with: job log $job" ;;
    status) docker ps -a --filter 'name=^cat-job-' --format 'table {{.Names}}\t{{.Status}}\t{{.RunningFor}}' ;;
    log)
      local name lines=${3:-100}; name=$(job_name "${2:-}")
      [[ $lines =~ ^[0-9]{1,4}$ ]] || refuse "job log: lines must be a number"
      docker logs --tail "$lines" "$name" 2>&1 ;;
    clean)
      local name; name=$(job_name "${2:-}")
      if docker ps -q --filter "name=^${name}$" | grep -q .; then refuse "$name is still running"; fi
      docker rm "$name" ;;
    *) refuse "job: start, status, log or clean" ;;
  esac
}

do_user() {
  local verb=${1:-} name=${2:-}
  [[ $name =~ ^[a-z][a-z0-9_-]{1,31}$ ]] || refuse "user: name must be lowercase letters, digits, - or _ (2-32)"
  lock
  case $verb in
    add|rotate)
      local cmd=add-user out key file="$KEY_DIR/$name.key"
      [[ $verb == rotate ]] && cmd=rotate-key
      # The Ingest prints the key once; it goes straight to a 600 file and never to this session's output.
      out=$(ingest "$cmd" "$name" 2>/dev/null) || { echo "$cmd failed for $name (output withheld: it may contain a key)"; exit 1; }
      key=$(grep -oE 'shown only now: sp_[A-Za-z0-9_-]+' <<< "$out" | head -1 | sed 's/^shown only now: //') || true
      [[ -n $key ]] || { echo "$cmd ran but no key was found in its output (output withheld)"; exit 1; }
      printf '%s\n' "$key" > "$file"; chmod 600 "$file"
      log "user $verb $name"
      echo "user $name: API key written to $file on the server (mode 600). Fetch it with your own login, not through this key." ;;
    disable) ingest disable-user "$name"; log "user disable $name"; echo "user $name disabled" ;;
    *) refuse "user: add, rotate or disable" ;;
  esac
}

case $ACTION in
  help) usage ;;
  status) do_status ;;
  deploy) do_deploy "${ARGS[1]:-}" ;;
  rollback) [[ -s $STATE_DIR/previous ]] || refuse "no previous commit recorded"; do_deploy "$(cat "$STATE_DIR/previous")" ;;
  logs) do_logs "${ARGS[1]:-}" "${ARGS[2]:-200}" ;;
  admin) do_admin "${ARGS[@]:1}" ;;
  job) do_job "${ARGS[@]:1}" ;;
  user) do_user "${ARGS[1]:-}" "${ARGS[2]:-}" ;;
  *) refuse "unknown action '$ACTION' (try: help)" ;;
esac
