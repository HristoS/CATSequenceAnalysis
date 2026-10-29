"""S4.3: time searches through the running API and its worker on the local catalogue (results/S43_WorkerTiming.csv).

POSTs each query, polls GET until Done, and prints the POST time and POST-to-Done time. Run from the repo root with the
API running (Release) and a key from the Ingest's add-user in a file:  python -I worker_smoke.py http://localhost:5080 key.txt data/genomes
The worker's own time per search is in the API's log ("Search ...: done in N ms"). Delete the smoke user's searches afterwards."""
import gzip, json, sys, time, urllib.request

base, key_path, genomes = sys.argv[1], sys.argv[2], sys.argv[3]
key = open(key_path).read().strip()
headers = {"Authorization": f"Bearer {key}"}


def request(method, path, body=None, content_type=None):
    h = dict(headers)
    if content_type:
        h["Content-Type"] = content_type
    req = urllib.request.Request(base + path, data=body, headers=h, method=method)
    with urllib.request.urlopen(req) as r:
        return r.status, json.loads(r.read() or b"null")


for _ in range(60):
    try:
        urllib.request.urlopen(base + "/openapi/v1.json")
        break
    except Exception:
        time.sleep(1)


def fasta(path, start=None, length=None):
    lines = gzip.open(path, "rt").read().split("\n")
    seq = "".join(l.strip() for l in lines[1:] if l and not l.startswith(">"))
    if start is not None:
        seq = seq[start:start + length]
    return f">q\n{seq}\n".encode()


def run(label, route, body):
    t0 = time.perf_counter()
    status, accepted = request("POST", route, body, "text/x-fasta")
    t_post = time.perf_counter() - t0
    while True:
        _, run = request("GET", accepted["href"] + "?limit=3")
        if run["status"] in ("Done", "Failed"):
            break
        time.sleep(0.02)
    total = time.perf_counter() - t0
    hits = (run.get("hits") or {}).get("items", [])
    top = ", ".join(f"{h['sequence']['name'][:40]} ({h.get('score', 0):.3f}{', off ' + str(h['offsetBp']) if 'offsetBp' in h else ''})" for h in hits[:2])
    print(f"{label:44} POST {t_post*1000:6.0f} ms  POST->{run['status']} {total*1000:6.0f} ms  | {top}")


sars = f"{genomes}/viral/sars_cov_2_NC_045512.2.fna.gz"
k12 = f"{genomes}/ecoli_k12_mg1655_GCF_000005845.2.fna.gz"
run("similar SARS-CoV-2 (cold)", "/api/searches/similar?top=10", fasta(sars))
run("similar SARS-CoV-2 (warm)", "/api/searches/similar?top=10", fasta(sars))
run("similar E. coli K-12 (cold band)", "/api/searches/similar?top=10", fasta(k12))
run("similar E. coli K-12 (warm)", "/api/searches/similar?top=10", fasta(k12))
run("localize K-12 10 kbp @1,000,000 (streamed)", "/api/searches/localize?kind=bacterium&top=5", fasta(k12, 1_000_000, 10_000))
run("localize K-12 50 kbp @2,500,000 (streamed)", "/api/searches/localize?kind=bacterium&top=5", fasta(k12, 2_500_000, 50_000))
