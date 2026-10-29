"""The API's end-to-end scenario (P4 S4.5), the one SequencePlatform.Api/Scenario.http walks through by hand.

  python -I scripts/platform-checks/e2e_scenario.py prepare                  writes the scenario's FASTA files
  python -I scripts/platform-checks/e2e_scenario.py run BASE_URL KEY_FILE    sends its requests and checks the answers

prepare writes plain FASTA from data/genomes/ into SequencePlatform.Api/scenario/ (git-ignored): the 29 bundled viral genomes,
SARS-CoV-2 alone as a query, its S gene, and E. coli K-12 MG1655 with that S gene planted at 2,000,000. run expects a fresh
database (migrated, one user), the API on BASE_URL, and the user's key in KEY_FILE; it prints each step and stops at the first
answer that isn't the expected one. Run from the repository root.
"""
import glob, gzip, json, os, sys, time, urllib.error, urllib.request

ROOT = os.getcwd()
OUT = os.path.join(ROOT, "SequencePlatform.Api", "scenario")
PLANT_AT = 2_000_000
# NC_045512.2 S gene, 21563..25384 (1-based, inclusive).
SPIKE = (21_562, 25_384)


def read_fasta(path):
    with gzip.open(path, "rt") as f:
        lines = f.read().split("\n")
    return lines[0][1:], "".join(l.strip() for l in lines[1:] if l and not l.startswith(">"))


def write_fasta(name, records):
    with open(os.path.join(OUT, name), "w", newline="\n") as f:
        for header, seq in records:
            f.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")


def prepare():
    os.makedirs(OUT, exist_ok=True)
    viruses = [read_fasta(p) for p in sorted(glob.glob(os.path.join(ROOT, "data", "genomes", "viral", "*.fna.gz")))]
    write_fasta("viruses.fna", viruses)
    sars = next(v for v in viruses if v[0].startswith("NC_045512"))
    write_fasta("sars_cov_2.fna", [sars])
    spike = sars[1][SPIKE[0]:SPIKE[1]]
    write_fasta("spike.fna", [("S SARS-CoV-2 surface glycoprotein gene", spike)])
    _, k12 = read_fasta(os.path.join(ROOT, "data", "genomes", "ecoli_k12_mg1655_GCF_000005845.2.fna.gz"))
    write_fasta("k12_planted.fna", [(f"k12_planted E. coli K-12 MG1655 with the SARS-CoV-2 S gene at {PLANT_AT:,}",
                                     k12[:PLANT_AT] + spike + k12[PLANT_AT:])])
    for name in sorted(os.listdir(OUT)):
        print(f"{name:20} {os.path.getsize(os.path.join(OUT, name)):>10,} bytes")


def run(base, key_file):
    key = open(key_file).read().strip()

    def request(method, path, body=None, content_type=None):
        headers = {"Authorization": f"Bearer {key}"}
        if content_type:
            headers["Content-Type"] = content_type
        req = urllib.request.Request(base + path, data=body, headers=headers, method=method)
        try:
            with urllib.request.urlopen(req) as r:
                return r.status, r.headers, json.loads(r.read() or b"null")
        except urllib.error.HTTPError as e:
            return e.code, e.headers, json.loads(e.read() or b"null")

    def fasta(name):
        return open(os.path.join(OUT, name), "rb").read()

    def check(step, condition, detail):
        print(f"{'ok  ' if condition else 'FAIL'} {step}: {detail}")
        if not condition:
            sys.exit(1)

    def poll(href):
        t0 = time.perf_counter()
        while True:
            _, _, run = request("GET", href)
            if run["status"] in ("Done", "Failed"):
                return run, time.perf_counter() - t0
            time.sleep(0.1)

    viruses = fasta("viruses.fna").count(b">")
    status, headers, body = request("POST", "/api/sequences?kind=virus", fasta("viruses.fna"), "text/x-fasta")
    check("1 upload viruses", status == 201 and len(body["records"]) == viruses, f"{status}, {len(body['records'])} records, Location {headers['Location']}")
    sars_id = next(r["id"] for r in body["records"] if r["recordId"].startswith("NC_045512"))

    status, _, body = request("POST", "/api/sequences?kind=bacterium&name=E.%20coli%20K-12%20with%20a%20planted%20S%20gene",
                              fasta("k12_planted.fna"), "text/x-fasta")
    planted = body["records"][0]
    check("2 upload planted K-12", status == 201, f"{status}, id {planted['id']}, {planted['lengthBp']:,} bp, {planted['windowCount']:,} windows")

    status, _, body = request("GET", "/api/sequences?provider=user%20upload&limit=100")
    check("3 list my uploads", status == 200 and len(body["items"]) == viruses + 1, f"{len(body['items'])} records")

    status, _, body = request("GET", f"/api/sequences/{planted['id']}")
    check("4 detail", status == 200 and body["owned"] and body["dataClass"] == "restricted",
          f"owned, {body['dataClass']}, levels {[l['windowSize'] for l in body['levels']]}")

    status, headers, accepted = request("POST", "/api/searches/similar?kind=virus&top=10", fasta("sars_cov_2.fna"), "text/x-fasta")
    check("5 find similar to SARS-CoV-2", status == 202 and headers["Location"] == accepted["href"], f"{status}, {accepted['href']}")
    run, seconds = poll(accepted["href"])
    hits = run["hits"]["items"]
    check("6 poll", run["status"] == "Done", f"Done after {seconds:.2f} s; band {run['summary']['bandMinBp']:,}-{run['summary']['bandMaxBp']:,} bp, "
          f"{run['summary']['rankedTotal']} ranked")
    print("       " + "; ".join(f"{h['position']}. {h['sequence']['name'][:45]} ({h['score']:.3f})" for h in hits[:5]))
    check("7 read the ranked list", hits[0]["sequenceId"] == sars_id and hits[0]["score"] > 0.999, "SARS-CoV-2 itself first")
    _, _, identical = request("GET", accepted["href"] + "?list=Identical")
    check("8 read the identical list", [h["sequenceId"] for h in identical["hits"]["items"]] == [sars_id], "SARS-CoV-2, by SHA-256")

    status, _, located = request("POST", f"/api/searches/localize?record={planted['id']}", fasta("spike.fna"), "text/x-fasta")
    run, seconds = poll(located["href"])
    top = run["hits"]["items"][0]
    check("9 localize the planted S gene", run["status"] == "Done" and top["sequenceId"] == planted["id"] and abs(top["offsetBp"] - PLANT_AT) <= 500
          and top["strand"] == "Forward", f"offset {top['offsetBp']:,} on the {top['strand']} strand, score {top['score']:.3f}, after {seconds:.2f} s")

    status, _, _ = request("DELETE", accepted["href"])
    _, _, expired = request("GET", accepted["href"])
    check("10 delete the search", status == 204 and expired["status"] == "Expired" and "hits" not in expired,
          f"{status}, now {expired['status']}, summary kept")


if __name__ == "__main__":
    if sys.argv[1:2] == ["prepare"]:
        prepare()
    elif sys.argv[1:2] == ["run"] and len(sys.argv) == 4:
        run(sys.argv[2], sys.argv[3])
    else:
        sys.exit(__doc__)
