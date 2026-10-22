import glob, os, re, sys

root = sys.argv[1]
os.chdir(root)

cols = []
for f in sorted(glob.glob("SequencePlatform.Core/Migrations/*.sql")):
    table = None
    for line in open(f, encoding="utf-8"):
        m = re.match(r"\s*create table (\w+)", line, re.I)
        if m:
            table = m.group(1)
            continue
        m = re.match(r"\s*alter table (\w+)", line, re.I)
        if m:
            table = m.group(1)
        m = re.match(r"\s*(?:add column\s+)?([a-z_][a-z0-9_]*)\s+(bigint|integer|smallint|text|bytea|real|double|timestamptz|boolean|jsonb|date|numeric|char|varchar|float)", line, re.I)
        if m and table and m.group(1).lower() not in ("primary", "unique", "check", "constraint", "create", "select", "from", "where"):
            cols.append((os.path.basename(f), table, m.group(1)))

code = ""
for f in glob.glob("SequencePlatform.*/**/*.cs", recursive=True):
    parts = f.replace("\\", "/").split("/")
    if "bin" in parts or "obj" in parts:
        continue
    code += open(f, encoding="utf-8", errors="replace").read()

for mig, table, col in cols:
    refs = len(re.findall(r"\b" + col + r"\b", code))
    if refs == 0:
        print(f"{mig:28} {table:18} {col:22} no code reference")
print(len(cols), "columns checked")
