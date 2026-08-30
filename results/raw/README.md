# Raw external-tool output (Mash, FastANI)

Not produced by `dotnet test` — these are the raw outputs of two external alignment-free baseline
tools, run once against the same 30-genome corpus in `data/genomes/` used throughout
`RealGenomeDistanceMatrixBenchmark.cs`. Committed as source data (like the genome `.fna.gz` files
themselves) so `ExternalBaselineComparisonBenchmark.cs` can parse/analyze them via `dotnet test`
without requiring Docker just to reproduce the analysis CSVs under `results/`. Regenerating the raw
files themselves does require Docker (Docker Desktop, tested with server 27.2.0) — exact commands
below reproduce them.

## mash_dist_raw.tsv

Tool: [Mash](https://github.com/marbl/Mash) 2.3 (Ondov et al., 2016), via
`staphb/mash:latest` (digest `sha256:44005a48dae918be4d338f73a8bec82b1a4d00327d8f787673917eda80afeae7`
at time of use — pull that digest explicitly to reproduce exactly, since `:latest` can move).

```
docker run --rm -v <repo>/data/genomes:/genomes:ro -v <out>:/out staphb/mash:latest sh -c '
  cd /genomes
  FILES=$(cat /out/genome_list.txt)          # the 30 bacterial *.fna.gz filenames, one per line,
                                               # excluding the 2 phage genomes (phix174/lambda_phage)
  mash sketch -k 21 -s 1000 -o /out/sketch $FILES
  mash dist /out/sketch.msh /out/sketch.msh > /out/mash_dist_raw.tsv
'
```

900 rows (30×30, all ordered pairs including self-pairs — `ExternalBaselineComparisonBenchmark.cs`
reduces this to the 435 unique unordered pairs). Columns: reference, query, Mash distance, p-value,
shared-hashes (`x/1000`). Wall-clock for the whole sketch+dist pass on this corpus: **6 s**.

## fastani_raw.tsv

Tool: [FastANI](https://github.com/ParBLiSS/FastANI) 1.33 (Jain et al., 2018), via
`staphb/fastani:latest` (digest `sha256:f3f1303e149f63908aeb4671f0982cb53c9e376b5d3852c7194f0f0fc8bea0a8`
at time of use).

```
docker run --rm -v <repo>/data/genomes:/genomes:ro -v <out>:/out staphb/fastani:latest \
  fastANI --ql /out/genome_paths.txt --rl /out/genome_paths.txt -o /out/fastani_raw.tsv
```

`genome_paths.txt` lists the same 30 files as `/genomes/<filename>`, one per line, both as query and
reference list (all-vs-all). Columns: query path, reference path, ANI (%), count of bidirectional
fragment mappings, total query fragments. **Only 58 of the possible 900 ordered combinations are
present** — FastANI does not emit a row at all when identity falls below its effective mapping
threshold (in practice, well below the corpus's many different-phylum pairs), rather than reporting a
low value. This coverage gap is itself part of what `ExternalBaselineComparisonBenchmark.cs` reports,
not a defect in how it was run. Wall-clock for the whole all-vs-all pass on this corpus: **93 s**
(vs. Mash's 6 s on the identical input).
