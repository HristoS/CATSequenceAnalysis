# SequenceAnalise.Benchmarks

The first command-line tool for the controlled Phase 1 localization dataset. It creates and validates a deterministic recipe for the 30 bacterial references documented in `data/genomes/README.md`. Query DNA is reconstructed from the source FASTA when needed; it is not copied into the manifest.

From the repository root:

```powershell
dotnet run --project SequenceAnalise.Benchmarks -- dataset generate --genomes-dir data/genomes --output data/localization/phase1-exact-v1.json
dotnet run --project SequenceAnalise.Benchmarks -- dataset validate --genomes-dir data/genomes --manifest data/localization/phase1-exact-v1.json
```

`generate` accepts an optional `--seed INT` (default `20261001`) and refuses to overwrite an existing output file. The manifest records each compressed FASTA SHA-256, selected primary record, deterministic sample seed, coordinates, query lengths, and hashes of the reconstructed query sequences. The validator checks the reference bytes, primary record, protocol metadata, deterministic locus coordinates, and query hashes.

The fixed protocol uses 20 non-overlapping, 500 bp-aligned loci per reference; primary query lengths are 2,000, 5,000, and 10,000 bp, with 500 bp and 1,000 bp stress cases. The 30-reference selection is explicit, so the two phage FASTAs at the root of `data/genomes/` and the nested viral corpus are not accidentally added by directory globbing.