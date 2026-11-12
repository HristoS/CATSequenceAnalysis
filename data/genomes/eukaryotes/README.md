# Eukaryote genomes for the cheap checks

Kept on this machine only (the files are git-ignored), so a check reads them from disk instead of downloading them again. Fetch them with

```
dotnet run scripts/platform-checks/fetch_eukaryotes.cs      # from the repo root, SEQUENCE_PLATFORM_DB set
```

and read them in a check through `scripts/platform-checks/shared/EukaryoteFiles.cs` (`#:include shared/EukaryoteFiles.cs`), which also
fetches anything missing. Each genome is the file the platform profiled: its SHA-256 must equal the stored `raw_sha256` (the CDS file's
the gene layer's `cds_sha256`), or it isn't kept. The platform database never holds these sequences.

Checks train on one genome of a pair and test on its relative, so the genomes come in pairs (`EukaryoteFiles.Pairs`, 2026-11-11):

| Pair | Records | Size | Why |
|---|---|---|---|
| *Podospora bellae-mahoneyi* / *P. pseudocomata* | 28677 / 31652 | ~35 Mbp each | compact fungus, genes ~1.5 kbp, short introns |
| *Nelumbo nucifera* / *N. lutea* | 61826 / 53150 | ~820 Mbp each | plant, long introns, ~5% of windows coding |

Per genome: `*_genomic.fna.gz` (NCBI genome FASTA), `*_cds_from_genomic.fna.gz` (annotated coding sequences, spliced) and
`*_assembly_report.txt` (which sequences are nuclear), all from NCBI's FTP tree.
