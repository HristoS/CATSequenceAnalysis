# Overnight catalogue run (2026-11-22): fills what the local catalogue lacks and rebuilds what depends on it. Windows PowerShell 5.1.
# Run from the repository root:  powershell -ExecutionPolicy Bypass -File scripts\platform-checks\overnight_catalogue.ps1
# Every step resumes where it stopped (each genome commits alone), so a step cut short is picked up by running the script again.
# Steps run in order and the run goes on after a failed step; the log and a summary go to .logs\overnight-<time>\ (git-ignored).
# The connection string comes from SequencePlatform.Ingest\appsettings.json, as for any Ingest run. No --contact is passed, so no
# email address is sent to BV-BRC or NCBI.
#
# Order, with times measured or extrapolated on 2026-11-22 (see docs/platform/README.md):
#  1. gene names (migration 0025): NCBI ~11 min (header-only CDS streams), then BV-BRC by kind from feature metadata, no genome download
#     (viruses ~25 min, archaea and plasmids minutes, bacteria ~2.5 h; extrapolated from 400 sampled genomes);
#  2. NCBI: the 11 failed and 5 pending crawl items, registered ones only (no new draw, no new animals), then the gene layer, rechecking
#     the 456 eukaryotes NCBI hadn't annotated;
#  3. BV-BRC: one genome per species not yet stored, per group (count unknown until listed), and their gene layer (names included);
#  4. rebuilds after ingest (plan rules): origins, score null, neighbour graph, gene clustering, protein null;
#  5. per-chromosome levels for the eukaryotes without them. Last, because nothing searches them; it resumes where it stopped,
#     smallest genomes first. The first run (2026-11-22) did all 582 (254.5 Gbp) in 4.0 h over 4 connections.
param(
    [int]$NewSpeciesPerGroup = 5000,
    [int]$Parallel = 4,
    [int]$ChromosomeConnections = 4,
    [switch]$SkipChromosomes
)

$ErrorActionPreference = "Continue"
$root = (Get-Location).Path
$stamp = Get-Date -Format "yyyy-MM-dd_HHmm"
$logs = Join-Path $root ".logs\overnight-$stamp"
New-Item -ItemType Directory -Force $logs | Out-Null
$summary = Join-Path $logs "summary.txt"

function Step([string]$name, [string[]]$arguments) {
    $log = Join-Path $logs ("{0:D2}_{1}.log" -f ($script:stepNumber++), $name)
    $started = Get-Date
    "{0:HH:mm:ss}  start  {1}  ({2})" -f $started, $name, ($arguments -join " ") | Tee-Object -FilePath $summary -Append
    & dotnet run --project SequencePlatform.Ingest -c Release --no-build -- @arguments *>&1 | Tee-Object -FilePath $log
    $code = $LASTEXITCODE
    "{0:HH:mm:ss}  end    {1}  exit {2}, {3:N0} min" -f (Get-Date), $name, $code, ((Get-Date) - $started).TotalMinutes | Tee-Object -FilePath $summary -Append
}
$script:stepNumber = 1

"Overnight catalogue run, started {0:yyyy-MM-dd HH:mm}" -f (Get-Date) | Tee-Object -FilePath $summary
& dotnet build SequencePlatform.Ingest -c Release *>&1 | Out-File (Join-Path $logs "00_build.log")
if ($LASTEXITCODE -ne 0) { "Build failed; see 00_build.log." | Tee-Object -FilePath $summary -Append; exit 1 }

Step "migrate" @("migrate")
Step "stats-before" @("stats")

# 1. Gene names.
Step "names-ncbi" @("backfill-gene-names", "--provider", "ncbi")
foreach ($kind in "virus", "archaeon", "plasmid", "bacterium") { Step "names-$kind" @("backfill-gene-names", "--provider", "bvbrc", "--kind", $kind) }

# 2. NCBI: registered items only, then the gene layer with a recheck of the unannotated.
foreach ($taxon in "fungi", "plants", "animals") {
    Step "ncbi-$taxon" @("ingest-ncbi", "--taxon", $taxon, "--registered-only", "--retry-failed", "--limit", "1000", "--parallel", "2")
}
Step "genes-ncbi" @("backfill-genes-ncbi", "--recheck", "--limit", "2000")

# 3. BV-BRC: new species, then their gene layer.
foreach ($group in "viruses", "archaea", "bacteria") {
    Step "bvbrc-$group" @("ingest-bvbrc", "--group", $group, "--per-species", "--limit", "$NewSpeciesPerGroup")
}
Step "genes-bvbrc" @("backfill-genes", "--limit", "100000", "--parallel", "$Parallel")
Step "names-left" @("backfill-gene-names")

# 4. Rebuilds after an ingest.
Step "build-origins" @("build-origins")
Step "build-null" @("build-null")
Step "build-neighbours" @("build-neighbours")
Step "cluster-genes" @("cluster-genes")
Step "build-protein-null" @("build-protein-null")

# 5. Per-chromosome levels, last.
if (-not $SkipChromosomes) { Step "chromosomes" @("backfill-chromosomes", "--limit", "1000", "--parallel", "$ChromosomeConnections") }

Step "stats-after" @("stats")
"Finished {0:yyyy-MM-dd HH:mm}. Logs in {1}" -f (Get-Date), $logs | Tee-Object -FilePath $summary -Append
