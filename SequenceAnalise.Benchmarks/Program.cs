using SequenceAnalise.Benchmarks;

const string usage =
    "usage: SequenceAnalise.Benchmarks dataset generate --genomes-dir PATH --output PATH [--seed INT] | " +
    "dataset validate --genomes-dir PATH --manifest PATH | " +
    "benchmark window-size --genomes-dir PATH --manifest PATH --output PATH [--window-size INT] [--query-length INT] | " +
    "benchmark window-size-sweep --genomes-dir PATH --manifest PATH --output PATH [--window-size INT] [--query-length INT] | " +
    "benchmark divergence --genomes-dir PATH --manifest PATH --output PATH [--window-size INT] [--query-length INT] [--divergence DOUBLE] | " +
    "benchmark summarize --input PATH --output PATH --kind divergence|window-size";

if (args.Length < 2 || (args[0] != "dataset" && args[0] != "benchmark"))
{
    Console.WriteLine(usage);
    return 1;
}

string? Option(string name)
{
    var index = Array.IndexOf(args, name);
    if (index < 0) return null;
    if (index + 1 >= args.Length || args[index + 1].StartsWith("--", StringComparison.Ordinal))
        throw new ArgumentException($"Option {name} requires a value.");
    return args[index + 1];
}

string[]? MultiOption(string name)
{
    var values = new List<string>();
    for (int i = 0; i < args.Length - 1; i++)
    {
        if (args[i] == name)
            values.Add(args[i + 1]);
    }
    return values.Count == 0 ? null : values.ToArray();
}

try
{
    switch (args[0])
    {
        case "dataset":
            switch (args[1])
            {
                case "generate":
                {
                    var datasetGenomesDirectory = Option("--genomes-dir");
                    var datasetOutputPath = Option("--output");
                    if (datasetGenomesDirectory is null || datasetOutputPath is null)
                    {
                        Console.WriteLine(usage);
                        return 1;
                    }

                    var seedText = Option("--seed");
                    int seed = PhaseOneDataset.defaultSeed;
                    if (seedText is not null && !int.TryParse(seedText, out seed))
                        throw new ArgumentException($"Invalid integer seed '{seedText}'.");
                    var manifest = PhaseOneDataset.Generate(datasetGenomesDirectory, datasetOutputPath, seed);
                    int queryCount = manifest.References.Sum(r => r.Loci.Sum(l => l.Queries.Length));
                    Console.WriteLine($"Wrote {manifest.References.Length} references and {queryCount} query recipes to {Path.GetFullPath(datasetOutputPath)}.");
                    return 0;
                }

                case "validate":
                {
                    var datasetGenomesDirectory = Option("--genomes-dir");
                    var datasetManifestPath = Option("--manifest");
                    if (datasetGenomesDirectory is null || datasetManifestPath is null)
                    {
                        Console.WriteLine(usage);
                        return 1;
                    }

                    var result = PhaseOneDataset.Validate(datasetGenomesDirectory, datasetManifestPath);
                    Console.WriteLine($"Valid: {result.References} references, {result.Queries} query recipes; no DNA sequences are stored in the manifest.");
                    return 0;
                }

                default:
                    Console.WriteLine(usage);
                    return 1;
            }

        case "benchmark":
            switch (args[1])
            {
                case "window-size":
                {
                    var benchmarkGenomesDirectory = Option("--genomes-dir");
                    var benchmarkManifestPath = Option("--manifest");
                    var benchmarkOutputPath = Option("--output");
                    if (benchmarkGenomesDirectory is null || benchmarkManifestPath is null || benchmarkOutputPath is null)
                    {
                        Console.WriteLine(usage);
                        return 1;
                    }

                    var windowSizesText = MultiOption("--window-size");
                    var queryLengthsText = MultiOption("--query-length");
                    var benchmarkWindowSizes = windowSizesText is null ? new[] { 500 } : windowSizesText.Select(int.Parse).ToArray();
                    var benchmarkQueryLengths = queryLengthsText is null ? new[] { 2_000, 5_000, 10_000 } : queryLengthsText.Select(int.Parse).ToArray();
                    WindowSizeLocalizationBenchmark.Run(benchmarkGenomesDirectory, benchmarkManifestPath, benchmarkOutputPath, benchmarkWindowSizes, benchmarkQueryLengths);
                    Console.WriteLine($"Wrote window-size benchmark rows to {Path.GetFullPath(benchmarkOutputPath)}.");
                    return 0;
                }

                case "window-size-sweep":
                {
                    var sweepGenomesDirectory = Option("--genomes-dir");
                    var sweepManifestPath = Option("--manifest");
                    var sweepOutputPath = Option("--output");
                    if (sweepGenomesDirectory is null || sweepManifestPath is null || sweepOutputPath is null)
                    {
                        Console.WriteLine(usage);
                        return 1;
                    }

                    var sweepWindowSizesText = MultiOption("--window-size");
                    var sweepQueryLengthsText = MultiOption("--query-length");
                    var sweepWindowSizes = sweepWindowSizesText is null ? new[] { 250, 500, 1000, 2000 } : sweepWindowSizesText.Select(int.Parse).ToArray();
                    var sweepQueryLengths = sweepQueryLengthsText is null ? new[] { 2_000, 5_000, 10_000 } : sweepQueryLengthsText.Select(int.Parse).ToArray();
                    WindowSizeRobustnessBenchmark.Run(sweepGenomesDirectory, sweepManifestPath, sweepOutputPath, sweepWindowSizes, sweepQueryLengths);
                    Console.WriteLine($"Wrote window-size robustness rows to {Path.GetFullPath(sweepOutputPath)}.");
                    return 0;
                }

                case "divergence":
                {
                    var divergenceGenomesDirectory = Option("--genomes-dir");
                    var divergenceManifestPath = Option("--manifest");
                    var divergenceOutputPath = Option("--output");
                    if (divergenceGenomesDirectory is null || divergenceManifestPath is null || divergenceOutputPath is null)
                    {
                        Console.WriteLine(usage);
                        return 1;
                    }

                    var divergenceWindowSizeText = Option("--window-size");
                    var divergenceQueryLengthsText = MultiOption("--query-length");
                    var divergenceRateText = MultiOption("--divergence");
                    var divergenceWindowSize = divergenceWindowSizeText is null ? 500 : int.Parse(divergenceWindowSizeText);
                    var divergenceQueryLengths = divergenceQueryLengthsText is null ? new[] { 2_000, 5_000, 10_000 } : divergenceQueryLengthsText.Select(int.Parse).ToArray();
                    var divergenceRates = divergenceRateText is null ? new[] { 0.00, 0.01, 0.03, 0.05, 0.10, 0.20, 0.30, 0.50 } : divergenceRateText.Select(double.Parse).ToArray();
                    DivergenceLocalizationBenchmark.Run(divergenceGenomesDirectory, divergenceManifestPath, divergenceOutputPath, divergenceWindowSize, divergenceQueryLengths, divergenceRates);
                    Console.WriteLine($"Wrote divergence benchmark rows to {Path.GetFullPath(divergenceOutputPath)}.");
                    return 0;
                }

                case "summarize":
                {
                    var inputPath = Option("--input");
                    var outputPath = Option("--output");
                    var kind = Option("--kind");
                    if (inputPath is null || outputPath is null || kind is null)
                    {
                        Console.WriteLine(usage);
                        return 1;
                    }

                    switch (kind)
                    {
                        case "divergence":
                            BenchmarkSummary.WriteDivergenceSummary(inputPath, outputPath);
                            break;
                        case "window-size":
                            BenchmarkSummary.WriteWindowSizeSummary(inputPath, outputPath);
                            break;
                        default:
                            Console.WriteLine(usage);
                            return 1;
                    }

                    Console.WriteLine($"Wrote benchmark summary to {Path.GetFullPath(outputPath)}.");
                    return 0;
                }

                default:
                    Console.WriteLine(usage);
                    return 1;
            }

        default:
            Console.WriteLine(usage);
            return 1;
    }
}
catch (Exception exception)
{
    Console.Error.WriteLine($"Error: {exception.Message}");
    return 2;
}