#:project ../../SequencePlatform.Experiments/SequencePlatform.Experiments.csproj
using Npgsql;
using SequencePlatform.Core;
using SequencePlatform.Experiments;

// S5.1 gate (P5, R2): does today's distance, 1 - Compare, separate taxa beyond GC content on pairs whose GC is nearly equal? The
// experiment itself lives in SequencePlatform.Experiments/GcStratum.cs since 2026-11-21, where the experiments page reruns it; this script
// runs the same code. Moving it changed no byte of the CSV on the local catalogue that day.
// Args: output path (default S51_GcStratumCheck.csv), window size (default 32000).
var output = args.Length > 0 ? args[0] : GcStratum.FileName;
int windowSize = args.Length > 1 ? int.Parse(args[1]) : GcStratum.DefaultWindowSize;
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

var result = await new GcStratum(windowSize).RunAsync(dataSource, ExperimentContext.Default with { Progress = Console.Error.WriteLine },
    CancellationToken.None);
var csv = result.Files.Single().Content;
await File.WriteAllBytesAsync(output, csv);
Console.Write(System.Text.Encoding.UTF8.GetString(csv));
