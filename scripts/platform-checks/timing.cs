#:project ../../SequencePlatform.Core/SequencePlatform.Core.csproj
using System.Diagnostics;
using System.Security.Cryptography;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

// S3.7: one query of each kind on the local catalogue, as the engine stands at the end of P3 (length band, neighbour re-ranking,
// origin rotation, localization). Per query: profiling the DNA, a first search on a fresh cached engine (it loads what it needs),
// the median of 5 warm searches on that engine, and the median of 3 searches streamed from the database on fresh uncached
// engines (as the CLI runs). Queries come from the bundled genomes; O157:H7 Sakai is stored in the catalogue, K-12 isn't.
// Pass an output path (default S37_Timing.csv).
var output = args.Length > 0 ? args[0] : "S37_Timing.csv";
await using var dataSource = NpgsqlDataSource.Create(ConnectionString.Resolve());

string Primary(string file) => FastaReader.ReadPrimaryRecord(Path.Combine("data/genomes", file)).Sequence;
var k12 = Primary("ecoli_k12_mg1655_GCF_000005845.2.fna.gz");
var o157 = Primary("ecoli_o157h7_sakai_GCF_000008865.2.fna.gz");
var sars = Primary("viral/sars_cov_2_NC_045512.2.fna.gz");

long o157Id;
await using (var command = dataSource.CreateCommand("select id from sequence_record where sequence_sha256 = $1"))
{
    command.Parameters.AddWithValue(SHA256.HashData(o157.Where(char.IsAsciiLetter).Select(c => (byte)char.ToUpperInvariant(c)).ToArray()));
    o157Id = (long)(await command.ExecuteScalarAsync())!;
}

double Median(List<double> values) => values.Order().ElementAt(values.Count / 2);

async Task<(double First, double Warm, double Streamed, string Top)> Time(Func<bool, Func<Task<string>>> engine)
{
    var cached = engine(true);
    var watch = Stopwatch.StartNew();
    var top = await cached();
    double first = watch.Elapsed.TotalMilliseconds;
    var warm = new List<double>();
    for (int i = 0; i < 5; i++)
    {
        watch.Restart();
        await cached();
        warm.Add(watch.Elapsed.TotalMilliseconds);
    }
    var streamed = new List<double>();
    for (int i = 0; i < 3; i++)
    {
        var fresh = engine(false);
        watch.Restart();
        await fresh();
        streamed.Add(watch.Elapsed.TotalMilliseconds);
    }
    return (first, Median(warm), Median(streamed), top);
}

Func<bool, Func<Task<string>>> Similar(Func<QueryProfile> query) => useCache =>
{
    var search = new SimilaritySearch(dataSource, useCache);
    var profile = query();
    return async () => { var r = await search.SearchAsync(profile); return $"{r.Ranked[0].SequenceId} ({r.Ranked[0].Score:F4})"; };
};
Func<bool, Func<Task<string>>> ByRecord(long id) => useCache =>
{
    var search = new SimilaritySearch(dataSource, useCache);
    return async () => { var r = await search.SearchRecordAsync(id); return $"{r.Ranked[0].SequenceId} ({r.Ranked[0].Score:F4})"; };
};
Func<bool, Func<Task<string>>> Localize(QueryProfile fragment, long? inRecord = null) => useCache =>
{
    var localization = new Localization(dataSource, useCache);
    return async () =>
    {
        var r = inRecord is long id ? await localization.LocalizeInRecordAsync(fragment, id) : await localization.LocalizeAsync(fragment);
        return $"{r.Hits[0].SequenceId} at {r.Hits[0].OffsetBp:N0} {r.Hits[0].Strand} ({r.Hits[0].Score:F4})";
    };
};

QueryProfile Profiled(string dna, out double ms)
{
    var watch = Stopwatch.StartNew();
    var profile = QueryProfile.FromSequence(dna);
    ms = watch.Elapsed.TotalMilliseconds;
    return profile;
}

var k12Query = Profiled(k12, out var k12Ms);
var sarsQuery = Profiled(sars, out var sarsMs);
var fragment10 = Profiled(o157.Substring(2_000_123, 10_000), out var fragment10Ms);
var fragment100 = Profiled(o157.Substring(3_000_321, 100_000), out var fragment100Ms);
var fragment2 = Profiled(o157.Substring(1_234_567, 2_000), out var fragment2Ms);

var queries = new (string Kind, string Query, long LengthBp, double ProfileMs, Func<bool, Func<Task<string>>> Engine)[]
{
    ("find-similar (FASTA)", "E. coli K-12 MG1655", k12.Length, k12Ms, Similar(() => k12Query)),
    ("find-similar (FASTA)", "SARS-CoV-2 Wuhan-Hu-1", sars.Length, sarsMs, Similar(() => sarsQuery)),
    ("find-similar (stored record)", $"E. coli O157:H7 Sakai, record {o157Id}", o157.Length, 0, ByRecord(o157Id)),
    ("localize (catalogue)", "O157:H7 fragment at 2,000,123", 10_000, fragment10Ms, Localize(fragment10)),
    ("localize (catalogue)", "O157:H7 fragment at 3,000,321", 100_000, fragment100Ms, Localize(fragment100)),
    ("localize (named record)", $"O157:H7 fragment at 1,234,567 in record {o157Id}", 2_000, fragment2Ms, Localize(fragment2, o157Id)),
};

using var report = File.CreateText(output);
report.WriteLine("Kind, Query, LengthBp, ProfileMs, CachedFirstMs, CachedWarmMs, StreamedMs, TopHit, Threads");
foreach (var (kind, name, lengthBp, profileMs, engine) in queries)
{
    var (first, warm, streamed, top) = await Time(engine);
    report.WriteLine($"{kind}, \"{name}\", {lengthBp}, {profileMs:F0}, {first:F0}, {warm:F0}, {streamed:F0}, \"{top}\", {Environment.ProcessorCount}");
    report.Flush();
    Console.WriteLine($"{kind,-30} {name,-50} profile {profileMs,6:F0} ms  first {first,6:F0}  warm {warm,6:F0}  streamed {streamed,6:F0}  top {top}");
}
