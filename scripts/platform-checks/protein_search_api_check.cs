#:property PublishAot=false
using System.Diagnostics;
using System.Globalization;
using System.Net.Http.Headers;
using System.Net.Http.Json;
using System.Text;
using System.Text.Json;

// Protein search through the API (2026-11-05), the card's end-to-end stop criteria: real BV-BRC proteins posted to
// POST /api/searches/protein one at a time, run by the API's worker on the local catalogue, read back from GET /api/searches/{id}.
// Queries: from the same seeded genomes as protein_search_check.cs (the first genome per kind it lists, by BV-BRC id below), 2 proteins of
// at least 500 aa and 1 of 333-499 aa per genome, as BV-BRC serves them; truth is the CDS's own record, start and PGfam.
// Rule fixed before the run: for proteins of at least 500 aa, the right PGfam at rank 1 for >= 95% of those with one, and the median
// query <= 60 s from POST to done.
// Args: API base URL, a file holding the API key (the Ingest's add-user output), output path (default ProteinSearch_Api.csv).
string[] genomes = ["2926711.3", "34094.20", "1628392.3", "1421.62", "3038546.9", "332614.177", "2866381.3", "2884875.3",
    "2903779.3", "3126101.3", "3056380.3", "2743999.6"];
var api = args[0].TrimEnd('/');
var key = File.ReadAllText(args[1]).Split("shown only now: ")[^1].Trim();
var output = args.Length > 2 ? args[2] : "ProteinSearch_Api.csv";

using var bvbrc = new HttpClient { Timeout = TimeSpan.FromMinutes(5) };
bvbrc.DefaultRequestHeaders.UserAgent.ParseAdd("CATSequenceAnalysis-platform-check/1.0");
bvbrc.DefaultRequestHeaders.Accept.ParseAdd("application/json");
using var client = new HttpClient { BaseAddress = new Uri(api), Timeout = TimeSpan.FromMinutes(5) };
client.DefaultRequestHeaders.Authorization = new AuthenticationHeaderValue("Bearer", key);

async Task<JsonElement> GetAsync(string url)
{
    await Task.Delay(1000);
    using var document = JsonDocument.Parse(await bvbrc.GetStringAsync(url));
    return document.RootElement.Clone();
}

var cases = new List<(string Genome, string Band, string Accession, int Start, string? Family, string Protein)>();
foreach (var genome in genomes)
{
    var features = (await GetAsync($"https://www.bv-brc.org/api/genome_feature/?and(eq(genome_id,{genome}),eq(feature_type,CDS),eq(annotation,PATRIC)," +
        "ge(aa_length,333))&select(accession,start,aa_length,pgfam_id,aa_sequence_md5)&sort(+feature_id)&limit(25000)"))
        .EnumerateArray().Where(f => f.TryGetProperty("aa_sequence_md5", out _)).ToList();
    var picker = new Random(20261105);
    var picked = features.Where(f => f.GetProperty("aa_length").GetInt32() >= 500).OrderBy(_ => picker.Next()).Take(2)
        .Concat(features.Where(f => f.GetProperty("aa_length").GetInt32() < 500).OrderBy(_ => picker.Next()).Take(1)).ToList();
    if (picked.Count == 0) continue;
    var md5s = string.Join(",", picked.Select(f => f.GetProperty("aa_sequence_md5").GetString()).Distinct());
    var proteins = (await GetAsync($"https://www.bv-brc.org/api/feature_sequence/?in(md5,({md5s}))&select(md5,sequence)&limit(100)"))
        .EnumerateArray().ToDictionary(p => p.GetProperty("md5").GetString()!, p => p.GetProperty("sequence").GetString()!);
    foreach (var f in picked)
        if (proteins.TryGetValue(f.GetProperty("aa_sequence_md5").GetString()!, out var protein) && !protein.TrimEnd('*').Contains('*'))
            cases.Add((genome, f.GetProperty("aa_length").GetInt32() >= 500 ? ">= 500 aa" : "333-499 aa", f.GetProperty("accession").GetString()!,
                f.GetProperty("start").GetInt32() - 1, f.TryGetProperty("pgfam_id", out var family) ? family.GetString() : null, protein));
}
Console.WriteLine($"{cases.Count} proteins from {genomes.Length} genomes.");

var rows = new List<(string Band, bool Exact, bool? FamilyTop, bool Weak, double Seconds, string Status)>();
foreach (var c in cases)
{
    var watch = Stopwatch.StartNew();
    var posted = await client.PostAsync("/api/searches/protein?top=10", new StringContent($">q\n{c.Protein}\n", Encoding.ASCII, "text/x-fasta"));
    if ((int)posted.StatusCode != 202)
    {
        rows.Add((c.Band, false, null, false, 0, $"HTTP {(int)posted.StatusCode}"));
        continue;
    }
    var href = (await posted.Content.ReadFromJsonAsync<JsonElement>()).GetProperty("href").GetString()!;
    JsonElement run;
    while (true)
    {
        run = await client.GetFromJsonAsync<JsonElement>(href);
        var status = run.GetProperty("status").GetString();
        if (status is "Done" or "Failed") break;
        await Task.Delay(500);
    }
    double seconds = watch.Elapsed.TotalSeconds;
    if (run.GetProperty("status").GetString() != "Done")
    {
        rows.Add((c.Band, false, null, false, seconds, "Failed: " + run.GetProperty("error").GetString()));
        continue;
    }
    var top = run.GetProperty("hits").GetProperty("items")[0];
    bool exact = top.GetProperty("sequence").GetProperty("recordId").GetString() == c.Accession && top.GetProperty("geneStartBp").GetInt32() == c.Start;
    bool? familyTop = c.Family is null ? null : exact || (top.TryGetProperty("familyId", out var f) && f.GetString() == c.Family);
    bool weak = run.GetProperty("summary").TryGetProperty("weakSignal", out var w) && w.GetBoolean();
    rows.Add((c.Band, exact, familyTop, weak, seconds, "Done"));
    Console.WriteLine($"  {c.Genome} {c.Band} {c.Protein.Length} aa: exact {exact}, family {familyTop?.ToString() ?? "-"}, weak {weak}, {seconds:F1} s");
}

string Percent(int part, int whole) => whole == 0 ? "" : (100.0 * part / whole).ToString("F1", CultureInfo.InvariantCulture) + "%";
var lines = new List<string> { "Band,Proteins,Done,Exact top-1,With PGfam,Family top-1,Flagged weak,Median s,Max s" };
foreach (var band in new[] { ">= 500 aa", "333-499 aa" })
{
    var cell = rows.Where(r => r.Band == band).ToList();
    var done = cell.Where(r => r.Status == "Done").ToList();
    var withFamily = done.Where(r => r.FamilyTop is not null).ToList();
    var times = done.Select(r => r.Seconds).Order().ToList();
    lines.Add(string.Join(",", band, cell.Count, done.Count, Percent(done.Count(r => r.Exact), done.Count), withFamily.Count,
        Percent(withFamily.Count(r => r.FamilyTop == true), withFamily.Count), Percent(done.Count(r => r.Weak), done.Count),
        times.Count == 0 ? "" : times[times.Count / 2].ToString("F1", CultureInfo.InvariantCulture),
        times.Count == 0 ? "" : times[^1].ToString("F1", CultureInfo.InvariantCulture)));
}
var long500 = rows.Where(r => r.Band == ">= 500 aa" && r.Status == "Done").ToList();
var familyLong = long500.Where(r => r.FamilyTop is not null).ToList();
double share = 100.0 * familyLong.Count(r => r.FamilyTop == true) / familyLong.Count;
var allTimes = rows.Where(r => r.Status == "Done").Select(r => r.Seconds).Order().ToList();
double median = allTimes[allTimes.Count / 2];
string verdict = share >= 95 && median <= 60 && rows.All(r => r.Status == "Done") ? "pass" : "stop";
lines.Add("");
lines.Add(string.Create(CultureInfo.InvariantCulture, $"Rule (>= 500 aa family top-1 >= 95% and median query <= 60 s; every job done),{share:F1}% and {median:F1} s: {verdict}"));
foreach (var failed in rows.Where(r => r.Status != "Done")) lines.Add($"Not done,{failed.Band},{failed.Status}");
File.WriteAllLines(output, lines);
Console.WriteLine(string.Join(Environment.NewLine, lines));
