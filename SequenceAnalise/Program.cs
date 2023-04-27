// See https://aka.ms/new-console-template for more information
using SequenceAnalyses;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Threading.Tasks;

Func<string, DNA, DNA, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms)> PrintComparision = (string prefix, DNA dna, DNA dna2) =>
{
    Stopwatch sw = new Stopwatch();
    sw.Start();

    //Console.WriteLine($"{prefix} CAT       {CATProfil.Compare(dna.CATProfil, dna2.CATProfil)} {sw.Elapsed.TotalMilliseconds}ms"); sw.Restart();
    //Console.WriteLine($"{prefix} CAT fixed {CATProfil.CompareFixed(dna.CATProfil, dna2.CATProfil)} {sw.Elapsed.TotalMilliseconds}ms"); sw.Restart();
    //Console.WriteLine($"{prefix} NEE       {dna.NeedlemanWunschResult(dna2.DnaString)} {sw.Elapsed.TotalMilliseconds}ms"); sw.Restart();

    double CATComparison = CATProfile.Compare(dna.CATProfile, dna2.CATProfile);
    double CATms = sw.Elapsed.TotalMilliseconds; sw.Restart();
    double CATComparisonFixed = CATProfile.CompareFixed(dna.CATProfile, dna2.CATProfile);
    double CATFixedms = sw.Elapsed.TotalMilliseconds; sw.Restart();
    double NeedelmanComparison = NeedlemanWunsch.Calculate(dna.DnaString, dna2.DnaString);
    double Needelmanms = sw.Elapsed.TotalMilliseconds; sw.Restart();

    return (CATComparison, CATms, CATComparisonFixed, CATFixedms, NeedelmanComparison, Needelmanms);
};

Func<(double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms), string> Format = ((double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms) res) =>
{
    return $"{res.CATComparison}, {res.CATms}, {res.CATComparisonFixed},{res.CATFixedms}, {res.NeedelmanComparison}, {res.NeedelmanComparison},{res.Needelmanms}";
};

Func<DNA, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms)> ExactSame = (DNA dna) => PrintComparision(nameof(ExactSame), dna, dna);

Func<DNA, double, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms)> SubsequenceFirstHalf = (DNA dna, double lenght) =>
{
    var dna2 = new DNA(dna.DnaString.Substring(0, (int)lenght));

    return PrintComparision(nameof(SubsequenceFirstHalf), dna, dna2);
};

Func<DNA, double, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms)> SubsequenceSecondHalf = (DNA dna, double lenght) =>
{
    var dna2 = new DNA(dna.DnaString.Substring((int)(dna.DnaString.Length - lenght)));
    return PrintComparision(nameof(SubsequenceSecondHalf), dna, dna2);
};

Func<DNA, double, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, double NeedelmanComparison, double Needelmanms)> SubsequenceInTheMiddle = (DNA dna, double lenght) =>
{
    var dna2 = new DNA(dna.DnaString.Substring((int)((dna.DnaString.Length - lenght) / 2), (int)lenght));
    return PrintComparision(nameof(SubsequenceInTheMiddle), dna, dna2);
};

Func<int, List<int>, List<int>, List<string>> generateResultsForHalfs = (int seq, List<int> halfs, List<int> random) =>
{
    List<string> result = new List<string>();
    for (int i = 0; i < 10; i++)
    {
        var dna = new DNA(Generator.GetRandomDNA(seq));
        var exactSame = ExactSame(dna);
        result.Add($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Format(exactSame)}");

        foreach (var half in halfs)
        {
            var firstHalf = SubsequenceFirstHalf(dna, half);
            var secondHalf = SubsequenceSecondHalf(dna, half);
            var middle = SubsequenceInTheMiddle(dna, half);

            result.Add($"{dna.DnaString.Length}, {half}, FirstHalf, {Format(firstHalf)}");
            result.Add($"{dna.DnaString.Length}, {half}, SecondHalf, {Format(secondHalf)}");
            result.Add($"{dna.DnaString.Length}, {half}, Middle, {Format(middle)}");
        }

        foreach (var rnd in random)
        {
            var dna2 = new DNA(Generator.GetRandomDNA(rnd));
            var aproximate = PrintComparision("Aproximate", dna, dna2);
            result.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Format(aproximate)}");
        }
    }
    return result;
};

Func<int, int, Task> examinCollisions = async (int sequensCount, int sequensLenght) =>
{
    ConcurrentDictionary<string, (int actualMatches, int catMatches)> results = new ConcurrentDictionary<string, (int actualMatches, int catMatches)>();
    List<CATProfile> profiles = new List<CATProfile>();
    for (int i = 0; i < sequensCount; i++)
    {
        var str = Generator.GetRandomDNA(sequensLenght);
        results[str] = (0, 0);

        profiles.Add(new CATProfile(str));
    }

    var totalRecords = 0l;
    foreach (var permutation in Search.GetPermutations(sequensLenght))
    {
        totalRecords++;
        var profile = new CATProfile(permutation);
        foreach (var random in profiles)
        {
            var catComparison = CATProfile.Compare(profile, random);
            if (catComparison == 1d)
            {
                var actualMatches = profile.AreEqual(random) ? 1 : 0;
                results.AddOrUpdate(random.DnaString, (key) => (actualMatches, 1), (key, val) => (val.actualMatches + actualMatches, val.catMatches + 1));
            }
        }
    }

    foreach (var res in results)
    {
        Console.WriteLine($"Total {res.Key} {sequensLenght} cat matches {res.Value.catMatches} of {totalRecords} {res.Value.catMatches / (double)totalRecords} Actual matches: {res.Value.actualMatches}");
    }

    var catMatchesAverage = results.Values.Average(x => x.catMatches);
    var actualMatchesAverage = results.Values.Average(x => x.actualMatches);

    Console.WriteLine($"Total {sequensLenght} cat matches {catMatchesAverage} of {totalRecords} {catMatchesAverage / (double)totalRecords} Actual matches: {actualMatchesAverage} of {actualMatchesAverage / (double)totalRecords}");
};

Func<int, int, Task> generateColisionReport = async (int sequensCount, int sequensLenght) =>
{
    ConcurrentDictionary<string, double> results = new ConcurrentDictionary<string, double>();
    ConcurrentDictionary<string, ConcurrentDictionary<string, int>> labaledResults = new ConcurrentDictionary<string, ConcurrentDictionary<string, int>>();
    var labels = new Dictionary<string, List<double>>()
                        {
                            { "0.0", new List<double>()},
                            { "0.1", new List<double>()},
                            { "0.2", new List<double>()},
                            { "0.3", new List<double>()},
                            { "0.4", new List<double>()},
                            { "0.5", new List<double>()},
                            { "0.6", new List<double>()},
                            { "0.7", new List<double>()},
                            { "0.8", new List<double>()},
                            { "0.9", new List<double>()},
                            { "1.0", new List<double>()},
                        };
    List<CATProfile> randomDNAProfiles = new List<CATProfile>();
    for (int i = 0; i < sequensCount; i++)
    {
        var str = Generator.GetRandomDNA(sequensLenght);
        results[str] = 0;
        labaledResults[str] = new ConcurrentDictionary<string, int>(new Dictionary<string, int>()
                        {
                            { "0.0", 0},
                            { "0.1", 0},
                            { "0.2", 0},
                            { "0.3", 0},
                            { "0.4", 0},
                            { "0.5", 0},
                            { "0.6", 0},
                            { "0.7", 0},
                            { "0.8", 0},
                            { "0.9", 0},
                            { "1.0", 0},
                        });
        randomDNAProfiles.Add(new CATProfile(str));
    }

    var totalRecords = 0l;
    foreach (var permutation in Search.GetPermutations(sequensLenght))
    {
        var profile = new CATProfile(permutation);
        foreach (var random in randomDNAProfiles)
        {
            var catComparison = CATProfile.Compare(profile, random);
            if (catComparison == 1d)
            {
                var res = NeedlemanWunsch.Calculate(random.DnaString, permutation);
                var label = res.ToString("0.000").Substring(0, 3);
                var labelResult = labaledResults.GetOrAdd(random.DnaString, (key) => new ConcurrentDictionary<string, int>());
                labelResult.AddOrUpdate(label, (key) => 1, (key, value) => value + 1);
            }
        }
    }

    using var result = File.CreateText($"SearchResults{sequensLenght}.csv");
    result.WriteLine($"dna string ,totalRecords, cat matches, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0");

    var totals = new List<double>();
    foreach (var res in labaledResults)
    {
        double total = res.Value.Values.Sum();
        totals.Add(total);
        result.WriteLine($"{res.Key}, {totalRecords}, {total}, {string.Join(",", res.Value.OrderBy(x => x.Key).Select(x => x.Value * 100 / total))}");
        foreach (var label in labels)
        {
            label.Value.Add(res.Value[label.Key] * 100 / total);
        }
        result.Flush();
        //Console.WriteLine($"{res.Key}, {totalRecords}, {res.Value.Values.Sum()}, {string.Join(",", res.Value.Values.Select(x => x * 100 / total))}");
        //Console.WriteLine($"Total {seqLength} cat matches {res.Value.catMatches} of {totalRecords} {res.Value.catMatches / (double)totalRecords} Actual matches: {res.Value.actualMatches}");
    }

    result.WriteLine($"Average:, {100 * totals.Average() / totalRecords}, {totals.Average()}, {string.Join(",", labels.OrderBy(x => x.Key).Select(x => x.Value.Average()))}, {labels.OrderBy(x => x.Key).Select(x => x.Value.Average()).Sum()}");
};

await generateColisionReport(100, 11);

var seqLength = new List<int>() { 100, 1000 };
var experimentNum = 90;
var tasks = new List<Task<KeyValuePair<string, List<string>>>>();
foreach (var seq in seqLength)
{
    var ssl10 = (int)(seq * 0.1);
    var ssl30 = (int)(seq * 0.3);
    var ssl50 = (int)(seq * 0.5);

    var ssl70 = (int)(seq * 0.7);
    var ssl90 = (int)(seq * 0.9);
    var ssl97 = (int)(seq * 0.97);
    var halfs = new List<int>() { ssl10, ssl30, ssl50, ssl70, ssl90, ssl97, };
    var random = new List<int> { ssl10, ssl30, ssl50, ssl70, ssl90, ssl97, seq };
    Console.WriteLine($"DNA 1 = {seq}");

    //var tasks = new Task<List<string>>[10];
    Console.WriteLine($"Start Tasks");
    for (int i = 0; i < experimentNum; i++)
    {
        tasks.Add(Task.Factory.StartNew(() =>
        {
            var res = generateResultsForHalfs(seq, halfs, random);
            return new KeyValuePair<string, List<string>>($"seq", res);
        }));
    }
}

var results = await Task.WhenAll(tasks.ToArray());
using (StreamWriter resultFile = File.AppendText("Examin.csv"))
{
    foreach (var res in results)
    {
        res.Value.ForEach(x => resultFile.WriteLine($"{x}"));
        resultFile.Flush();
    }
}