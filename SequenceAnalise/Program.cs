// See https://aka.ms/new-console-template for more information
using MathNet.Numerics;
using SequenceAnalyses;
using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection.Emit;
using System.Runtime.Intrinsics.X86;
using System.Runtime.Serialization.Formatters.Binary;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

Console.WriteLine("Hello, World!");

IEnumerable<string> GetPermutations(int sequenceLenght)
{
    var values = new char[] { 'A', 'C', 'G', 'T' };
    // var values = new Dictionary<int, object[]>(parameters.Select((p, index) => new KeyValuePair<int, object[]>(index, p.GetCachedValue(handler).ToArray())).OrderBy(k => k.Key));
    //var lenghts = values.Select(v => v.Value.Length);
    var permutations = Math.Pow(values.Length, sequenceLenght);// values.Any() ? lenghts.Aggregate(1, (acc, v) => acc * v) : 0;
    var counter = Enumerable.Repeat(0, sequenceLenght).ToArray(); //  values.Select(v => 0).ToArray();

    for (int i = 0; i < permutations; i++)
    {
        //var parameterValues = new Dictionary<string, object>();
        var sequence = new char[sequenceLenght];
        for (int j = 0; j < sequenceLenght; j++)
        {
            sequence[j] = values[counter[j]];
        }

        yield return new string(sequence);
        //foreach (var val in values)
        //{
        //    parameterValues[parameters[val.Key].Name] = val.Value[counter[val.Key]];
        //}

        //yield return parameterValues;
        // update counter
        counter[counter.Length - 1]++;
        for (int j = counter.Length - 1; j > 0; j--)
        {
            if (counter[j] >= values.Length)
            {
                counter[j] = 0;
                counter[j - 1]++;
            }
        }
    }
}

Dictionary<string, int> columns = new Dictionary<string, int>
{
    { "DNA 1", 0 },
    { "DNA 2", 1 },
    { "DNA 1 Length", 2},
    { "DNA 2 Length", 3},
    { "Coparison Type", 4},
    { "CAT", 5},
    { "CAT Fixed", 6},
    { "Needelman", 7},
    { "delta CAT", 8},
    { "delta CAT fixed", 9}
};

Func<string, DNA, DNA, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms)> PrintComparision = (string prefix, DNA dna, DNA dna2) =>
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
    (double exactMatcha, double gaps) NeedelmanComparison = dna.NeedlemanWunschResult(dna2.DnaString);
    double Needelmanms = sw.Elapsed.TotalMilliseconds; sw.Restart();

    return (CATComparison, CATms, CATComparisonFixed, CATFixedms, NeedelmanComparison, Needelmanms);
};

Func<(double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms), string> Format = ((double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms) res) =>
{
    return $"{res.CATComparison}, {res.CATms}, {res.CATComparisonFixed},{res.CATFixedms}, {res.NeedelmanComparison.exactMatcha}, {res.NeedelmanComparison.gaps},{res.Needelmanms}";
};

Func<DNA, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms)> ExactSame = (DNA dna) => PrintComparision(nameof(ExactSame), dna, dna);

Func<DNA, double, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms)> SubsequenceFirstHalf = (DNA dna, double lenght) =>
{
    var dna2 = new DNA(dna.DnaString.Substring(0, (int)lenght));

    return PrintComparision(nameof(SubsequenceFirstHalf), dna, dna2);
};

Func<DNA, double, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms)> SubsequenceSecondHalf = (DNA dna, double lenght) =>
{
    var dna2 = new DNA(dna.DnaString.Substring((int)(dna.DnaString.Length - lenght)));
    return PrintComparision(nameof(SubsequenceSecondHalf), dna, dna2);
};

Func<DNA, double, (double CATComparison, double CATms, double CATComparisonFixed, double CATFixedms, (double exactMatcha, double gaps) NeedelmanComparison, double Needelmanms)> SubsequenceInTheMiddle = (DNA dna, double lenght) =>
{
    var dna2 = new DNA(dna.DnaString.Substring((int)((dna.DnaString.Length - lenght) / 2), (int)lenght));
    return PrintComparision(nameof(SubsequenceInTheMiddle), dna, dna2);
};

//using (StreamWriter resultFile = File.AppendText($"Test.csv"))
//    foreach (var sequensLenght in Search.GetPermutations(15))
//    {
//        try
//        {
//            var catProfile = new CATProfile(sequensLenght);
//            if (catProfile.a.Cos > 1 || catProfile.a.H == 0 || catProfile.t.Cos > 1 || catProfile.t.H == 0)
//                resultFile.WriteLine(sequensLenght);
//        }
//        catch (Exception ex)
//        {
//            resultFile.WriteLine(sequensLenght);
//        }
//    }

//var refSeq = new DNA("GAATTCAGTTA");
//var alineSeq = new DNA("GGATCGA");
//refSeq.NeedlemanWunschResult(alineSeq.DnaString);

//Dictionary<string, DNA> profiles = new Dictionary<string, DNA>();
//foreach (string sequence in GetPermutations(6))
//{
//    profiles[sequence] = new DNA(sequence);
//}

//using (StreamWriter resultFile = File.CreateText("Test.csv"))
//{
//    foreach (var dna in profiles)
//    {
//        resultFile.WriteLine($"{dna.Key}, {dna.Value.GpsProfil}, {profiles.Values.Any(x => x != dna.Value && dna.Value.GpsProfil == x.GpsProfil)}");
//    }
//}

//var aaga = new DNA("AAGA");
//var tatc = new DNA("TATC");

// generate results in some table

//Action<int, int> ex = (int count, int sequenceLenght) =>
//{
//    using (StreamWriter resultFile = File.AppendText($"Test_{count}_{sequenceLenght}.csv"))
//    {
//        var ssl10 = (int)(sequenceLenght * 0.1);
//        var ssl30 = (int)(sequenceLenght * 0.3);
//        var ssl50 = (int)(sequenceLenght * 0.5);
//        var ssl70 = (int)(sequenceLenght * 0.7);
//        var ssl90 = (int)(sequenceLenght * 0.9);

//        var halfs = new List<int>() { ssl50, ssl70, ssl90 };//{ ssl10, ssl30, ssl50 };

//        var ssl97 = (int)(sequenceLenght * 0.97);
//        var random = new List<int> { ssl70, ssl90, ssl97, sequenceLenght };//{ ssl10, ssl30, ssl50, ssl70, ssl90, ssl97, sequensLenght };

//        var i = 0;
//        foreach (var profile in Search.GetProfile(count, sequenceLenght))
//        {
//            Console.WriteLine($"Experiment = {i} {count}_{sequenceLenght}");
//            var dna = new DNA(profile.DnaString);
//            var exactSame = ExactSame(dna);
//            resultFile.WriteLine($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Format(exactSame)}");

//            foreach (var half in halfs)
//            {
//                var firstHalf = SubsequenceFirstHalf(dna, half);
//                var secondHalf = SubsequenceSecondHalf(dna, half);
//                var middle = SubsequenceInTheMiddle(dna, half);

//                resultFile.WriteLine($"{dna.DnaString.Length}, {half}, FirstHalf, {Format(firstHalf)}");
//                resultFile.WriteLine($"{dna.DnaString.Length}, {half}, SecondHalf, {Format(secondHalf)}");
//                resultFile.WriteLine($"{dna.DnaString.Length}, {half}, Middle, {Format(middle)}");
//            }

//            foreach (var rnd in random)
//            {
//                var dna2 = new DNA(Generator.GetRandomDNA(rnd));
//                var aproximate = PrintComparision("Aproximate", dna, dna2);
//                resultFile.WriteLine($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Format(aproximate)}");
//            }
//            i++;
//        }
//    }
//};

//var tsks = new List<Task>()
//{
//    Task.Run(()=> Search.GeneratePermutationData(10)),
//    Task.Run(()=> Search.GeneratePermutationData(11)),
//    Task.Run(()=> Search.GeneratePermutationData(12)),
//    Task.Run(()=> Search.GeneratePermutationData(13)),
//};
////{
////    Task.Run(() => Search.GenerateRandomSeq(200, 1 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 2 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 3 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 4 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 5 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 6 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 7 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 8 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 9 * 10000)),
////    Task.Run(() => Search.GenerateRandomSeq(200, 10 * 10000)),
////};
//{
//    Task.Run(() => ex(200, 1 * 10000)),
//    Task.Run(() => ex(200, 2 * 10000)),
//    Task.Run(() => ex(200, 3 * 10000)),
//    Task.Run(() => ex(200, 4 * 10000)),
//    Task.Run(() => ex(200, 5 * 10000)),
//    Task.Run(() => ex(200, 6 * 10000)),
//    Task.Run(() => ex(200, 7 * 10000)),
//    Task.Run(() => ex(200, 8 * 10000)),
//    Task.Run(() => ex(200, 9 * 10000)),
//    Task.Run(() => ex(200, 10 * 10000)),
//};

//Task.WaitAll(tsks.ToArray());
//return;
//Search.GeneratePermutationData(12);

//return;
//var seqLength = 10;
//using (StreamWriter resultFile = File.CreateText("TestAT.csv"))
//{
//    foreach (var profile1 in Search.GetProfile(12, "file1.dat"))
//    {
//        foreach (var profile2 in Search.GetProfile(12, "file2.dat"))
//        {
//            var compareTest = CATProfile.CompareTest(profile1, profile2);
//            if (compareTest > 1)
//            {
//                resultFile.WriteLine($"{profile1.DnaString},{profile2.DnaString},{compareTest}");
//            }
//        }
//    }
//    resultFile.Flush();
//}

//return;
//NeedlemanWunsch.Calculate("AAACCGACCA", "CTCATTTTGT");
//return;

Func<int, Action<string>, Task<double>> iterateInParalel = async (int sequensLenght, Action<string> iterator) =>
{
    var taskAtoC = Task.Run<double>(() =>
    {
        var totalRecords = 0d;

        foreach (var permutation in Search.GetPermutationsAtoC(sequensLenght, false))
        {
            totalRecords++;
            iterator(permutation);
        }

        return totalRecords;
    });
    var taskCtoG = Task.Run<double>(() =>
    {
        var totalRecords = 0d;

        foreach (var permutation in Search.GetPermutationsCtoG(sequensLenght, false))
        {
            totalRecords++;
            iterator(permutation);
        }

        return totalRecords;
    });
    var taskGtoT = Task.Run<double>(() =>
    {
        var totalRecords = 0d;

        foreach (var permutation in Search.GetPermutationsGtoT(sequensLenght, true))
        {
            totalRecords++;
            iterator(permutation);
        }

        return totalRecords;
    });

    return (await Task.WhenAll(taskAtoC, taskCtoG, taskGtoT)).Sum();
};

Func<int, Task> examinContribution = async (int sequensLenght) =>
{
    ConcurrentDictionary<string, bool> results = new ConcurrentDictionary<string, bool>();

    var totalRecords = await iterateInParalel(sequensLenght, (permutation) =>
    {
        var profile = new CATProfile(permutation);
        if (profile.CanContribute() == false)
        {
            results.AddOrUpdate(permutation, (key) => false, (key, val) => false);
        }
    });

    foreach (var res in results)
    {
        Console.WriteLine($"{res.Key} cannot contribute");
    }
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

    var totalRecords = await iterateInParalel(sequensLenght, (permutation) =>
    {
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
    });

    foreach (var res in results)
    {
        Console.WriteLine($"Total {res.Key} {sequensLenght} cat matches {res.Value.catMatches} of {totalRecords} {res.Value.catMatches / (double)totalRecords} Actual matches: {res.Value.actualMatches}");
    }

    var catMatchesAverage = results.Values.Average(x => x.catMatches);
    var actualMatchesAverage = results.Values.Average(x => x.actualMatches);

    Console.WriteLine($"Total {sequensLenght} cat matches {catMatchesAverage} of {totalRecords} {catMatchesAverage / (double)totalRecords} Actual matches: {actualMatchesAverage} of {actualMatchesAverage / (double)totalRecords}");
};

Func<string, Task<ConcurrentDictionary<string, int>>> examinSequenseCollision = async (string dnaString) =>
{
    ConcurrentDictionary<string, int> collisions = new ConcurrentDictionary<string, int>(new Dictionary<string, int>()
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
    var examinDNA = new CATProfile(dnaString);
    Console.WriteLine($"{dnaString}");
    var totalRecords = await iterateInParalel(dnaString.Length, (permutation) =>
    {
        var profile = new CATProfile(permutation);

        var catComparison = CATProfile.Compare(profile, examinDNA);
        if (catComparison == 1d)
        {
            var res = NeedlemanWunsch.Calculate(dnaString, permutation);
            var label = res.ToString("0.000").Substring(0, 3);
            collisions.AddOrUpdate(label, (key) => 1, (key, value) => value + 1);
        }
    });
    return collisions;
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

    double totalRecords = await iterateInParalel(sequensLenght, (permutation) =>
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
    });

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

Action<int, int> comparison = (int sequensCount, int sequensLenght) =>
{
    for (int i = 0; i < sequensCount; i++)
    {
        //var baseDna = Generator.GetRandomDNA(sequensLenght);
        //var dna1 = "AGT" + baseDna;
        //var dna2 = "GTC" + baseDna;

        var baseDna = Generator.GetRandomDNA(sequensLenght);
        var dna1 = Generator.GetRandomDNA(sequensLenght);
        var dna2 = Generator.GetRandomDNA(sequensLenght);
        var p1 = new CATProfile(dna1);
        var p2 = new CATProfile(dna2);
        var r = CATProfile.Compare(p1, p2);
        var res = NeedlemanWunsch.Calculate(dna1, dna2);
        Console.WriteLine($"{r} expected {res}");
    }
};

Func<int, int, List<int>, List<int>, List<string>> doExperiment = (int sequensLenght, int count, List<int> halfs, List<int> random) =>
{
    Task<List<string>>[] experiments = new Task<List<string>>[count];
    for (int i = 0; i < count; i++)
    {
        experiments[i] = Task.Factory.StartNew(() =>
        {
            List<string> result = new List<string>();
            Console.WriteLine($"sequensLenght {sequensLenght}, experiment {Task.CurrentId} {DateTime.Now.ToString("hh:mm:ss:ffff")}");
            var dna = new DNA(Generator.GetRandomDNA(sequensLenght));
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
            Console.WriteLine($"sequensLenght {sequensLenght}, experiment {Task.CurrentId} {DateTime.Now.ToString("hh:mm:ss:ffff")}");
            return result;
        });
    }
    Task.WhenAll(experiments).Wait();
    return experiments.SelectMany(x => x.Result).ToList();
};

Action<int, int> evaluatePrecision = (int sequensLenght, int permutations) =>
{
    var mainPart = Generator.GetRandomDNA(sequensLenght);
    foreach (var diffPart1 in Search.GetPermutations(permutations))
    {
        foreach (var diffPart2 in Search.GetPermutations(permutations))
        {
            var equalInDiffPart = 0;
            for (int i = 0; i < permutations; i++)
            {
                equalInDiffPart += diffPart1[i] == diffPart2[i] ? 1 : 0;
            }
            var profile1 = new CATProfile(diffPart1 + mainPart);
            var profile2 = new CATProfile(diffPart2 + mainPart);
            var catcopmarison = CATProfile.Compare(profile1, profile2);
            var actual = (equalInDiffPart + mainPart.Length) / (double)(mainPart.Length + permutations);
            Console.WriteLine($"{profile1.DnaString} {profile2.DnaString} expected % {actual} CAT {catcopmarison} delta {actual - catcopmarison}");
        }
    }
};
//var profile1 = new CATProfile("CGCAGGAGCA");
//CATProfile.Compare(profile1, profile1);
//foreach (var r in await examinSequenseCollision("GCCGGTATTG"))
//{
//    Console.WriteLine($"{r.Key} {r.Value}");
//}
//await examinContribution(10);
//await examinCollisions(100, 11);

var tasksgenerateColisionReport = new List<Task>()
{
//Task.Run(() =>{ Console.WriteLine($"Started {10}");  return generateColisionReport(1000, 10);}),
Task.Run(() =>{ Console.WriteLine($"Started {11}");  return generateColisionReport(1000, 11);}),
Task.Run(() =>{ Console.WriteLine($"Started {12}");  return generateColisionReport(1000, 12);}),
Task.Run(() =>{ Console.WriteLine($"Started {13}");  return generateColisionReport(1000, 13);}),
Task.Run(() =>{ Console.WriteLine($"Started {15}");  return generateColisionReport(1000, 15);}),
};
Task.WaitAll(tasksgenerateColisionReport.ToArray());

//var profile2 = new CATProfile("CCAG");
//var profile3 = new CATProfile("CCTG");
//var profile4 = new CATProfile("ACGT");
//CATProfile.Compare(profile1, profile2);
//CATProfile.Compare(profile1, profile3);
//CATProfile.Compare(profile2, profile3);
//for (int i = 0; i < 1; i++)
//{
//    evaluatePrecision(0, 4);
//}

//return;

//var seqLength = new List<int>() { 1_000 }; //,, 5_000  10_000  25_000 , 50_000 };
//var experimentNum = 100;
//var tasks = new List<Task<KeyValuePair<string, List<string>>>>();

//foreach (var seq in seqLength)
//{
//    var ssl10 = (int)(seq * 0.1);
//    var ssl30 = (int)(seq * 0.3);
//    var ssl50 = (int)(seq * 0.5);

//    var ssl70 = (int)(seq * 0.7);
//    var ssl90 = (int)(seq * 0.9);
//    var ssl97 = (int)(seq * 0.97);
//    var halfs = new List<int>() { ssl10, ssl30, ssl50, ssl70, ssl90, ssl97, };
//    var random = new List<int> { ssl90, ssl97, seq };
//    var sequensLenght = seq * 1;
//    Console.WriteLine($"DNA 1 = {seq}");

//    tasks.Add(Task.Factory.StartNew(() =>
//    {
//        var res = doExperiment(sequensLenght, experimentNum, halfs, random);
//        return new KeyValuePair<string, List<string>>($"sequensLenght", res);
//    }));
//}
//var results = await Task.WhenAll(tasks.ToArray());
//using (StreamWriter resultFile = File.AppendText("Examin.csv"))
//{
//    foreach (var res in results)
//    {
//        res.Value.ForEach(x => resultFile.WriteLine($"{x}"));
//        resultFile.Flush();
//    }
//}

//ConcurrentDictionary<string, int> sequensLenght = new ConcurrentDictionary<string, int>();
//var total = await iterateInParalel(5, (permutation) => { sequensLenght.AddOrUpdate(permutation, (key) => 1, (key, value) => value++); });
//Console.WriteLine($"{total} {sequensLenght.Count(x => x.Value > 1)}");
//foreach (var res in results)
//{
//    result.WriteLine($"{res.Key}, {totalRecords}, {res.Value.catMatches}, {res.Value.actualMatches}");
//    result.Flush();
//    Console.WriteLine($"Total {seqLength} cat matches {res.Value.catMatches} of {totalRecords} {res.Value.catMatches / (double)totalRecords} Actual matches: {res.Value.actualMatches}");
//}
//await examinContribution(10);
//await examinCollisions(100, 10);

//var result = await examinSequenseCollision("GTACGCCCTT");
//foreach (var res in result)
//    Console.WriteLine($"{res.Key} {res.Value}");
//return;
//ConcurrentDictionary<string, int> sequensLenght = new ConcurrentDictionary<string, int>();
//var tasks = new List<Task>()
//{
//Task.Run(() =>{ Console.WriteLine($"Started {10}");  return generateColisionReport(1000, 10);}),
//Task.Run(() =>{ Console.WriteLine($"Started {11}");  return generateColisionReport(1000, 11);}),
//Task.Run(() =>{ Console.WriteLine($"Started {12}");  return generateColisionReport(1000, 12);}),
//Task.Run(() =>{ Console.WriteLine($"Started {13}");  return generateColisionReport(1000, 13);}),
//Task.Run(() =>{ Console.WriteLine($"Started {15}");  return generateColisionReport(1000, 15);}),
//};
//Task.WaitAll(tasks.ToArray());

//return;

//var totalRecords = Search.GetTotalRecords(seqLength);
//var dnaAt18205 = Search.GetDNAAt(seqLength, 18205);
//var ctaAt18205 = new CATProfil(dnaAt18205);
//var ctaAt18205Slim = ctaAt18205.GetCATProfilSlim();
//Stopwatch strStopwatch = new Stopwatch();
//Stopwatch catStopwatch = new Stopwatch();
//var catMatches = 0;
//var catActual = 0;
//using var result = File.CreateText($"SearchResults{seqLength}.txt");
//Console.WriteLine(dnaAt18205);
//result.WriteLine($"Search: {dnaAt18205}");
//foreach (var profile in Search.GetProfile(10))
//{
//    var slim = profile.GetCATProfilSlim();

//    //var strComparison = dnaAt18205 == profile.DnaString;

//    //if (strComparison)
//    //{
//    //    profile.ToString();
//    //    ctaAt18205.ToString();
//    //    var res = CATProfil.Compare(ctaAt18205, profile);
//    //}
//    //catStopwatch.Start();
//    //var catComparison = CATProfil.Compare(ctaAt18205, profile);
//    //catStopwatch.Stop();
//    //catMatches += catComparison == 1d ? 1 : 0;
//    //if (catComparison == 1d)
//    //{
//    //    var slim = profile.GetCATProfilSlim();
//    //    catActual += slim.AreEqual(ctaAt18205Slim) ? 1 : 0;
//    //    result.WriteLine($"{profile.DnaString} {ctaAt18205Slim.ToString()} {slim.ToString()} {slim.AreEqual(ctaAt18205Slim)}");
//    //Console.WriteLine($" {ctaAt18205.DnaString} {profile.DnaString} {ctaAt18205.DnaString == profile.DnaString}");
//    //}
//}
////result.WriteLine($"Total {seqLength} cat matches {catMatches} of {totalRecords} {catMatches / (double)totalRecords} Actual matches: {catActual}");
////result.Flush();
////Console.WriteLine($"Total string comparisons {strStopwatch.Elapsed.TotalMilliseconds} ms");
////Console.WriteLine($"Total cat comparisons {catStopwatch.Elapsed.TotalMilliseconds} ms");
////Console.WriteLine($"Total {seqLength} cat matches {catMatches} of {totalRecords} {catMatches / (double)totalRecords} Actual matches: {catActual}");
//return;

//Func<int, List<int>, List<int>, List<string>> doExperiment = (int sequensLenght, List<int> halfs, List<int> random) =>
//{
//    List<string> result = new List<string>();
//    for (int i = 0; i < 10; i++)
//    {
//        var dna = new DNA(Generator.GetRandomDNA(sequensLenght));
//        var exactSame = ExactSame(dna);
//        result.Add($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Format(exactSame)}");

//        foreach (var half in halfs)
//        {
//            var firstHalf = SubsequenceFirstHalf(dna, half);
//            var secondHalf = SubsequenceSecondHalf(dna, half);
//            var middle = SubsequenceInTheMiddle(dna, half);

//            result.Add($"{dna.DnaString.Length}, {half}, FirstHalf, {Format(firstHalf)}");
//            result.Add($"{dna.DnaString.Length}, {half}, SecondHalf, {Format(secondHalf)}");
//            result.Add($"{dna.DnaString.Length}, {half}, Middle, {Format(middle)}");
//        }

//        foreach (var rnd in random)
//        {
//            var dna2 = new DNA(Generator.GetRandomDNA(rnd));
//            var aproximate = PrintComparision("Aproximate", dna, dna2);
//            result.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Format(aproximate)}");
//        }
//    }
//    return result;
//};

//using (StreamWriter resultFile = File.AppendText("Test1.csv"))
//{
//    var seqLength = new List<int>() { 50000 };
//    var experimentNum = 90;

//    foreach (var sequensLenght in seqLength)
//    {
//        var ssl10 = (int)(sequensLenght * 0.1);
//        var ssl30 = (int)(sequensLenght * 0.3);
//        var ssl50 = (int)(sequensLenght * 0.5);

//        var ssl70 = (int)(sequensLenght * 0.7);
//        var ssl90 = (int)(sequensLenght * 0.9);
//        var ssl97 = (int)(sequensLenght * 0.97);
//        var halfs = new List<int>() { ssl10, ssl30, ssl50, ssl70, ssl90, ssl97, };
//        var random = new List<int> { ssl10, ssl30, ssl50, ssl70, ssl90, ssl97, sequensLenght };
//        Console.WriteLine($"DNA 1 = {sequensLenght}");

//        //var tasks = new Task<List<string>>[10];
//        Console.WriteLine($"Start Tasks");
//        for (int i = 0; i < experimentNum; i++)
//        {
//            //tasks[i] = Task.Factory.StartNew(() => doExperiment(sequensLenght, halfs, random));
//            Console.WriteLine($"Experiment = {i}");
//            var dna = new DNA(Generator.GetRandomDNA(sequensLenght));
//            var exactSame = ExactSame(dna);
//            resultFile.WriteLine($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Format(exactSame)}");

//            foreach (var half in halfs)
//            {
//                var firstHalf = SubsequenceFirstHalf(dna, half);
//                var secondHalf = SubsequenceSecondHalf(dna, half);
//                var middle = SubsequenceInTheMiddle(dna, half);

//                resultFile.WriteLine($"{dna.DnaString.Length}, {half}, FirstHalf, {Format(firstHalf)}");
//                resultFile.WriteLine($"{dna.DnaString.Length}, {half}, SecondHalf, {Format(secondHalf)}");
//                resultFile.WriteLine($"{dna.DnaString.Length}, {half}, Middle, {Format(middle)}");
//            }

//            foreach (var rnd in random)
//            {
//                var dna2 = new DNA(Generator.GetRandomDNA(rnd));
//                var aproximate = PrintComparision("Aproximate", dna, dna2);
//                resultFile.WriteLine($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Format(aproximate)}");
//            }
//        }
//        //Console.WriteLine($"Wait for results");
//        //var results = await Task.WhenAll(tasks);
//        //Console.WriteLine($"write results");
//        //foreach (var res in results)
//        //{
//        //    for (int i = 0; i < res.Count; i++)
//        //    {
//        //        resultFile.WriteLine(res[i]);
//        //    }
//        //}
//    }
//}

//DNA dna1 = new DNA("GCAGCTGGATTATGGCCATGCGACTGAGGGTGGACTTGAGATACACAGGTTGGTGTGTACCCAATGTTGGTGCTTGCGTGGTTAGAGAGATAGACAAGAGTGCTCGTGACCGTGTCCCTTCTTGAATTGGGGCCCTACAAGTGGACATGTTCGCCTAAACCCTGCCAGCAGGTGTATGAATTTGCCTAGTCAGTGGGTACAGCCGGAGCGCCAGAGCTTTAGTGTCCACCGATTAGTGGGGGCGTCATGTGTGCCCGAGGCGAGTGACAAATAGCTTACTAGGAGAAGGATGGAAGGAGGACTCGCGCTAATCCGCGCGCGGTCAACTTCTGGGCCACACTCGCGACATCGAGCGAGGGTACAAGCGTTTCTTACTTTAACAAATGCCTTACGAAGCCTGTTGGGGCCTCCGGCCGAGAATCCCTAGCCATGCGGGCTGAATACGGCTCGGGACCGACAGCGACTACAATAATAGTACCTACAACTACCAGCTTCTCTCAGTATATTTCGAGTGCCAAACGCTGGACAACTTAGATATCATACTGTCGCGCGGGCCTCCAAATTGGGGCGCATTTATCTATATTCTACCAGCCAGATTACTTCAAATCGCTGTACCCTAGTATAAGATGCTGCCGAGATAAGAGCTTGGTGGTAGCATCGCCATCAG");//Generator.GetRandomDNA(25546));
//DNA dna2 = new DNA(dna1.DnaString.Substring(0, (dna1.DnaString.Length / 4)), 4);//Generator.GetRandomDNA(25546));
//var res = DNA.Profil.Compare(dna2.GpsProfil, dna1.GpsProfil);
//var resNidelman = dna2.NeedlemanWunschResult(dna1.DnaString);
//Console.WriteLine(res);//dna1.GpsProfil.ToString())
//Console.WriteLine(string.Format("Needelman: {1} {0}CATAnalis: {2} {0}", Environment.NewLine, resNidelman, res));

//dna1 = new DNA("GCAGCTGGATTATGGCCATGCGACTGAGGGTGGACTTGAGATACACAGGTTGGTGTGTACCCAATGTTGGTGCTTGCGTGGTTAGAGAGATAGACAAGAGTGCTCGTGACCGTGTCCCTTCTTGAATTGGGGCCCTACAAGTGGACATGTTCGCCTAAACCCTGCCAGCAGGTGTATGAATTTGCCTAGTCAGTGGGTACAGCCGGAGCGCCAGAGCTTTAGTGTCCACCGATTAGTGGGGGCGTCATGTGTGCCCGAGGCGAGTGACAAATAGCTTACTAGGAGAAGGATGGAAGGAGGACTCGCGCTAATCCGCGCGCGGTCAACTTCTGGGCCACACTCGCGACATCGAGCGAGGGTACAAGCGTTTCTTACTTTAACAAATGCCTTACGAAGCCTGTTGGGGCCTCCGGCCGAGAATCCCTAGCCATGCGGGCTGAATACGGCTCGGGACCGACAGCGACTACAATAATAGTACCTACAACTACCAGCTTCTCTCAGTATATTTCGAGTGCCAAACGCTGGACAACTTAGATATCATACTGTCGCGCGGGCCTCCAAATTGGGGCGCATTTATCTATATTCTACCAGCCAGATTACTTCAAATCGCTGTACCCTAGTATAAGATGCTGCCGAGATAAGAGCTTGGTGGTAGCATCGCCATCAG");//Generator.GetRandomDNA(25546));
//dna2 = new DNA(dna1.DnaString.Substring((dna1.DnaString.Length / 4)), 4);
//res = DNA.Profil.Compare(dna2.GpsProfil, dna1.GpsProfil);

//Console.WriteLine(res);//dna1.GpsProfil.ToString())

//dna1 = new DNA("GCAGCTGGATTATGGCCATGCGACTGAGGGTGGACTTGAGATACACAGGTTGGTGTGTACCCAATGTTGGTGCTTGCGTGGTTAGAGAGATAGACAAGAGTGCTCGTGACCGTGTCCCTTCTTGAATTGGGGCCCTACAAGTGGACATGTTCGCCTAAACCCTGCCAGCAGGTGTATGAATTTGCCTAGTCAGTGGGTACAGCCGGAGCGCCAGAGCTTTAGTGTCCACCGATTAGTGGGGGCGTCATGTGTGCCCGAGGCGAGTGACAAATAGCTTACTAGGAGAAGGATGGAAGGAGGACTCGCGCTAATCCGCGCGCGGTCAACTTCTGGGCCACACTCGCGACATCGAGCGAGGGTACAAGCGTTTCTTACTTTAACAAATGCCTTACGAAGCCTGTTGGGGCCTCCGGCCGAGAATCCCTAGCCATGCGGGCTGAATACGGCTCGGGACCGACAGCGACTACAATAATAGTACCTACAACTACCAGCTTCTCTCAGTATATTTCGAGTGCCAAACGCTGGACAACTTAGATATCATACTGTCGCGCGGGCCTCCAAATTGGGGCGCATTTATCTATATTCTACCAGCCAGATTACTTCAAATCGCTGTACCCTAGTATAAGATGCTGCCGAGATAAGAGCTTGGTGGTAGCATCGCCATCAG");//Generator.GetRandomDNA(25546));
//dna2 = new DNA(dna1.DnaString.Substring((dna1.DnaString.Length / 4), dna1.DnaString.Length / 4), 4);
//res = DNA.Profil.Compare(dna2.GpsProfil, dna1.GpsProfil);

//Console.WriteLine(res);//dna1.GpsProfil.ToString())

//dna1 = new DNA("GCAGCTGGATTATGGCCATGCGACTGAGGGTGGACTTGAGATACACAGGTTGGTGTGTACCCAATGTTGGTGCTTGCGTGGTTAGAGAGATAGACAAGAGTGCTCGTGACCGTGTCCCTTCTTGAATTGGGGCCCTACAAGTGGACATGTTCGCCTAAACCCTGCCAGCAGGTGTATGAATTTGCCTAGTCAGTGGGTACAGCCGGAGCGCCAGAGCTTTAGTGTCCACCGATTAGTGGGGGCGTCATGTGTGCCCGAGGCGAGTGACAAATAGCTTACTAGGAGAAGGATGGAAGGAGGACTCGCGCTAATCCGCGCGCGGTCAACTTCTGGGCCACACTCGCGACATCGAGCGAGGGTACAAGCGTTTCTTACTTTAACAAATGCCTTACGAAGCCTGTTGGGGCCTCCGGCCGAGAATCCCTAGCCATGCGGGCTGAATACGGCTCGGGACCGACAGCGACTACAATAATAGTACCTACAACTACCAGCTTCTCTCAGTATATTTCGAGTGCCAAACGCTGGACAACTTAGATATCATACTGTCGCGCGGGCCTCCAAATTGGGGCGCATTTATCTATATTCTACCAGCCAGATTACTTCAAATCGCTGTACCCTAGTATAAGATGCTGCCGAGATAAGAGCTTGGTGGTAGCATCGCCATCAG");//Generator.GetRandomDNA(25546));
//dna2 = new DNA("TCGAGAAGAGATAACCTCGAGACTTACGTCTCAGAGGCACCCGGTTACGTCAGATGATTAGCTGTCGGCGTATCTGCCTGATCACCGATCAATAGCTGTAAGGGGACCTAAATCAGGCAATCGGGTGTAAAGCAGGAGTCGGACTTTGATGTCGCGAGTTTACACAACGCGCAGGCTATGGATCCTCTTGAAAGATAGATGGCTTACCAAGGTCCTGCCTACGCTATCTTCTTTGGTGGAAGTTCCATAAGTGAGCTAGCGCGTCGACATATAGTAGTATTCGACGACAGAGTCGATGGACCTGGCCGCATTCGGGCGCTAGTTACTGTATATCGTATTTGGAACACGGAAATTGCGCGGCCTACCAGATTGTCTAAT");//Generator.GetRandomDNA(25546));
//res = DNA.Profil.Compare(dna1.GpsProfil, dna2.GpsProfil);

//Console.WriteLine(res);//dna1.GpsProfil.ToString())

//dna1 = new DNA("GATTCTTTCGGATGGGGCCTTCAATGTGCTCTTGTATGCAGTAGAATTCATGATCAGCCCGACTCAGGCAGCACCCCTTGGCATTTTCCCCTGGCAGAGTGGGCTCAGTCATGCGCGCCCCCGAGCATGCGCGGTCTGCCTCATAGGAATGTAAGGTATAATGCTCCGAATCTCTCATGCCGCGTATCGTTTGCTTATCGACTAAACCTGTAAAGATCGTGGCACTTTGGTTGGTCCTTAAACACAAGTGATATTATGGTAGTTATGCGCGCGTGAGCGCTACCTTGGTAAATCGTGCGCGCAAAGCCTGCTGCTAACGGGATCTGCGCATTGGGAACGACGGGTTACGACACCAGCCATTTTGAGGACCGGTTCTCCAGTGGACCTCCAGCTAGCCAGAGAGCAAATTCCGGGAGCTGTGCGTAGCATCGCTCTGATATGATAGCCGTGAAACATCACACCATGTATATAAGCACCCAAAAAAGCATCGTATGGACCGCACTTCTAATCGTGGACCCAGGTCGGTTGCACCGGCGCCTGATTACCGTTCCATCAGCTATGCACTTTCATATACTGCTAGATATACTGAAAATTGATCTAGATAGAAAATACACATTCGATTTTTAGACATGACAGTAAATCTCAATAGCCGGCATGGATTGTAGGCGGTTGAGACCGCCAACGCGATTTGGTAACATTGATCCTAGCGAGTGCTGCGCCAATGTTATACGATGACCTAAAT");//Generator.GetRandomDNA(25546));
//dna2 = new DNA("AAGCGGAGGGTTTCAAGATTAAACAAGTTCATTTCGCCCACTAAAACTGGTCCTCTAGCCCCATACACGATATCTACTTTGAGTTGTCTGGATCGCCGGAAAACCGGTCGGGCTGCGGTTAACTGGATGAAATATTAGCAGGTCGAGGTCCAAGGACGTTACTTATTAGGGTTGGCTTCCCCCCCCGGAGAGTTTTCCTTTTAAATGCGGCATGTAGGTGCACTCGTCCTGAAGGAGATAGTGAGCTGGGGGTATGGCGCTGCTTTGCCTCAAGTTGCCTGACGTGTGCGAACCCCCGTTATTACCCGCGCGGGCGCTAAAGGAGCGGTATTTATCCTGTGTGGAC");//Generator.GetRandomDNA(25546));
//res = DNA.Profil.Compare(dna1.GpsProfil, dna2.GpsProfil);

//Console.WriteLine(res);

//dna1 = new DNA("GAATAGCGTCGACGCTATTCACAGCATCAATTGCATGGCCGTTTGCACGGGACGCATCACTGTCGCGTTCGCATGGAAAGTGACCTGGCAAGAAGCTCGGTGGTTTTATGGAGAACGCTACCATGGAGTCAGCTCATCGACAAGCTGAACGGGGTTTCACCCCTTAAATCGGTCTACGACATCTTGTGCGGACGCTATACACCGCATCAGATCTCTTGTTTGTCCCATCCGACGTGGGTAGGTTAATACATCCAGAGACTGATGGACTGTTACGGTCATCGGACCTGTTCTTAGATCCTTTGCCTCGTTTGGATACTGTGAGGTCGTTACTATATTTGAGGGTAGAGGGCGCGTAACGTCCAAGTGGTCATAAGACTTGATCTGTTCGTGCATGTCCGATGCGGACCCTATCTACTCTGAACTCTTCTTGAAGGTGTCTGCCTCGCGAAGTTATATCCGATTCAATAGGGGAGTTGCCGGACTTTTACGGGCACACAACCGCC");//Generator.GetRandomDNA(25546));
//dna2 = new DNA("CGCATAAGGTTTGTTCGGTTCCGATTGCCATCAACCGCCCGACCCGGACATCCCTTCGAGTAAGGAGAGCTGATTTAAATTCCTGGGTATTAGTTGAACCTTAGCGGGACGTTCAACGCTGCCTGGTCTATCTAAGTTAAACCCTACGAGAGTGAATCATAAGGTGGTACGCAAACGACTTTTCGTAAGGCTTGGACCCGAATCTGAAGTCGCCAATAATGCCCACGCGATCACAGCATGCAGTCCTAAGAGCGCTCTTCGACGAGCACTCTTTCCGAGCCTT");//Generator.GetRandomDNA(25546));
//res = DNA.Profil.Compare(dna1.GpsProfil, dna2.GpsProfil);

//Console.WriteLine(res);

//DNA dna10 = new DNA("GAATTCAGTTA");//Generator.GetRandomDNA(25546));
//DNA dna20 = new DNA("GGATCGA");//Generator.GetRandomDNA(25546));

//res = DNA.Profil.Compare(dna10.GpsProfil, dna20.GpsProfil);
//Console.WriteLine(res);
//var resNidelman0 = dna10.NeedlemanWunschResult(dna20.DnaString);
//Debug.WriteLine(string.Format("Needelman: {1} {0}CATAnalis: {2} {0}", Environment.NewLine, resNidelman0, res0));
////DNA dna1 = new DNA("GAATTCAGTTA");//Generator.GetRandomDNA(25546));
////DNA dna2 = new DNA("GGATCGA");//Generator.GetRandomDNA(25546));

////var res = DNA.Profil.Compare(dna1.GpsProfil, dna2.GpsProfil);
////var res1 = DNA.Profil.Compare(dna1.GpsProfil, dna1.GpsProfil);
////var res2 = DNA.Profil.Compare(dna2.GpsProfil, dna2.GpsProfil);
////this.chart1.ChartAreas[0].TransformPoints();

//Debug.WriteLine("END EXPERIMENT");//dna1.GpsProfil.ToString());