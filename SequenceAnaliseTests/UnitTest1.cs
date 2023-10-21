using SequenceAnalyses;
using System;
using System.Collections.Concurrent;

namespace SequenceAnaliseTests
{
    public class Tests
    {
        public static IEnumerable<int> SequenceLengthProvider()
        {
            yield return 10;
            yield return 11;
            yield return 12;
            yield return 13;
            yield return 15;
        }

        public static IEnumerable<(int count, int secuenceLength)> CompareCATwithNWProvider()
        {
            yield return (50, 100);
            yield return (50, 1_000);
            yield return (50, 10_000);
            yield return (1, 50_000);
        }

        public static IEnumerable<(int count, int secuenceLength)> CompareCATwithNWLongProvider()
        {
            for (int i = 0; i < 50; i++)
            {
                yield return (i, 50_000);
            }
        }

        [SetUp]
        public void Setup()
        {
        }

        [Test]
        public void SecialCase()
        {
            //var dna = "TACACGCACG";
            var dna = "CATTTGGGAAAGTSC";
            var catProfile = new CATProfile(dna);

            CATProfile.Compare(catProfile, catProfile);
        }

        [Test, TestCaseSource(nameof(SequenceLengthProvider))]
        public async Task CanContrinuteForFormigTringle(int sequenceLength)
        {
            ConcurrentDictionary<string, bool> cannotContrubute = new ConcurrentDictionary<string, bool>();

            var totalRecords = await Search.IterateInParalel(sequenceLength, (permutation) =>
            {
                var profile = new CATProfile(permutation);
                if (profile.CanContribute() == false)
                {
                    cannotContrubute.AddOrUpdate(permutation, (key) => false, (key, val) => false);
                }
            });

            Assert.That(cannotContrubute.Count, Is.EqualTo(0));
        }

        [Test, TestCaseSource(nameof(SequenceLengthProvider))]
        public async Task ColisionReportFor_1000(int sequenceLength)
        {
            using var result = File.CreateText($"SearchResults{sequenceLength}.csv");
            await Search.GenerateColisionReport(100, sequenceLength, result);
            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(CompareCATwithNWProvider))]
        public async Task CompareCATwithNW((int count, int secuenceLength) input)
        {
            var ssl50 = (int)(input.secuenceLength * 0.5);
            var ssl70 = (int)(input.secuenceLength * 0.7);
            var ssl90 = (int)(input.secuenceLength * 0.9);
            var ssl97 = (int)(input.secuenceLength * 0.97);

            var halfs = new List<int>() { ssl50, ssl70, ssl90, ssl97 };

            var randomSecuenceLength = new List<int> { ssl50, ssl70, ssl90, ssl97, input.secuenceLength };

            Task<List<string>>[] experiments = new Task<List<string>>[input.count];
            for (int i = 0; i < input.count; i++)
            {
                experiments[i] = Task.Factory.StartNew(() =>
                {
                    List<string> result = new List<string>();

                    var dna = new DNA(Generator.GetRandomDNA(input.secuenceLength));
                    var exactSame = Search.ExactSame(dna);
                    result.Add($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Search.Format(Search.CompareCATWithNW(dna, exactSame))}");

                    foreach (var half in halfs)
                    {
                        var firstHalf = Search.SubsequenceFirstHalf(dna, half);
                        var secondHalf = Search.SubsequenceSecondHalf(dna, half);
                        var middle = Search.SubsequenceInTheMiddle(dna, half);

                        result.Add($"{dna.DnaString.Length}, {half}, FirstHalf, {Search.Format(Search.CompareCATWithNW(dna, firstHalf))}");
                        result.Add($"{dna.DnaString.Length}, {half}, SecondHalf, {Search.Format(Search.CompareCATWithNW(dna, secondHalf))}");
                        result.Add($"{dna.DnaString.Length}, {half}, Middle, {Search.Format(Search.CompareCATWithNW(dna, middle))}");
                    }

                    foreach (var rnd in randomSecuenceLength)
                    {
                        var dna2 = new DNA(Generator.GetRandomDNA(rnd));
                        var aproximate = Search.CompareCATWithNW(dna, dna2);
                        result.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Search.Format(aproximate)}");
                    }

                    return result;
                });
            }

            await Task.WhenAll(experiments);
            var results = experiments.SelectMany(x => x.Result).ToList();

            using var report = File.CreateText($"CompareCATWithNW_{input.secuenceLength}_{input.count}.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, NW Result, NW ms");
            foreach (var item in results)
            {
                report.WriteLine(item);
            }
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(CompareCATwithNWLongProvider))]
        public void CompareCATwithNWLong((int experiment, int secuenceLength) input)
        {
            var ssl50 = (int)(input.secuenceLength * 0.5);
            var ssl70 = (int)(input.secuenceLength * 0.7);
            var ssl90 = (int)(input.secuenceLength * 0.9);
            var ssl97 = (int)(input.secuenceLength * 0.97);

            var halfs = new List<int>() { ssl50, ssl70, ssl90, ssl97 };

            var randomSecuenceLength = new List<int> { ssl50, ssl70, ssl90, ssl97, input.secuenceLength };

            List<string> results = new List<string>();

            var dna = new DNA(Generator.GetRandomDNA(input.secuenceLength));
            var exactSame = Search.ExactSame(dna);
            results.Add($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Search.Format(Search.CompareCATWithNW(dna, exactSame))}");

            foreach (var half in halfs)
            {
                var firstHalf = Search.SubsequenceFirstHalf(dna, half);
                var secondHalf = Search.SubsequenceSecondHalf(dna, half);
                var middle = Search.SubsequenceInTheMiddle(dna, half);

                results.Add($"{dna.DnaString.Length}, {half}, FirstHalf, {Search.Format(Search.CompareCATWithNW(dna, firstHalf))}");
                results.Add($"{dna.DnaString.Length}, {half}, SecondHalf, {Search.Format(Search.CompareCATWithNW(dna, secondHalf))}");
                results.Add($"{dna.DnaString.Length}, {half}, Middle, {Search.Format(Search.CompareCATWithNW(dna, middle))}");
            }

            foreach (var rnd in randomSecuenceLength)
            {
                var dna2 = new DNA(Generator.GetRandomDNA(rnd));
                var aproximate = Search.CompareCATWithNW(dna, dna2);
                results.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Search.Format(aproximate)}");
            }

            using var report = File.CreateText($"CompareCATWithNW_{input.secuenceLength}_{input.experiment}.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, NW Result, NW ms");
            foreach (var item in results)
            {
                report.WriteLine(item);
            }
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(CompareCATwithNWProvider))]
        public async Task CompareCATWithKMP((int count, int secuenceLength) input)
        {
            var ssl50 = (int)(input.secuenceLength * 0.5);
            var ssl70 = (int)(input.secuenceLength * 0.7);
            var ssl90 = (int)(input.secuenceLength * 0.9);
            var ssl97 = (int)(input.secuenceLength * 0.97);

            var halfs = new List<int>() { ssl50, ssl70, ssl90, ssl97 };

            var randomSecuenceLength = new List<int> { ssl50, ssl70, ssl90, ssl97, input.secuenceLength };

            Task<List<string>>[] experiments = new Task<List<string>>[input.count];
            for (int i = 0; i < input.count; i++)
            {
                experiments[i] = Task.Factory.StartNew(() =>
                {
                    List<string> result = new List<string>();

                    var dna = new DNA(Generator.GetRandomDNA(input.secuenceLength));
                    var exactSame = Search.ExactSame(dna);
                    result.Add($"{dna.DnaString.Length}, {dna.DnaString.Length}, WithItself, {Search.Format(Search.CompareCATWithKMP(dna, exactSame))}");

                    foreach (var half in halfs)
                    {
                        var firstHalf = Search.SubsequenceFirstHalf(dna, half);
                        var secondHalf = Search.SubsequenceSecondHalf(dna, half);
                        var middle = Search.SubsequenceInTheMiddle(dna, half);

                        result.Add($"{dna.DnaString.Length}, {half}, FirstHalf, {Search.Format(Search.CompareCATWithKMP(dna, firstHalf))}");
                        result.Add($"{dna.DnaString.Length}, {half}, SecondHalf, {Search.Format(Search.CompareCATWithKMP(dna, secondHalf))}");
                        result.Add($"{dna.DnaString.Length}, {half}, Middle, {Search.Format(Search.CompareCATWithKMP(dna, middle))}");
                    }

                    foreach (var rnd in randomSecuenceLength)
                    {
                        var dna2 = new DNA(Generator.GetRandomDNA(rnd));
                        var aproximate = Search.CompareCATWithKMP(dna, dna2);
                        result.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Search.Format(aproximate)}");
                    }

                    return result;
                });
            }

            await Task.WhenAll(experiments);
            var results = experiments.SelectMany(x => x.Result).ToList();

            using var report = File.CreateText($"CompareCATWithKMP_{input.secuenceLength}_{input.count}.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, NW Result, NW ms");
            foreach (var item in results)
            {
                report.WriteLine(item);
            }
            report.Flush();

            Assert.Pass();
        }

        [Test]
        public async Task CombineCompareCATWithKMP()
        {
            using var report = File.CreateText($"ReportCATWithKMP.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, KMP Result, KMP ms");
            foreach (var file in Directory.GetFiles(Directory.GetCurrentDirectory(), "CompareCATWithKMP*.csv"))
            {
                string[] lines = File.ReadAllLines(file);
                foreach (var line in lines.Skip(1).ToArray())
                    report.WriteLine(line);

                report.Flush();
            }

            Assert.Pass();
        }

        [Test]
        public async Task CombineCompareCATWithNW()
        {
            using var report = File.CreateText($"ReportCATWithNW.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, NW Result, NW ms");
            foreach (var file in Directory.GetFiles(Directory.GetCurrentDirectory(), "CompareCATWithNW*.csv"))
            {
                string[] lines = File.ReadAllLines(file);
                foreach (var line in lines.Skip(1).ToArray())
                    report.WriteLine(line);

                report.Flush();
            }

            Assert.Pass();
        }
    }
}