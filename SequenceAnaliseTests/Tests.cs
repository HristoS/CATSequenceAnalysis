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

        public static IEnumerable<(int count, int sequenceLength)> CompareCATwithNWProvider()
        {
            yield return (50, 100);
            yield return (50, 1_000);
            yield return (50, 10_000);
            yield return (1, 50_000);
        }

        public static IEnumerable<(int count, int sequenceLength)> CompareCATwithNWLongProvider()
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

        [Test, TestCaseSource(nameof(SequenceLengthProvider))]
        public async Task CanContributeForFormingTriangle(int sequenceLength)
        {
            ConcurrentDictionary<string, bool> cannotContribute = new ConcurrentDictionary<string, bool>();

            var totalRecords = await Search.IterateInParallel(sequenceLength, (permutation) =>
            {
                var profile = new CATProfile(permutation);
                if (profile.CanContribute() == false)
                {
                    cannotContribute.AddOrUpdate(permutation, (key) => false, (key, val) => false);
                }
            });

            Assert.That(cannotContribute.Count, Is.EqualTo(0));
        }

        [Test, TestCaseSource(nameof(SequenceLengthProvider))]
        public async Task CollisionReportFor_1000(int sequenceLength)
        {
            using var result = File.CreateText($"SearchResults{sequenceLength}.csv");
            await Search.GenerateCollisionReport(1000, sequenceLength, result);
            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(CompareCATwithNWProvider))]
        public async Task CompareCATwithNW((int count, int sequenceLength) input)
        {
            var ssl50 = (int)(input.sequenceLength * 0.5);
            var ssl70 = (int)(input.sequenceLength * 0.7);
            var ssl90 = (int)(input.sequenceLength * 0.9);
            var ssl97 = (int)(input.sequenceLength * 0.97);

            var halfs = new List<int>() { ssl50, ssl70, ssl90, ssl97 };

            var randomSequenceLength = new List<int> { ssl50, ssl70, ssl90, ssl97, input.sequenceLength };

            Task<List<string>>[] experiments = new Task<List<string>>[input.count];
            for (int i = 0; i < input.count; i++)
            {
                experiments[i] = Task.Factory.StartNew(() =>
                {
                    List<string> result = new List<string>();

                    var dna = new DNA(Generator.GetRandomDNA(input.sequenceLength));
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

                    foreach (var rnd in randomSequenceLength)
                    {
                        var dna2 = new DNA(Generator.GetRandomDNA(rnd));
                        var approximate = Search.CompareCATWithNW(dna, dna2);
                        result.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Search.Format(approximate)}");
                    }

                    return result;
                });
            }

            await Task.WhenAll(experiments);
            var results = experiments.SelectMany(x => x.Result).ToList();

            using var report = File.CreateText($"CompareCATWithNW_{input.sequenceLength}_{input.count}.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, NW Result, NW ms");
            foreach (var item in results)
            {
                report.WriteLine(item);
            }
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(CompareCATwithNWLongProvider))]
        public void CompareCATwithNWLong((int experiment, int sequenceLength) input)
        {
            var ssl50 = (int)(input.sequenceLength * 0.5);
            var ssl70 = (int)(input.sequenceLength * 0.7);
            var ssl90 = (int)(input.sequenceLength * 0.9);
            var ssl97 = (int)(input.sequenceLength * 0.97);

            var halfs = new List<int>() { ssl50, ssl70, ssl90, ssl97 };

            var randomSequenceLength = new List<int> { ssl50, ssl70, ssl90, ssl97, input.sequenceLength };

            List<string> results = new List<string>();

            var dna = new DNA(Generator.GetRandomDNA(input.sequenceLength));
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

            foreach (var rnd in randomSequenceLength)
            {
                var dna2 = new DNA(Generator.GetRandomDNA(rnd));
                var approximate = Search.CompareCATWithNW(dna, dna2);
                results.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Search.Format(approximate)}");
            }

            using var report = File.CreateText($"CompareCATWithNW_{input.sequenceLength}_{input.experiment}.csv");
            report.WriteLine($"Dna1 Length, Dna2 Length, Kind, CAT Result, CAT ms, NW Result, NW ms");
            foreach (var item in results)
            {
                report.WriteLine(item);
            }
            report.Flush();

            Assert.Pass();
        }

        [Test, TestCaseSource(nameof(CompareCATwithNWProvider))]
        public async Task CompareCATWithKMP((int count, int sequenceLength) input)
        {
            var ssl50 = (int)(input.sequenceLength * 0.5);
            var ssl70 = (int)(input.sequenceLength * 0.7);
            var ssl90 = (int)(input.sequenceLength * 0.9);
            var ssl97 = (int)(input.sequenceLength * 0.97);

            var halfs = new List<int>() { ssl50, ssl70, ssl90, ssl97 };

            var randomSequenceLength = new List<int> { ssl50, ssl70, ssl90, ssl97, input.sequenceLength };

            Task<List<string>>[] experiments = new Task<List<string>>[input.count];
            for (int i = 0; i < input.count; i++)
            {
                experiments[i] = Task.Factory.StartNew(() =>
                {
                    List<string> result = new List<string>();

                    var dna = new DNA(Generator.GetRandomDNA(input.sequenceLength));
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

                    foreach (var rnd in randomSequenceLength)
                    {
                        var dna2 = new DNA(Generator.GetRandomDNA(rnd));
                        var approximate = Search.CompareCATWithKMP(dna, dna2);
                        result.Add($"{dna.DnaString.Length}, {dna2.DnaString.Length}, Random, {Search.Format(approximate)}");
                    }

                    return result;
                });
            }

            await Task.WhenAll(experiments);
            var results = experiments.SelectMany(x => x.Result).ToList();

            using var report = File.CreateText($"CompareCATWithKMP_{input.sequenceLength}_{input.count}.csv");
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