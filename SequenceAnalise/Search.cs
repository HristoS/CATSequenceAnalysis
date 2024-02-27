using Accord.Math;
using MathNet.Numerics;
using SequenceAnalise;
using System;
using System.Collections.Concurrent;
using System.Diagnostics;
using System.Net;

namespace SequenceAnalyses
{
    public static class Search
    {
        public static IEnumerable<string> GetPermutations(int sequenceLength)
        {
            var values = new char[] { 'A', 'C', 'G', 'T' };
            // var values = new Dictionary<int, object[]>(parameters.Select((p, index) => new KeyValuePair<int, object[]>(index, p.GetCachedValue(handler).ToArray())).OrderBy(k => k.Key));
            //var lengths = values.Select(v => v.Value.Length);
            var permutations = Math.Pow(values.Length, sequenceLength);// values.Any() ? lengths.Aggregate(1, (acc, v) => acc * v) : 0;
            var counter = Enumerable.Repeat(0, sequenceLength).ToArray(); //  values.Select(v => 0).ToArray();

            for (double i = 0; i < permutations; i++)
            {
                //var parameterValues = new Dictionary<string, object>();
                var sequence = new char[sequenceLength];
                for (int j = 0; j < sequenceLength; j++)
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

        public static IEnumerable<string> GetPermutationsQ(string start, string end, bool include)
        {
            var values = new char[] { 'A', 'C', 'G', 'T' };

            if (start.Length != end.Length)
                throw new ArgumentException();

            var startValue = start.Select((x, i) => values.IndexOf(x) * Math.Pow(values.Length, i)).Sum(x => x);
            var endValue = end.Select((x, i) => values.IndexOf(x) * Math.Pow(values.Length, i)).Sum(x => x);

            if (startValue >= endValue)
                throw new ArgumentException();

            var sequenceLength = start.Length;

            var permutations = endValue - startValue;
            permutations = include ? permutations + 1 : permutations;
            var counter = start.Select(x => values.IndexOf(x)).ToArray();
            //var counterEnd = end.Select(x => values.IndexOf(x)).ToArray();

            for (double i = 0; i < permutations; i++)
            {
                var sequence = new char[sequenceLength];
                for (int j = 0; j < sequenceLength; j++)
                {
                    sequence[j] = values[counter[j]];
                }

                yield return new string(sequence);

                // update counter
                counter[0]++;
                for (int j = 0; j < counter.Length - 1; j++)
                {
                    if (counter[j] >= values.Length)
                    {
                        counter[j] = 0;
                        counter[j + 1]++;
                    }
                }
            }
        }

        public static IEnumerable<string> GetPermutationsAtoC(int sequenceLength, bool inclusive)
        {
            return Search.GetPermutationsQ(new string(Enumerable.Repeat('A', sequenceLength).ToArray()), new string(Enumerable.Repeat('C', sequenceLength).ToArray()), inclusive);
        }

        public static IEnumerable<string> GetPermutationsCtoG(int sequenceLength, bool inclusive)
        {
            return Search.GetPermutationsQ(new string(Enumerable.Repeat('C', sequenceLength).ToArray()), new string(Enumerable.Repeat('G', sequenceLength).ToArray()), inclusive);
        }

        public static IEnumerable<string> GetPermutationsGtoT(int sequenceLength, bool inclusive)
        {
            return Search.GetPermutationsQ(new string(Enumerable.Repeat('G', sequenceLength).ToArray()), new string(Enumerable.Repeat('T', sequenceLength).ToArray()), inclusive);
        }

        public static async Task<double> IterateInParallel(int sequenceLength, Action<string> iterator)
        {
            var taskAtoC = Task.Run<double>(() =>
            {
                var totalRecords = 0d;

                foreach (var permutation in Search.GetPermutationsAtoC(sequenceLength, false))
                {
                    totalRecords++;
                    iterator(permutation);
                }

                return totalRecords;
            });
            var taskCtoG = Task.Run<double>(() =>
            {
                var totalRecords = 0d;

                foreach (var permutation in Search.GetPermutationsCtoG(sequenceLength, false))
                {
                    totalRecords++;
                    iterator(permutation);
                }

                return totalRecords;
            });
            var taskGtoT = Task.Run<double>(() =>
            {
                var totalRecords = 0d;

                foreach (var permutation in Search.GetPermutationsGtoT(sequenceLength, true))
                {
                    totalRecords++;
                    iterator(permutation);
                }

                return totalRecords;
            });

            return (await Task.WhenAll(taskAtoC, taskCtoG, taskGtoT)).Sum();
        }

        /// <summary>
        /// Generates the collision report.
        /// Generates sequenceCount random DNA sequens with length sequenceLength. Finds all collisions from all possible permutation for the generated sequence.
        /// Apples NeedlemanWunsch aliment on sequence which are in collision.
        /// </summary>
        /// <param name="sequenceCount">The sequence count.</param>
        /// <param name="sequenceLength">The sequence length.</param>
        /// <param name="result">The result file.</param>
        public static async Task GenerateCollisionReport(int sequenceCount, int sequenceLength, StreamWriter result)
        {
            ConcurrentDictionary<string, double> results = new ConcurrentDictionary<string, double>();
            ConcurrentDictionary<string, ConcurrentDictionary<string, int>> labeledResults = new ConcurrentDictionary<string, ConcurrentDictionary<string, int>>();
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
            for (int i = 0; i < sequenceCount; i++)
            {
                var str = Generator.GetRandomDNA(sequenceLength);
                results[str] = 0;
                labeledResults[str] = new ConcurrentDictionary<string, int>(new Dictionary<string, int>()
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

            double totalRecords = await Search.IterateInParallel(sequenceLength, (permutation) =>
            {
                var profile = new CATProfile(permutation);
                foreach (var random in randomDNAProfiles)
                {
                    var catComparison = CATProfile.Compare(profile, random);
                    if (catComparison == 1d)
                    {
                        var res = NeedlemanWunsch.Calculate(random.DnaString, permutation);
                        var label = res.ToString("0.000").Substring(0, 3);
                        var labelResult = labeledResults.GetOrAdd(random.DnaString, (key) => new ConcurrentDictionary<string, int>());
                        labelResult.AddOrUpdate(label, (key) => 1, (key, value) => value + 1);
                    }
                }
            });

            result.WriteLine($"dna string ,totalRecords, cat matches, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0");

            var totals = new List<double>();
            foreach (var res in labeledResults)
            {
                double total = res.Value.Values.Sum();
                totals.Add(total);
                result.WriteLine($"{res.Key}, {totalRecords}, {total}, {string.Join(",", res.Value.OrderBy(x => x.Key).Select(x => x.Value * 100 / total))}");
                foreach (var label in labels)
                {
                    label.Value.Add(res.Value[label.Key] * 100 / total);
                }
                result.Flush();
            }

            result.WriteLine($"Average:, {100 * totals.Average() / totalRecords}, {totals.Average()}, {string.Join(",", labels.OrderBy(x => x.Key).Select(x => x.Value.Average()))}, {labels.OrderBy(x => x.Key).Select(x => x.Value.Average()).Sum()}");
        }

        public static (double CATComparison, double CATms, double NWComparison, double NWms) CompareCATWithNW(DNA dna, DNA dna2)
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            double CATComparison = CATProfile.Compare(dna.CATProfile, dna2.CATProfile);
            sw.Stop();
            double CATms = sw.Elapsed.TotalMilliseconds;
            sw.Restart();

            var NWComparison = NeedlemanWunsch.Calculate(dna.DnaString, dna2.DnaString);
            sw.Stop();
            double NWms = sw.Elapsed.TotalMilliseconds;

            return (CATComparison, CATms, NWComparison, NWms);
        }

        public static (double CATComparison, double CATms, double NWComparison, double NWms) CompareCATWithKMP(DNA dna, DNA dna2)
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            double CATComparison = CATProfile.Compare(dna.CATProfile, dna2.CATProfile);
            sw.Stop();
            double CATms = sw.Elapsed.TotalMilliseconds;
            sw.Restart();

            var NWComparison = PatternSearchingKMP.KMPSearch(dna2.DnaString, dna.DnaString);
            sw.Stop();
            double NWms = sw.Elapsed.TotalMilliseconds;

            return (CATComparison, CATms, NWComparison.Count, NWms);
        }

        public static string Format((double CATComparison, double CATms, double NWComparison, double NWms) result)
        {
            return $"{result.CATComparison}, {result.CATms}, {result.NWComparison}, {result.NWms}";
        }

        public static DNA ExactSame(DNA dna)
        {
            return dna;
        }

        public static DNA SubsequenceFirstHalf(DNA dna, double length)
        {
            return new DNA(dna.DnaString.Substring(0, (int)length));
        }

        public static DNA SubsequenceSecondHalf(DNA dna, double length)
        {
            return new DNA(dna.DnaString.Substring((int)(dna.DnaString.Length - length)));
        }

        public static DNA SubsequenceInTheMiddle(DNA dna, double length)
        {
            return new DNA(dna.DnaString.Substring((int)((dna.DnaString.Length - length) / 2), (int)length));
        }
    }
}