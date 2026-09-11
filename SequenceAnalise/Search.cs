using Accord.Math;
using MathNet.Numerics;
using SequenceAnalise;
using System;
using System.Collections.Concurrent;
using System.Diagnostics;

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

        // Inverse of the value computation in GetPermutationsQ: decodes a base-4 index
        // back into the sequence it corresponds to (position 0 is the least significant digit).
        private static string ValueToSequence(double value, int sequenceLength)
        {
            var values = new char[] { 'A', 'C', 'G', 'T' };
            var sequence = new char[sequenceLength];
            for (int i = 0; i < sequenceLength; i++)
            {
                sequence[i] = values[(int)(value % values.Length)];
                value = Math.Floor(value / values.Length);
            }

            return new string(sequence);
        }

        public static Task<double> IterateInParallel(int sequenceLength, Action<string> iterator)
        {
            return IterateInParallel(sequenceLength, Environment.ProcessorCount, iterator);
        }

        public static async Task<double> IterateInParallel(int sequenceLength, int degreeOfParallelism, Action<string> iterator)
        {
            double totalPermutations = Math.Pow(4, sequenceLength);
            int shardCount = (int)Math.Max(1, Math.Min(degreeOfParallelism, totalPermutations));

            double baseSize = Math.Floor(totalPermutations / shardCount);
            double remainder = totalPermutations - baseSize * shardCount;

            var tasks = new Task<double>[shardCount];
            double cursor = 0;
            for (int shard = 0; shard < shardCount; shard++)
            {
                double shardSize = baseSize + (shard < remainder ? 1 : 0);
                var startSequence = ValueToSequence(cursor, sequenceLength);
                cursor += shardSize;
                var endSequence = ValueToSequence(cursor - 1, sequenceLength);

                tasks[shard] = Task.Run<double>(() =>
                {
                    // GetPermutationsQ requires a strictly wider start/end range, so a
                    // single-sequence shard (only possible when totalPermutations is small
                    // relative to the requested parallelism) has to be handled directly.
                    if (startSequence == endSequence)
                    {
                        iterator(startSequence);
                        return 1d;
                    }

                    var totalRecords = 0d;

                    foreach (var permutation in Search.GetPermutationsQ(startSequence, endSequence, true))
                    {
                        totalRecords++;
                        iterator(permutation);
                    }

                    return totalRecords;
                });
            }

            return (await Task.WhenAll(tasks)).Sum();
        }

        /// <summary>
        /// Generates the collision report.
        /// Generates sequenceCount random DNA sequens with length sequenceLength. Finds all collisions from all possible permutation for the generated sequence.
        /// Apples NeedlemanWunsch aliment on sequence which are in collision.
        /// </summary>
        /// <param name="sequenceCount">The sequence count.</param>
        /// <param name="sequenceLength">The sequence length.</param>
        /// <param name="seed">Seed for the random reference sequences, so the report is reproducible.</param>
        /// <param name="result">The result file.</param>
        public static async Task GenerateCollisionReport(int sequenceCount, int sequenceLength, int seed, StreamWriter result)
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
                var str = Generator.GetRandomDNA(sequenceLength, new Random(HashCode.Combine(seed, i)));
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

            double catComparison = CATProfile.Compare(dna.CATProfile, dna2.CATProfile);
            sw.Stop();
            double catMs = sw.Elapsed.TotalMilliseconds;
            sw.Restart();

            var nwComparison = NeedlemanWunsch.Calculate(dna.DnaString, dna2.DnaString);
            sw.Stop();
            double nwMs = sw.Elapsed.TotalMilliseconds;

            return (catComparison, catMs, nwComparison, nwMs);
        }

        public static (double CATComparison, double CATms, double NWComparison, double NWms) CompareCATWithKMP(DNA dna, DNA dna2)
        {
            Stopwatch sw = new Stopwatch();
            sw.Start();

            double catComparison = CATProfile.Compare(dna.CATProfile, dna2.CATProfile);
            sw.Stop();
            double catMs = sw.Elapsed.TotalMilliseconds;
            sw.Restart();

            var nwComparison = PatternSearchingKMP.KMPSearch(dna2.DnaString, dna.DnaString);
            sw.Stop();
            double nwMs = sw.Elapsed.TotalMilliseconds;

            return (catComparison, catMs, nwComparison.Count, nwMs);
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