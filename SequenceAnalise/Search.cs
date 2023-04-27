using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using static SequenceAnalyses.CATProfile;

namespace SequenceAnalyses
{
    public static class Search
    {
        public static IEnumerable<string> GetPermutations(int sequenceLenght)
        {
            var values = new char[] { 'A', 'C', 'G', 'T' };
            // var values = new Dictionary<int, object[]>(parameters.Select((p, index) => new KeyValuePair<int, object[]>(index, p.GetCachedValue(handler).ToArray())).OrderBy(k => k.Key));
            //var lenghts = values.Select(v => v.Value.Length);
            var permutations = Math.Pow(values.Length, sequenceLenght);// values.Any() ? lenghts.Aggregate(1, (acc, v) => acc * v) : 0;
            var counter = Enumerable.Repeat(0, sequenceLenght).ToArray(); //  values.Select(v => 0).ToArray();

            for (double i = 0; i < permutations; i++)
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

        public static void GeneratePermutationData(int length)
        {
            var ascii = Encoding.ASCII;
            var bufferSize = length + CATProfileSlim.GetSize();
            using (var resultFile = File.Create($"Permutation_{length}.dat", bufferSize))
            {
                foreach (var permutation in GetPermutations(length))
                {
                    var dna = new DNA(permutation);

                    var cat = dna.CATProfile.GetCATProfilSlim();
                    var asciiStr = ascii.GetBytes(permutation);

                    resultFile.Write(asciiStr);
                    resultFile.Write(cat.Serialize());
                }

                resultFile.Flush();
            }
        }

        public static IEnumerable<CATProfile> GetProfile(int length)
        {
            var ascii = Encoding.ASCII;
            var bufferSize = length + CATProfileSlim.GetSize();
            var offset = 0;
            using (var resultFile = File.OpenRead($"Permutation_{length}.dat"))
            {
                resultFile.Seek(0, SeekOrigin.Begin);
                var buffer = new byte[bufferSize];
                while (resultFile.Read(buffer, offset, bufferSize) > 0)
                {
                    var dnaString = ascii.GetString(buffer, 0, length);
                    var cat = new CATProfileSlim();
                    cat.Deserialize(buffer.Skip(length).ToArray());
                    yield return cat.GetCATProfile(dnaString);
                }
            }
        }

        public static IEnumerable<CATProfile> GetProfile(int length, string file)
        {
            var ascii = Encoding.ASCII;
            var bufferSize = length + CATProfileSlim.GetSize();
            var offset = 0;
            using (var resultFile = File.OpenRead($"{file}"))
            {
                resultFile.Seek(0, SeekOrigin.Begin);
                var buffer = new byte[bufferSize];
                while (resultFile.Read(buffer, offset, bufferSize) > 0)
                {
                    var dnaString = ascii.GetString(buffer, 0, length);
                    var cat = new CATProfileSlim();
                    cat.Deserialize(buffer.Skip(length).ToArray());
                    yield return cat.GetCATProfile(dnaString);
                }
            }
        }

        public static long GetTotalRecords(int length)
        {
            FileInfo fi = new FileInfo($"Permutation_{length}.dat");
            var bufferSize = length + CATProfileSlim.GetSize();
            return fi.Length / bufferSize;
        }

        public static string GetDNAAt(int length, int index)
        {
            var ascii = Encoding.ASCII;
            var bufferSize = length + CATProfileSlim.GetSize();
            var offset = 0;
            using (var resultFile = File.OpenRead($"Permutation_{length}.dat"))
            {
                var buffer = new byte[bufferSize];
                resultFile.Seek(index * bufferSize, SeekOrigin.Begin);
                resultFile.Read(buffer, offset, bufferSize);
                return ascii.GetString(buffer, 0, length);
            }
        }

        public static void GenerateRandomSeq(int count, int sequenceLenght)
        {
            Console.WriteLine($"Generate {count} seq");
            var ascii = Encoding.ASCII;
            var bufferSize = sequenceLenght + CATProfileSlim.GetSize();
            using (var resultFile = File.Create($"{count}_Random_{sequenceLenght}.dat", bufferSize))
            {
                for (int i = 0; i < count; i++)
                {
                    var dnaStr = Generator.GetRandomDNA(sequenceLenght);
                    var dna = new DNA(dnaStr);

                    var cat = dna.CATProfile.GetCATProfilSlim();
                    var asciiStr = ascii.GetBytes(dnaStr);

                    resultFile.Write(asciiStr);
                    resultFile.Write(cat.Serialize());
                }

                resultFile.Flush();
            }
            Console.WriteLine($"END Generate {count} seq");
        }

        public static IEnumerable<CATProfile> GetProfile(int count, int sequenceLenght)
        {
            var ascii = Encoding.ASCII;
            var bufferSize = sequenceLenght + CATProfileSlim.GetSize();
            var offset = 0;
            using (var resultFile = File.OpenRead($"{count}_Random_{sequenceLenght}.dat"))
            {
                resultFile.Seek(0, SeekOrigin.Begin);
                var buffer = new byte[bufferSize];
                while (resultFile.Read(buffer, offset, bufferSize) > 0)
                {
                    var dnaString = ascii.GetString(buffer, 0, sequenceLenght);
                    var cat = new CATProfileSlim();
                    cat.Deserialize(buffer.Skip(sequenceLenght).ToArray());
                    yield return cat.GetCATProfile(dnaString);
                }
            }
        }
    }
}