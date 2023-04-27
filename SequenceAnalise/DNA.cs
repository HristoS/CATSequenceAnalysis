using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Linq;
using static SequenceAnalyses.BenchmarkProfile;

namespace SequenceAnalyses
{
    public enum Base
    { A, T, C, G };

    public class DNA
    {
        public CATProfile CATProfile { get; }

        public string DnaString { get; }

        public DNA(string dnaString)
        {
            this.DnaString = dnaString;

            this.CATProfile = new CATProfile(dnaString);
        }
    }

    public class CATProfile
    {
        public class CATProfileSlim
        {
            public BenchmarkProfilSlim c;
            public BenchmarkProfilSlim a;
            public BenchmarkProfilSlim t;

            public CATProfile GetCATProfile(string dnsString)
            {
                var result = new CATProfile();
                result.a = a.GetBenchmarkProfilе(BenchmarkRepository.A);
                result.c = c.GetBenchmarkProfilе(BenchmarkRepository.C);
                result.t = t.GetBenchmarkProfilе(BenchmarkRepository.T);
                result.DnaString = dnsString;
                return result;
            }

            public static int GetSize()
            {
                return BenchmarkProfilSlim.GetSize() * 3;
            }

            public byte[] Serialize()
            {
                return c.Serialize()
                    .Concat(a.Serialize())
                    .Concat(t.Serialize())
                    .ToArray();
            }

            public int Deserialize(byte[] bytes)
            {
                var index = 0;
                this.c = new BenchmarkProfilSlim();
                index += this.c.Deserialize(bytes, index);
                this.a = new BenchmarkProfilSlim();
                index += this.a.Deserialize(bytes, index);
                this.t = new BenchmarkProfilSlim();
                index += this.t.Deserialize(bytes, index);
                return index;
            }

            public bool AreEqual(CATProfileSlim profileSlim)
            {
                return this.c.AreEqual(profileSlim.c) && this.a.AreEqual(profileSlim.a) && this.t.AreEqual(profileSlim.t);
            }

            public override string ToString()
            {
                return $"C.H: {c.h} C.D: {c.d} A.H: {a.h} A.D: {a.d} T.H: {t.h} T.D: {t.d}";
            }
        }

        public string DnaString { get; private set; }

        private BenchmarkProfile c;
        private BenchmarkProfile a;
        private BenchmarkProfile t;

        public CATProfile(string dnaString)
        {
            this.DnaString = dnaString;
            a = new BenchmarkProfile(BenchmarkRepository.A, dnaString);
            c = new BenchmarkProfile(BenchmarkRepository.C, dnaString);
            t = new BenchmarkProfile(BenchmarkRepository.T, dnaString);

            a.Calculate(t);
            t.Calculate(a);
            c.Calculate(a);
        }

        private CATProfile()
        {
        }

        public CATProfileSlim GetCATProfilSlim()
        {
            return new CATProfileSlim()
            {
                a = this.a.GetBenchmarkProfilSlim(),
                t = this.t.GetBenchmarkProfilSlim(),
                c = this.c.GetBenchmarkProfilSlim()
            };
        }

        public static double CompareFixed(CATProfile x, CATProfile y)
        {
            CATProfile longer = x.DnaString.Length > y.DnaString.Length ? x : y;
            CATProfile shorter = longer == x ? y : x;

            var delta = shorter.DnaString.Length / (double)longer.DnaString.Length;

            var resC = 1 - Math.Sqrt(Math.Pow(delta * longer.c.D - shorter.c.D, 2) + Math.Pow(delta * longer.c.H - shorter.c.H, 2));
            var resA = 1 - Math.Sqrt(Math.Pow(delta * longer.a.D - shorter.a.D, 2) + Math.Pow(delta * longer.a.H - shorter.a.H, 2));
            var resT = 1 - Math.Sqrt(Math.Pow(delta * longer.t.D - shorter.t.D, 2) + Math.Pow(delta * longer.t.H - shorter.t.H, 2));

            return ((resC + resA + resT) / 3d);
        }

        public static double Compare(CATProfile x, CATProfile y)
        {
            var resC = 1 - Math.Sqrt(Math.Pow(x.c.D - y.c.D, 2) + Math.Pow(x.c.H - y.c.H, 2));
            var resA = 1 - Math.Sqrt(Math.Pow(x.a.D - y.a.D, 2) + Math.Pow(x.a.H - y.a.H, 2));
            var resT = 1 - Math.Sqrt(Math.Pow(x.t.D - y.t.D, 2) + Math.Pow(x.t.H - y.t.H, 2));

            return ((resC + resA + resT) / 3d);
        }

        public bool AreEqual(CATProfile profile)
        {
            return this.a.AreEqual(profile.a) && this.t.AreEqual(profile.t) && this.c.AreEqual(profile.c);
        }
    }

    public class BenchmarkProfile
    {
        [Serializable]
        public class BenchmarkProfilSlim
        {
            public double h;

            public double d;

            public BenchmarkProfile GetBenchmarkProfilе(Benchmark benchmark)
            {
                var result = new BenchmarkProfile();

                result.h = h;
                result.d = d;
                result.benchmark = benchmark;
                return result;
            }

            public byte[] Serialize()
            {
                return BitConverter.GetBytes(h)
                    .Concat(BitConverter.GetBytes(d))
                    .ToArray();
            }

            public int Deserialize(byte[] bytes, int index = 0)
            {
                this.h = BitConverter.ToDouble(bytes, index);
                this.d = BitConverter.ToDouble(bytes, index + 8);
                return 16;
            }

            public static int GetSize()
            {
                return 16;
            }

            public bool AreEqual(BenchmarkProfilSlim profilSlim)
            {
                return this.h == profilSlim.h && this.d == profilSlim.d;
            }
        }

        private double dnaDistance = 0;
        private double dnaMatches = 0;
        private Benchmark benchmark;
        private string dnaString;
        private double cos;
        private double h;
        private double d;

        public double DNAMatches { get => dnaMatches; }

        public double DNADistance { get => dnaDistance; }

        public double Cos { get => this.cos; }

        public double H { get => this.h; }

        public double D { get => this.d; }

        public BenchmarkProfile(Benchmark benchmark, string dnaString)
        {
            this.benchmark = benchmark;
            this.dnaString = dnaString;

            for (int i = 0; i < dnaString.Length; i++)
            {
                dnaMatches += benchmark.Match(i, dnaString[i]);
            }

            this.dnaDistance = 1 - dnaMatches / (double)dnaString.Length;
        }

        private BenchmarkProfile()
        {
        }

        public BenchmarkProfilSlim GetBenchmarkProfilSlim()
        {
            return new BenchmarkProfilSlim()
            {
                d = this.D,
                h = this.h
            };
        }

        public void Calculate(BenchmarkProfile profile)
        {
            this.cos = BenchmarkCOS(profile);
            this.h = BenchmarkH(this.cos);
            this.d = BenchmarkD(this.cos);
        }

        private double BenchmarkDistance(Benchmark benchmark)
        {
            var quotient = dnaString.Length % 4;
            var remainder = dnaString.Length / 4;

            double matches = 0;
            for (int i = 0; i < 4; i++)
            {
                matches += this.benchmark.Match(i, benchmark[i]);
            }
            matches *= quotient;
            for (int i = 0; i < remainder; i++)
            {
                matches += this.benchmark.Match(i, benchmark[i]);
            }

            return 1 - matches / (double)Math.Max(dnaString.Length, 4);
        }

        private double BenchmarkCOS(BenchmarkProfile benchmarkProfile)
        {
            var benchmarkDistance = BenchmarkDistance(benchmarkProfile.benchmark);
            if (this.dnaDistance * benchmarkDistance == 0)
                return 0;
            return (Math.Pow(this.dnaDistance, 2) + Math.Pow(benchmarkDistance, 2) - Math.Pow(benchmarkProfile.dnaDistance, 2)) / (2 * this.dnaDistance * benchmarkDistance);
        }

        private double BenchmarkD(double benchmarkCOS)
        {
            return this.dnaDistance * benchmarkCOS;
        }

        private double BenchmarkH(double benchmarkCOS)
        {
            return Math.Sqrt(Math.Pow(this.dnaDistance, 2) - Math.Pow(this.dnaDistance * benchmarkCOS, 2));
        }

        public bool AreEqual(BenchmarkProfile profil)
        {
            return this.h == profil.h && this.d == profil.d;
        }
    }

    public class Benchmark
    {
        private double[] baseDistance = new double[4] { 1d, 0d, 0d, 0d };
        private string baseSequence;

        public Benchmark(string sequence)
        {
            baseSequence = sequence;
        }

        public char this[int index]
        {
            get { return this.baseSequence[index % baseSequence.Length]; }
        }

        public double Match(int index, char @base)
        {
            //return this.baseSequence[index % baseSequence.Length] == @base ? 1 : 0;
            if (this.baseSequence[(index + 0) % baseSequence.Length] == @base)
                return baseDistance[0];
            if (this.baseSequence[(index + 1) % baseSequence.Length] == @base)
                return baseDistance[1];
            if (this.baseSequence[(index + 2) % baseSequence.Length] == @base)
                return baseDistance[2];
            if (this.baseSequence[(index + 3) % baseSequence.Length] == @base)
                return baseDistance[3];

            return 0d;
        }
    }

    public static class BenchmarkRepository
    {
        private static Benchmark a = new Benchmark("ACGT");
        private static Benchmark t = new Benchmark("TGCA");
        private static Benchmark c = new Benchmark("CGAT");
        private static Benchmark g = new Benchmark("GCTA");

        public static Benchmark C => c;

        public static Benchmark A => a;

        public static Benchmark T => t;

        public static Benchmark G => g;
    }
}