using Clipper2Lib;
using System;

namespace SequenceAnalyses
{
    public enum Base
    { A, T, C, G };

    public class DNA
    {
        //private static NeedlemanWunsch NeedlemanWunschCalculation = new NeedlemanWunsch();

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
        private const double acDistance = 0.75d; // 1 full match on each of the other benchmarks
        private BenchmarkProfile a;
        private BenchmarkProfile t;
        private BenchmarkProfile c;
        private string dnaString;
        private readonly PathsD path = new PathsD();
        private readonly double area;

        public string DnaString => dnaString;

        public CATProfile(string dnaString)
            : this(BenchmarkRepository.A, BenchmarkRepository.T, BenchmarkRepository.C, dnaString)
        {
        }

        public CATProfile(Benchmark benchmarkA, Benchmark benchmarkT, Benchmark benchmarkC, string dnaString)
        {
            this.dnaString = dnaString;
            double nearMatchesA = 0d;
            double nearMatchesT = 0d;
            double nearMatchesC = 0d;

            double dnaDistanceA = 0d;
            double dnaDistanceT = 0d;
            double dnaDistanceC = 0d;

            double exactMatchesA = 0d;
            double exactMatchesT = 0d;
            double exactMatchesC = 0d;

            double prevMatchesA = 0d;
            double prevMatchesT = 0d;
            double prevMatchesC = 0d;

            double bonusTotalA = 0;
            double bonusTotalT = 0;
            double bonusTotalC = 0;

            double sequencePrevLengthC = 0;
            double sequencePrevLengthA = 0;
            double sequencePrevLengthT = 0;

            for (int i = 0; i < dnaString.Length; i++)
            {
                var nearMatchC = benchmarkC.NearMatch(i, dnaString[i]);
                var nearMatchA = benchmarkA.NearMatch(i, dnaString[i]);
                var nearMatchT = benchmarkT.NearMatch(i, dnaString[i]);

                var exactMatchC = benchmarkC.ExactMatch(i, dnaString[i]);
                var exactMatchA = benchmarkA.ExactMatch(i, dnaString[i]);
                var exactMatchT = benchmarkT.ExactMatch(i, dnaString[i]);

                nearMatchesC += nearMatchC + prevMatchesC * nearMatchC;
                nearMatchesA += nearMatchA + prevMatchesA * nearMatchA;
                nearMatchesT += nearMatchT + prevMatchesT * nearMatchT;

                exactMatchesC += exactMatchC + prevMatchesC * exactMatchC;
                exactMatchesA += exactMatchA + prevMatchesA * exactMatchA;
                exactMatchesT += exactMatchT + prevMatchesT * exactMatchT;

                prevMatchesC = sequencePrevLengthC / (i + 1) + exactMatchC + nearMatchC - Benchmark.minPoint;
                prevMatchesA = sequencePrevLengthA / (i + 1) + exactMatchA + nearMatchA - Benchmark.minPoint;
                prevMatchesT = sequencePrevLengthT / (i + 1) + exactMatchT + nearMatchT - Benchmark.minPoint;

                bonusTotalC += i == dnaString.Length - 1 ? 0 : prevMatchesC;
                bonusTotalA += i == dnaString.Length - 1 ? 0 : prevMatchesA;
                bonusTotalT += i == dnaString.Length - 1 ? 0 : prevMatchesT;

                sequencePrevLengthC += prevMatchesC;
                sequencePrevLengthA += prevMatchesA;
                sequencePrevLengthT += prevMatchesT;
            }

            dnaDistanceC = (nearMatchesC + exactMatchesC) / (double)(bonusTotalC + dnaString.Length);
            dnaDistanceA = (nearMatchesA + exactMatchesA) / (double)(bonusTotalA + dnaString.Length);
            dnaDistanceT = (nearMatchesT + exactMatchesT) / (double)(bonusTotalT + dnaString.Length);

            a = new BenchmarkProfile(dnaDistanceA, 1, dnaString.Length);
            t = new BenchmarkProfile(dnaDistanceT, 1, dnaString.Length);
            c = new BenchmarkProfile(dnaDistanceC, acDistance, dnaString.Length);

            a.Calculate(t);
            t.Calculate(a);
            c.Calculate(a);

            path.Add(Clipper.MakePath(new double[] { c.D * 100, c.H * 100, a.D * 100, a.H * 100, t.D * 100, t.H * 100 }));
            area = Math.Abs(Clipper.Area(path));
        }

        public static double CompareAreas(CATProfile x, CATProfile y)
        {
            var intersect = Clipper.Intersect(x.path, y.path, FillRule.NonZero, 8);
            return Clipper.Area(intersect) / Math.Min(x.area, y.area);

            var areaX = Math.Abs(x.c.D * (x.a.H - x.t.H) + x.a.D * (x.t.H - x.c.H) + x.t.D * (x.c.H - x.a.H)) / 2;
            var areaY = Math.Abs(y.c.D * (y.a.H - y.t.H) + y.a.D * (y.t.H - y.c.H) + y.t.D * (y.c.H - y.a.H)) / 2;

            return Math.Min(x.dnaString.Length, y.dnaString.Length) / Math.Max(x.dnaString.Length, y.dnaString.Length) + Math.Min(areaX, areaY) / Math.Max(areaX, areaY);
        }

        public static double CompareDistances(CATProfile x, CATProfile y)
        {
            var resC = Math.Sqrt(Math.Pow(x.c.D - y.c.D, 2) + Math.Pow(x.c.H - y.c.H, 2));// + Math.Pow((x.c.DNADistance - y.c.DNADistance), 2));
            var resA = Math.Sqrt(Math.Pow(x.a.D - y.a.D, 2) + Math.Pow(x.a.H - y.a.H, 2));// + Math.Pow((x.a.DNADistance - y.a.DNADistance), 2));
            var resT = Math.Sqrt(Math.Pow(x.t.D - y.t.D, 2) + Math.Pow(x.t.H - y.t.H, 2));// + Math.Pow((x.t.DNADistance - y.t.DNADistance), 2));

            //var resA = Math.Sqrt(Math.Pow(x.a.D - y.a.D, 2) + Math.Pow((x.a.DNADistance - y.a.DNADistance), 2));
            //var resT = Math.Sqrt(Math.Pow(x.t.D - y.t.D, 2) + Math.Pow((x.t.DNADistance - y.t.DNADistance), 2));
            //var resC = Math.Sqrt(Math.Pow(x.c.D - y.c.D, 2) + Math.Pow((x.c.DNADistance - y.c.DNADistance), 2));

            return 1 - (resC + resA + resT) / 3d;
        }

        public static double Compare(CATProfile x, CATProfile y)
        {
            //return CompareAreas(x, y);
            return CompareDistances(x, y);
        }

        public bool AreEqual(CATProfile profile)
        {
            return this.a.AreEqual(profile.a) && this.t.AreEqual(profile.t) && this.c.AreEqual(profile.c);
        }

        public bool CanContribute()
        {
            return a.DNADistance + t.DNADistance > 1
                && a.DNADistance + 1 > t.DNADistance
                && t.DNADistance + 1 > a.DNADistance
                && a.DNADistance + c.DNADistance > acDistance
                && a.DNADistance + acDistance > c.DNADistance
                && c.DNADistance + acDistance > a.DNADistance;
        }
    }

    public class BenchmarkProfile
    {
        private double dnaDistance = 0;
        private double dnaMatches = 0;
        private double benchmarkDistance = 0;
        private double dnaLength = 0;

        private double cos;
        private double h;
        private double d;

        public double DNAMatches { get => dnaMatches; }

        public double DNADistance { get => dnaDistance; }

        public double DNALength { get => dnaLength; }

        public double Cos { get => this.cos; }

        public double H { get => this.h; }

        public double D { get => this.d; }

        public BenchmarkProfile(double dnaDistance, double benchmarkDistance)
        {
            this.dnaDistance = dnaDistance;
            this.benchmarkDistance = benchmarkDistance;
        }

        public BenchmarkProfile(double dnaDistance, double benchmarkDistance, double dnaLength)
        {
            this.dnaDistance = dnaDistance;
            this.benchmarkDistance = benchmarkDistance;
            this.dnaLength = dnaLength;
        }

        public void Calculate(BenchmarkProfile profile)
        {
            this.cos = BenchmarkCOS(profile);
            this.h = BenchmarkH(this.cos);
            this.d = BenchmarkD(this.cos);
        }

        private double BenchmarkCOS(BenchmarkProfile profile)
        {
            return (Math.Pow(this.dnaDistance, 2) + Math.Pow(benchmarkDistance, 2) - Math.Pow(profile.dnaDistance, 2)) / (2 * this.dnaDistance * benchmarkDistance);
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
        public const double maxPoint = 1d;
        public const double minPoint = 0.4d;
        private double[] baseDistance = new double[4] { 0d, 0.6d, minPoint, 0.6d };
        private string baseSequence;

        public Benchmark(string sequence)
        {
            baseSequence = sequence;
        }

        public char this[int index]
        {
            get { return this.baseSequence[index % baseSequence.Length]; }
        }

        public double NearMatch(int index, char @base)
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

        public double ExactMatch(int index, char @base)
        {
            return this.baseSequence[index % baseSequence.Length] == @base ? maxPoint : 0;
        }

        //public double Accumulate(int index, char @base)
    }

    public static class BenchmarkRepository
    {
        private static Benchmark a = new Benchmark("ACGT");
        private static Benchmark t = new Benchmark("GTAC");
        private static Benchmark c = new Benchmark("AGTC");//"CTGA"
        public static Benchmark A => a;

        public static Benchmark T => t;

        public static Benchmark C => c;
    }
}