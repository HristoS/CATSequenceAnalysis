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

        public (double exactMatch, double gaps) NeedlemanWunschResult(string alignSeq)
        {
            return (NeedlemanWunsch.Calculate(this.DnaString, alignSeq), 0);
        }
    }

    public class CATProfile
    {
        private const double acDistance = 0.86d;
        private BenchmarkProfile a;
        private BenchmarkProfile t;
        private BenchmarkProfile c;
        private string dnaString;

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
            double bonus = 1; // 2 - Benchmark.maxPoint;
            double bonusTotalA = 0;
            double bonusTotalT = 0;
            double bonusTotalC = 0;
            for (int i = 0; i < dnaString.Length; i++)
            {
                var nearMatchA = benchmarkA.NearMatch(i, dnaString[i]);
                var nearMatchT = benchmarkT.NearMatch(i, dnaString[i]);

                var exacatMatchA = benchmarkA.ExactMatch(i, dnaString[i]);
                var exacatMatchT = benchmarkT.ExactMatch(i, dnaString[i]);

                var exacatMatchC = benchmarkC.ExactMatch(i, dnaString[i]);//(nearMatchA + nearMatchT);//
                var nearMatchC = benchmarkC.NearMatch(i, dnaString[i]);//(exacatMatchA + exacatMatchT) / 2; //

                nearMatchesA += (nearMatchA + prevMatchesA * (nearMatchA - 0d));
                nearMatchesT += (nearMatchT + prevMatchesT * (nearMatchT - 0d));
                nearMatchesC += (nearMatchC + prevMatchesC * (nearMatchC - 0d));

                exactMatchesA += (exacatMatchA + prevMatchesA * exacatMatchC);
                exactMatchesT += (exacatMatchT + prevMatchesT * exacatMatchT);
                exactMatchesC += (exacatMatchC + prevMatchesC * exacatMatchC);

                prevMatchesA = (exacatMatchA + nearMatchA - Benchmark.minPoint) / Benchmark.maxPoint; // [bonus, 0]
                prevMatchesT = (exacatMatchT + nearMatchT - Benchmark.minPoint) / Benchmark.maxPoint; // [bonus, 0]
                prevMatchesC = (exacatMatchC + nearMatchC - Benchmark.minPoint) / Benchmark.maxPoint; // [bonus, 0]

                bonusTotalA += prevMatchesA;
                bonusTotalT += prevMatchesT;
                bonusTotalC += prevMatchesC;

                if (i == dnaString.Length - 1)
                {
                    nearMatchesA += (nearMatchA + prevMatchesA * nearMatchA);
                    nearMatchesT += (nearMatchT + prevMatchesT * nearMatchT);
                    nearMatchesC += (nearMatchC + prevMatchesC * nearMatchC);

                    exactMatchesA += (exacatMatchA + prevMatchesA * exacatMatchC);
                    exactMatchesT += (exacatMatchT + prevMatchesT * exacatMatchT);
                    exactMatchesC += (exacatMatchC + prevMatchesC * exacatMatchC);
                }
            }
            var bonusTotal = (bonusTotalA + bonusTotalT + bonusTotalC) / 3d;

            dnaDistanceA = (nearMatchesA + exactMatchesA) / (double)(bonusTotal + dnaString.Length);
            dnaDistanceT = (nearMatchesT + exactMatchesT) / (double)(bonusTotal + dnaString.Length);
            dnaDistanceC = (nearMatchesC + exactMatchesC) / (double)(bonusTotal + dnaString.Length);

            a = new BenchmarkProfile(dnaDistanceA, 1, dnaString.Length);
            t = new BenchmarkProfile(dnaDistanceT, 1, dnaString.Length);
            c = new BenchmarkProfile(dnaDistanceC, acDistance, dnaString.Length);
            if (CanContribute() == false)
            {
                bonusTotal.ToString();
            }
            a.Calculate(t);
            t.Calculate(a);
            c.Calculate(a);
        }

        public static double CompareFixed(CATProfile x, CATProfile y)
        {
            var resA = Math.Sqrt(Math.Pow(x.a.D - y.a.D, 2) + Math.Pow(x.a.H - y.a.H, 2));
            var resT = Math.Sqrt(Math.Pow(x.t.D - y.t.D, 2) + Math.Pow(x.t.H - y.t.H, 2));

            return ((resA + resT) / 2d);
        }

        public static double Compare(CATProfile x, CATProfile y)
        {
            var resA = Math.Sqrt(Math.Pow(x.a.D - y.a.D, 2) + Math.Pow(x.a.H - y.a.H, 2));// + Math.Pow((x.a.DNADistance - y.a.DNADistance), 2));
            var resT = Math.Sqrt(Math.Pow(x.t.D - y.t.D, 2) + Math.Pow(x.t.H - y.t.H, 2));// + Math.Pow((x.t.DNADistance - y.t.DNADistance), 2));
            var resC = Math.Sqrt(Math.Pow(x.c.D - y.c.D, 2) + Math.Pow(x.c.H - y.c.H, 2));// + Math.Pow((x.c.DNADistance - y.c.DNADistance), 2));

            return 1 - ((resC + resA + resT) / 3d);
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

        private bool CanContribute(double distanceA, double distanceT)
        {
            return distanceA + distanceT > 1
               && distanceA + 1 > distanceT
               && distanceT + 1 > distanceA;
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
        public const double minPoint = 0.36d;
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