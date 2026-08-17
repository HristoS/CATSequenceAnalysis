using Clipper2Lib;
using System;

namespace SequenceAnalyses
{
    public enum Base
    { A, T, C, G };

    public class BucketDistributor
    {
        /// <summary>
        /// Distributes a total sum of 1.0 into a set number of buckets,
        /// where each previous bucket is multiplied by a specific growth factor.
        /// </summary>
        /// <param name="bucketCount">The number of buckets (e.g., 8, 16, 20, 32)</param>
        /// <param name="multiplier">The multiplier factor between adjacent buckets (e.g., 1.2, 1.5, 2.0)</param>
        /// <returns>An array of doubles representing the distributed values summing up to 1.0</returns>
        public static double[] DistributeValues(int bucketCount, double multiplier, double totalSum)
        {
            double[] buckets = new double[bucketCount];

            // Handle edge case where multiplier is 1 (equal distribution)
            if (Math.Abs(multiplier - 1.0) < 0.000001)
            {
                double evenShare = totalSum / bucketCount;
                for (int i = 0; i < bucketCount; i++)
                {
                    buckets[i] = evenShare;
                }
                return buckets;
            }

            // Calculate the value of the final, smallest bucket (Bucket N)
            // Formula: B_n = T * (r - 1) / (r^N - 1)
            double smallestBucket = totalSum * (multiplier - 1.0) / (Math.Pow(multiplier, bucketCount) - 1.0);

            // Fill the array from Bucket 1 (largest) down to Bucket N (smallest)
            for (int i = 0; i < bucketCount; i++)
            {
                int exponent = bucketCount - 1 - i;
                buckets[i] = smallestBucket * Math.Pow(multiplier, exponent);
            }

            return buckets;
        }
    }

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
        // Benchmark separations, derived from aligned-position matches between the benchmark
        // strings themselves (BenchmarkRepository.A="ACGT", T="GTAC", C="AGTC"): distance = 1 -
        // (shared positions / 4). A-T share 0 of 4 positions -> 1.0. A-C and T-C each share
        // exactly 1 of 4 (position 0 and position 3 respectively) -> 0.75. A-T-C therefore form a
        // fixed, valid triangle (0.75+0.75>1, 1+0.75>0.75), used below for true trilateration.
        private const double atDistance = 1d;

        private const double acDistance = 0.75d;
        private const double tcDistance = 0.75d;

        // Fixed reference positions for A, T, C in one shared plane: A at the origin, T on the
        // x-axis at distance atDistance, and C placed so |A-C|=acDistance and |T-C|=tcDistance.
        private static readonly double cRefX = (acDistance * acDistance - tcDistance * tcDistance + atDistance * atDistance) / (2 * atDistance);

        private static readonly double cRefY = Math.Sqrt(acDistance * acDistance - cRefX * cRefX);

        private BenchmarkProfile a;
        private BenchmarkProfile t;
        private BenchmarkProfile c;
        private string dnaString;

        // The DNA sequence's location, trilaterated from its three measured distances
        // (a.DNADistance, t.DNADistance, c.DNADistance) against the fixed A/T/C triangle above -
        // one point in a single shared frame, rather than three separately-framed projections.
        public (double X, double Y) Point { get; private set; }

        public string DnaString => dnaString;

        public CATProfile(string dnaString)
            : this(BenchmarkRepository.A, BenchmarkRepository.T, BenchmarkRepository.C, dnaString)
        {
        }

        public CATProfile(Benchmark benchmarkA, Benchmark benchmarkT, Benchmark benchmarkC, string dnaString)
        {
            this.dnaString = dnaString;
            double matchPointsC = 0d;
            double matchPointsA = 0d;
            double matchPointsT = 0d;

            for (int i = 0; i < dnaString.Length; i++)
            {
                var @base = dnaString[i];

                matchPointsC += benchmarkC.NearMatch(i, @base) + benchmarkC.ExactMatch(i, @base);
                matchPointsA += benchmarkA.NearMatch(i, @base) + benchmarkA.ExactMatch(i, @base);
                matchPointsT += benchmarkT.NearMatch(i, @base) + benchmarkT.ExactMatch(i, @base);
            }

            double dnaDistanceC = matchPointsC / (double)(dnaString.Length);
            double dnaDistanceA = matchPointsA / (double)(dnaString.Length);
            double dnaDistanceT = matchPointsT / (double)(dnaString.Length);

            a = new BenchmarkProfile(dnaDistanceA, atDistance, dnaString.Length);
            t = new BenchmarkProfile(dnaDistanceT, atDistance, dnaString.Length);
            c = new BenchmarkProfile(dnaDistanceC, acDistance, dnaString.Length);

            a.Calculate(t);
            t.Calculate(a);
            c.Calculate(a);

            Point = Trilaterate(dnaDistanceA, dnaDistanceT, dnaDistanceC);
        }

        // Standard 2D trilateration (the same linearization GPS receivers use): with A/T/C fixed
        // at known positions, subtracting the |P-A|^2=dA^2 equation from |P-T|^2=dT^2 and from
        // |P-C|^2=dC^2 eliminates the quadratic terms and leaves two linear equations for (x,y).
        private static (double X, double Y) Trilaterate(double dA, double dT, double dC)
        {
            double x = (dA * dA - dT * dT + atDistance * atDistance) / (2 * atDistance);
            double y = (dA * dA - dC * dC + acDistance * acDistance - 2 * cRefX * x) / (2 * cRefY);
            return (x, y);
        }

        private PathsD BuildPath()
        {
            var path = new PathsD();
            path.Add(Clipper.MakePath(new double[] { c.D * 100, c.H * 100, a.D * 100, a.H * 100, t.D * 100, t.H * 100 }));
            return path;
        }

        public static double CompareAreas(CATProfile x, CATProfile y)
        {
            var xPath = x.BuildPath();
            var yPath = y.BuildPath();
            var intersect = Clipper.Intersect(xPath, yPath, FillRule.NonZero, 8);
            return Clipper.Area(intersect) / Math.Min(Math.Abs(Clipper.Area(xPath)), Math.Abs(Clipper.Area(yPath)));
        }

        public static double CompareDistances(CATProfile x, CATProfile y)
        {
            var dC = x.c.D - y.c.D;
            var hC = x.c.H - y.c.H;
            var dA = x.a.D - y.a.D;
            var hA = x.a.H - y.a.H;
            var dT = x.t.D - y.t.D;
            var hT = x.t.H - y.t.H;

            var resC = Math.Sqrt(dC * dC + hC * hC);
            var resA = Math.Sqrt(dA * dA + hA * hA);
            var resT = Math.Sqrt(dT * dT + hT * hT);

            return 1 - (resC + resA + resT) / 3d;
        }

        public static double CompareTrilateration(CATProfile x, CATProfile y)
        {
            var dx = x.Point.X - y.Point.X;
            var dy = x.Point.Y - y.Point.Y;
            return 1 - Math.Sqrt(dx * dx + dy * dy);
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
            return (this.dnaDistance * this.dnaDistance + benchmarkDistance * benchmarkDistance - profile.dnaDistance * profile.dnaDistance) / (2 * this.dnaDistance * benchmarkDistance);
        }

        private double BenchmarkD(double benchmarkCOS)
        {
            return this.dnaDistance * benchmarkCOS;
        }

        private double BenchmarkH(double benchmarkCOS)
        {
            var d = this.dnaDistance * benchmarkCOS;
            return Math.Sqrt(this.dnaDistance * this.dnaDistance - d * d);
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