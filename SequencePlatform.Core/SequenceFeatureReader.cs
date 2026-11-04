using System.Security.Cryptography;

namespace SequencePlatform.Core
{
    /// <summary>
    /// The D8 opt-ins (2026-10-15), computed in the same streaming pass that feeds
    /// <c>ZCurveWindowedProfile.FromStream</c>, so the DNA is never held or read twice. Wrap the record's
    /// sequence reader in this, hand it to <c>FromStream</c>, then call <see cref="Finish"/>.
    ///
    /// Normalization, so line width and letter case can't change the values: only ASCII letters count,
    /// uppercased; whitespace, digits, gaps and other symbols are skipped as if absent.
    ///  - SHA-256 over those letters as ASCII bytes. Two records are "identical" when this matches.
    ///  - Dinucleotide counts over adjacent pairs of A, C, G, T, index 4 * first + second in A, C, G, T
    ///    order (AA, AC, AG, AT, CA, ..., TT). Any other letter (N, IUPAC codes) breaks the pair.
    ///  - 6-mer counts the same way (D8 widened 2026-11-04): one strand, index the base-4 number of the six bases with the
    ///    first most significant (AAAAAA = 0, ..., TTTTTT = 4095); any other letter breaks the run. Every shorter k sums out
    ///    of them except near run ends, so they aren't summed into the dinucleotides; those stay counted directly.
    /// </summary>
    public sealed class SequenceFeatureReader : TextReader
    {
        public const string BaseOrder = "ACGT";
        public const int HexamerLength = 6;
        public const int HexamerCount = 1 << (2 * HexamerLength);
        private const int hexamerMask = HexamerCount - 1;

        private readonly TextReader inner;
        private readonly IncrementalHash sha256 = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);
        private readonly long[] dinucleotideCounts = new long[16];
        private readonly int[] hexamerCounts = new int[HexamerCount];
        private byte[] letterBytes = Array.Empty<byte>();
        private int previousBase = -1;
        private int hexamer;
        private int run;
        private bool reachedEnd;
        private bool finished;

        public SequenceFeatureReader(TextReader inner) => this.inner = inner;

        /// <summary>Every character handed on to the consumer, before normalization (the record's length_bp).</summary>
        public long CharactersRead { get; private set; }

        public override int Read(char[] buffer, int index, int count)
        {
            if (finished) throw new InvalidOperationException("Finish was already called.");
            int read = inner.Read(buffer, index, count);
            if (read == 0) reachedEnd = true;
            else Accumulate(buffer.AsSpan(index, read));
            CharactersRead += read;
            return read;
        }

        private void Accumulate(ReadOnlySpan<char> chunk)
        {
            if (letterBytes.Length < chunk.Length) letterBytes = new byte[chunk.Length];
            int letters = 0;
            foreach (var raw in chunk)
            {
                if (!char.IsAsciiLetter(raw)) continue;
                var c = char.ToUpperInvariant(raw);
                letterBytes[letters++] = (byte)c;

                int baseIndex = c switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => -1 };
                if (baseIndex >= 0 && previousBase >= 0) dinucleotideCounts[4 * previousBase + baseIndex]++;
                previousBase = baseIndex;
                if (baseIndex < 0) run = 0;
                else
                {
                    hexamer = ((hexamer << 2) | baseIndex) & hexamerMask;
                    if (run < HexamerLength) run++;   // capped, so a run longer than int.MaxValue can't wrap
                    if (run == HexamerLength) hexamerCounts[hexamer] = checked(hexamerCounts[hexamer] + 1);
                }
            }
            sha256.AppendData(letterBytes, 0, letters);
        }

        /// <summary>
        /// The record's features. Only valid once the consumer has read to the end, which FromStream does;
        /// a partly read record would give a hash of a prefix, so that throws instead.
        /// </summary>
        public (byte[] Sha256, long[] DinucleotideCounts) Finish()
        {
            if (!reachedEnd) throw new InvalidOperationException("The sequence wasn't read to the end.");
            if (finished) throw new InvalidOperationException("Finish was already called.");
            finished = true;
            return (sha256.GetHashAndReset(), (long[])dinucleotideCounts.Clone());
        }

        /// <summary>The record's 6-mer counts, once <see cref="Finish"/> has been called.</summary>
        public int[] HexamerCounts => finished ? (int[])hexamerCounts.Clone()
            : throw new InvalidOperationException("Call Finish first: the counts are only complete at the end of the sequence.");

        protected override void Dispose(bool disposing)
        {
            if (disposing)
            {
                inner.Dispose();
                sha256.Dispose();
            }
            base.Dispose(disposing);
        }
    }
}
