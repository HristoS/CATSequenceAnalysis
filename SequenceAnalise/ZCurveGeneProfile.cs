using System;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Linq;

namespace SequenceAnalyses
{
    /// <summary>
    /// Per-gene Z-curve profile for protein search (platform steps S7.R, decided 2026-11-04, and S7.S, 2026-11-05): each coding sequence,
    /// in coding orientation, keeps only codon positions 1 and 2 (the bases a protein fixes). Per window and axis it keeps three Legendre
    /// moments of the base's position t in [-1, 1]: the mean (P0), the trend (P1) and the curvature (P2), each scaled by sqrt(2k + 1) so all
    /// three carry the same noise on random sequence. Windows anchored at the start codon line up exactly with a query built from the
    /// protein, which is what makes the moments informative; on unanchored genome windows they hurt (results/GenomeMoments_*.csv), so the
    /// genome ladder keeps means.
    ///
    /// Window rule, by the protein's length L (codons without the stop): from <see cref="LongResidues"/> aa, windows of
    /// <see cref="WindowBases"/> such bases (500 bp of gene) from the start codon, 97.0% top-1 among 225,827 genes from the protein alone
    /// (results/T3_Pool300.csv); from <see cref="MinResidues"/> aa, <see cref="SplitWindows"/> equal windows over the whole gene, which a
    /// query is compared with only among genes of similar length (results/ShortProtein_*.csv); shorter genes aren't kept. Features are
    /// stored as 16-bit halves, which changed no measured number (results/CdsPrecision_Search.csv); 8 bits cost up to 6 points.
    /// </summary>
    public sealed class ZCurveGeneProfile
    {
        public const int WindowBases = 333;
        public const int SplitWindows = 5;
        public const int MinResidues = 50;
        public const int LongResidues = 500;
        public const int FeatureCount = 9;
        public const int NoFamily = -1;
        // NCBI translation tables in TCAG order (index = 16 * first + 4 * second + third): table 11 (bacteria, archaea; internal codons as
        // table 1) and table 4 (Mycoplasma and relatives, TGA = Trp).
        public const string standardCode = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        public const string mycoplasmaCode = "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        // Table 12, the CUG-Ser yeasts' nuclear code: CTG reads as Ser. Table 1, every other nuclear eukaryote here, is standardCode.
        public const string alternativeYeastCode = "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
        private const string codonBases = "TCAG";
        private const uint binaryMagic = 0x5047435A; // "ZCGP" read as a little-endian uint32
        private const ushort singleSequenceVersion = 2;
        private const ushort splicedVersion = 3;
        private const int binaryHeaderSize = 16;
        // Every float16 bit pattern as the double it equals (Half to double is exact), so reading a stored feature is one load (ScoreRow).
        private static readonly double[] halfValues = Enumerable.Range(0, 1 << 16).Select(bits => (double)BitConverter.UInt16BitsToHalf((ushort)bits)).ToArray();

        private readonly Gene[] genes;
        // Structure of arrays over every gene's windows in gene order: features[f][window], f = 3 * order + axis.
        private readonly Half[][] features;

        public IReadOnlyList<Gene> Genes => genes;

        public int WindowCount { get; }

        /// <summary>The ZCGP version <see cref="ToBinary"/> writes: 2 while every gene is one unspliced stretch of the record itself, else 3.</summary>
        public int BinaryFormatVersion => genes.Any(g => g.Chromosome != 0 || g.Introns != 0) ? splicedVersion : singleSequenceVersion;

        private ZCurveGeneProfile(Gene[] genes, Half[][] features)
        {
            this.genes = genes;
            this.features = features;
            WindowCount = features[0].Length;
        }

        /// <summary>
        /// One coding sequence on the record: [Start, End) in record coordinates, its strand, BV-BRC PGfam number and protein length
        /// (codons without the stop), and where its windows are. A record that is a whole assembly (an NCBI eukaryote) numbers its
        /// sequences: <paramref name="Chromosome"/> is the gene's, and Start and End are on it. <paramref name="Introns"/> is how much of
        /// [Start, End) a spliced gene leaves out, so its coding length is End - Start - Introns.
        /// </summary>
        public readonly record struct Gene(int Start, int End, bool Reverse, int Family, int FirstWindow, int WindowCount, int Residues,
            int Chromosome = 0, int Introns = 0);

        /// <summary>
        /// Profiles every coding sequence of one record of at least <see cref="MinResidues"/> aa, in order of start. Coordinates are 0-based,
        /// end exclusive; a reverse-strand gene is read as the reverse complement. A final stop codon is left out, so annotations with or
        /// without it give the same windows. Bases other than A, C, G, T contribute nothing, as in <see cref="ZCurveProfile"/>.
        /// </summary>
        public static ZCurveGeneProfile FromRecord(string record, IEnumerable<(int Start, int End, bool Reverse, int Family)> codingSequences,
            string geneticCode = standardCode) =>
            FromCodingSequences(codingSequences.Select(c => (0, c.Start, c.End, c.Reverse, c.Family, CodingSequence(record, c.Start, c.End, c.Reverse))),
                geneticCode);

        /// <summary>
        /// Profiles coding sequences given already in coding orientation and spliced, each with where it lies: its sequence's number,
        /// [Start, End) on it spanning every exon, and its strand. Genes under <see cref="MinResidues"/> aa are left out, as in
        /// <see cref="FromRecord"/>; the rest are kept in order of chromosome, then start.
        /// </summary>
        public static ZCurveGeneProfile FromCodingSequences(
            IEnumerable<(int Chromosome, int Start, int End, bool Reverse, int Family, string Gene)> codingSequences, string geneticCode = standardCode)
        {
            var genes = new List<(Gene Gene, double[] Windows)>();
            foreach (var (chromosome, start, end, reverse, family, gene) in codingSequences.OrderBy(c => c.Chromosome).ThenBy(c => c.Start))
            {
                int introns = end - start - gene.Length;
                if (chromosome < 0 || introns < 0) throw new ArgumentException($"A coding sequence of {gene.Length} bp can't lie on [{start}, {end}).");
                int residues = Residues(gene, geneticCode);
                if (residues < MinResidues) continue;
                var windows = QueryFromGene(gene, geneticCode);
                genes.Add((new Gene(start, end, reverse, family, 0, windows.Length / FeatureCount, residues, chromosome, introns), windows));
            }
            return Build(genes.Select(g => (g.Gene, (Func<int, int, double>)((w, f) => g.Windows[w * FeatureCount + f]))).ToList());
        }

        /// <summary>The genes at <paramref name="indices"/> (a storage bucket), with their windows, kept in order of start.</summary>
        public ZCurveGeneProfile Select(IEnumerable<int> indices) =>
            Build(indices.Order().Select(g => (genes[g], (Func<int, int, double>)((w, f) => (double)features[f][genes[g].FirstWindow + w]))).ToList());

        /// <summary>The genes of several profiles of one record (its storage buckets) as one profile, in order of chromosome, then start.</summary>
        public static ZCurveGeneProfile Combine(IEnumerable<ZCurveGeneProfile> parts) =>
            Build(parts.SelectMany(p => p.genes.Select(g => (g, (Func<int, int, double>)((w, f) => (double)p.features[f][g.FirstWindow + w]))))
                .OrderBy(g => g.Item1.Chromosome).ThenBy(g => g.Item1.Start).ToList());

        private static ZCurveGeneProfile Build(List<(Gene Gene, Func<int, int, double> Feature)> parts)
        {
            var genes = new Gene[parts.Count];
            int total = parts.Sum(p => p.Gene.WindowCount);
            var features = Enumerable.Range(0, FeatureCount).Select(_ => new Half[total]).ToArray();
            int first = 0;
            for (int g = 0; g < parts.Count; g++)
            {
                var (gene, feature) = parts[g];
                genes[g] = gene with { FirstWindow = first };
                for (int w = 0; w < gene.WindowCount; w++)
                    for (int f = 0; f < FeatureCount; f++)
                        features[f][first + w] = (Half)feature(w, f);
                first += gene.WindowCount;
            }
            return new ZCurveGeneProfile(genes, features);
        }

        /// <summary>A coding sequence of a record in coding orientation: [start, end), the reverse complement when <paramref name="reverse"/>.</summary>
        public static string CodingSequence(string record, int start, int end, bool reverse)
        {
            if (start < 0 || end > record.Length || end <= start) throw new ArgumentOutOfRangeException(nameof(start));
            return reverse ? ReverseComplement(record.AsSpan(start, end - start)) : record.Substring(start, end - start);
        }

        /// <summary>
        /// The protein of a coding sequence: every complete codon, a final stop left out; an internal stop stays as '*' and a codon
        /// with a base other than A, C, G, T becomes 'X', so callers can see both.
        /// </summary>
        public static string Translate(string gene, string geneticCode = standardCode)
        {
            int codons = gene.Length / 3;
            var protein = new char[codons];
            for (int c = 0; c < codons; c++)
                protein[c] = CodonIndex(gene, 3 * c) is int index && index >= 0 ? geneticCode[index] : 'X';
            return codons > 0 && protein[^1] == '*' ? new string(protein, 0, codons - 1) : new string(protein);
        }

        /// <summary>A gene's windows as a query, window by window (<see cref="FeatureCount"/> values each), from its DNA in coding orientation.</summary>
        public static double[] QueryFromGene(string gene, string geneticCode = standardCode)
        {
            int residues = Residues(gene, geneticCode);
            return Windows(2 * residues, WindowOf(residues), i => Axes(gene[3 * (i / 2) + i % 2]));
        }

        /// <summary>
        /// A protein's windows as a query: each codon position 1 or 2 gets the mean axes of the amino acid's synonymous codons, weighted
        /// by <paramref name="codonUsage"/> (64 counts in TCAG order; uniform when null or when it has none of the synonyms). For 17
        /// amino acids those bases are fixed, so this is exact; for Leu, Arg and Ser it keeps the axes their codons share.
        /// </summary>
        public static double[] QueryFromProtein(string protein, long[]? codonUsage = null, string geneticCode = standardCode)
        {
            var table = ExpectedAxes(codonUsage, geneticCode);
            return Windows(2 * protein.Length, WindowOf(protein.Length), i =>
            {
                int letter = char.ToUpperInvariant(protein[i / 2]) - 'A';
                return letter is >= 0 and < 26 ? table[2 * letter + i % 2] : (0, 0, 0);
            });
        }

        /// <summary>Per amino acid letter and codon position 1 or 2, the usage-weighted mean axes of its synonyms; indexed 2 * (letter - 'A') + position.</summary>
        public static (double X, double Y, double Z)[] ExpectedAxes(long[]? codonUsage = null, string geneticCode = standardCode)
        {
            var table = new (double X, double Y, double Z)[2 * 26];
            foreach (var codons in Enumerable.Range(0, 64).Where(c => geneticCode[c] != '*').GroupBy(c => geneticCode[c]))
            {
                long total = codonUsage is null ? 0 : codons.Sum(c => codonUsage[c]);
                int synonyms = codons.Count();
                foreach (var codon in codons)
                {
                    double weight = total > 0 ? codonUsage![codon] / (double)total : 1.0 / synonyms;
                    for (int position = 0; position < 2; position++)
                    {
                        var (x, y, z) = Axes(codonBases[position == 0 ? codon / 16 : codon / 4 % 4]);
                        ref var cell = ref table[2 * (codons.Key - 'A') + position];
                        cell = (cell.X + weight * x, cell.Y + weight * y, cell.Z + weight * z);
                    }
                }
            }
            return table;
        }

        /// <summary>Gene <paramref name="gene"/>'s windows as stored, window by window: <see cref="FeatureCount"/> values each, 3 * order + axis.</summary>
        public double[] WindowsOf(int gene)
        {
            var (first, count) = (genes[gene].FirstWindow, genes[gene].WindowCount);
            var windows = new double[count * FeatureCount];
            for (int w = 0; w < count; w++)
                for (int f = 0; f < FeatureCount; f++)
                    windows[w * FeatureCount + f] = (double)features[f][first + w];
            return windows;
        }

        /// <summary>Each feature's mean over every window of every gene (3 * order + axis); zeros when there's no window.</summary>
        public double[] MeanWindow()
        {
            var mean = new double[FeatureCount];
            for (int f = 0; f < FeatureCount && WindowCount > 0; f++)
            {
                double sum = 0;
                foreach (var value in features[f]) sum += (double)value;
                mean[f] = sum / WindowCount;
            }
            return mean;
        }

        /// <summary>
        /// 1 - the mean over the common windows of the per-window Euclidean distance between a query (from <see cref="QueryFromGene"/>
        /// or <see cref="QueryFromProtein"/>) and gene <paramref name="gene"/>, the convention of <see cref="ZCurveWindowedProfile.BestAlignment"/>.
        /// NaN when either has no window. A query split into <see cref="SplitWindows"/> means something only against genes of similar
        /// length, which the caller chooses.
        /// </summary>
        public double Score(double[] query, int gene)
        {
            int count = Math.Min(query.Length / FeatureCount, genes[gene].WindowCount);
            if (count == 0) return double.NaN;
            int first = genes[gene].FirstWindow;
            double sum = 0;
            for (int w = 0; w < count; w++)
            {
                double squared = 0;
                for (int f = 0; f < FeatureCount; f++)
                {
                    double difference = query[w * FeatureCount + f] - (double)features[f][first + w];
                    squared += difference * difference;
                }
                sum += Math.Sqrt(squared);
            }
            return 1 - sum / count;
        }

        /// <summary>
        /// Version 2: "ZCGP" | uint16 version | uint16 reserved | int32 geneCount | int32 windowCount | per gene, in order of start, three
        /// unsigned LEB128 varints: start minus the previous gene's start, End - Start, and (Family + 1) * 4 + 2 * (a stop codon ends it) +
        /// (reverse strand) | then the <see cref="FeatureCount"/> feature blocks, each windowCount float16 values, little-endian. A gene's
        /// residues and window count follow from its length and stop bit. Never any bases.
        ///
        /// Version 3 (<see cref="BinaryFormatVersion"/>, genes on several sequences or spliced) is version 2 with two more varints per gene,
        /// genes in order of chromosome, then start: first the chromosome minus the previous gene's, and the start counts from the previous
        /// gene's start only on the same chromosome, from 0 on a new one; after the flags, the introns. End - Start is the span, and the
        /// coding length End - Start - Introns is what residues and windows follow from.
        /// </summary>
        public byte[] ToBinary()
        {
            int n = WindowCount;
            int version = BinaryFormatVersion;
            var header = new List<byte>(8 * genes.Length);
            int previous = 0, previousChromosome = 0;
            foreach (var gene in genes)
            {
                if (gene.Chromosome < previousChromosome || (gene.Chromosome == previousChromosome && gene.Start < previous))
                    throw new InvalidOperationException("Genes must be in order of chromosome, then start.");
                if (gene.Chromosome != previousChromosome) previous = 0;
                int length = gene.End - gene.Start;
                bool stop = (length - gene.Introns) / 3 != gene.Residues;
                if (version == splicedVersion) WriteVarint(header, (ulong)(gene.Chromosome - previousChromosome));
                WriteVarint(header, (ulong)(gene.Start - previous));
                WriteVarint(header, (ulong)length);
                WriteVarint(header, ((ulong)(gene.Family + 1) << 2) | (stop ? 2UL : 0) | (gene.Reverse ? 1UL : 0));
                if (version == splicedVersion) WriteVarint(header, (ulong)gene.Introns);
                (previous, previousChromosome) = (gene.Start, gene.Chromosome);
            }
            var bytes = new byte[binaryHeaderSize + header.Count + 2 * FeatureCount * n];
            var span = bytes.AsSpan();
            BinaryPrimitives.WriteUInt32LittleEndian(span, binaryMagic);
            BinaryPrimitives.WriteUInt16LittleEndian(span[4..], (ushort)version);
            BinaryPrimitives.WriteInt32LittleEndian(span[8..], genes.Length);
            BinaryPrimitives.WriteInt32LittleEndian(span[12..], n);
            header.CopyTo(bytes, binaryHeaderSize);
            int offset = binaryHeaderSize + header.Count;
            for (int f = 0; f < FeatureCount; f++)
                for (int w = 0; w < n; w++)
                    BinaryPrimitives.WriteHalfLittleEndian(span[(offset + 2 * (f * n + w))..], features[f][w]);
            return bytes;
        }

        /// <summary>Inverse of <see cref="ToBinary"/>.</summary>
        public static ZCurveGeneProfile FromBinary(ReadOnlySpan<byte> bytes)
        {
            var genes = new Gene[GeneCountOf(bytes)];
            var (offset, n) = ReadGenes(bytes, genes);
            var features = new Half[FeatureCount][];
            for (int f = 0; f < FeatureCount; f++)
            {
                features[f] = new Half[n];
                for (int w = 0; w < n; w++)
                    features[f][w] = BinaryPrimitives.ReadHalfLittleEndian(bytes[(offset + 2 * (f * n + w))..]);
            }
            return new ZCurveGeneProfile(genes, features);
        }

        /// <summary>
        /// <see cref="Score"/> for every gene of a stored row (<see cref="ToBinary"/> bytes) that <paramref name="candidate"/> accepts, read in
        /// place: no profile is built and only the candidates' windows are read, so a search scanning many rows decodes what it compares and
        /// nothing else. Scores equal <see cref="FromBinary"/> then <see cref="Score"/> bit for bit (same values, same order of operations).
        /// </summary>
        /// <param name="genes">Scratch for the row's genes, at least <see cref="GeneCountOf"/> long; reused across rows by the caller.</param>
        /// <param name="hit">Called with each candidate gene and its score.</param>
        public static void ScoreRow(ReadOnlySpan<byte> bytes, double[] query, Gene[] genes, Func<Gene, bool> candidate, Action<Gene, double> hit)
        {
            int geneCount = GeneCountOf(bytes);
            var (offset, n) = ReadGenes(bytes, genes.AsSpan(0, geneCount));
            int queryWindows = query.Length / FeatureCount;
            for (int g = 0; g < geneCount; g++)
            {
                var gene = genes[g];
                if (!candidate(gene)) continue;
                int count = Math.Min(queryWindows, gene.WindowCount);
                if (count == 0)
                {
                    hit(gene, double.NaN);
                    continue;
                }
                double sum = 0;
                for (int w = 0; w < count; w++)
                {
                    double squared = 0;
                    int queryStart = w * FeatureCount;
                    int geneWindow = gene.FirstWindow + w;
                    for (int f = 0; f < FeatureCount; f++)
                    {
                        double difference = query[queryStart + f] -
                            halfValues[BinaryPrimitives.ReadUInt16LittleEndian(bytes[(offset + 2 * (f * n + geneWindow))..])];
                        squared += difference * difference;
                    }
                    sum += Math.Sqrt(squared);
                }
                hit(gene, 1 - sum / count);
            }
        }

        /// <summary>How many genes a stored row holds: the size <see cref="ScoreRow"/>'s scratch needs.</summary>
        public static int GeneCountOf(ReadOnlySpan<byte> bytes)
        {
            if (bytes.Length < binaryHeaderSize || BinaryPrimitives.ReadUInt32LittleEndian(bytes) != binaryMagic)
                throw new FormatException("Not a ZCGP profile.");
            int geneCount = BinaryPrimitives.ReadInt32LittleEndian(bytes[8..]);
            if (geneCount < 0) throw new FormatException("ZCGP header does not match the payload length.");
            return geneCount;
        }

        // The gene list of a stored row into genes (exactly its gene count long), checked against the window count and payload length.
        // Returns where the feature blocks start and the window count.
        private static (int FeatureOffset, int WindowCount) ReadGenes(ReadOnlySpan<byte> bytes, Span<Gene> genes)
        {
            var version = BinaryPrimitives.ReadUInt16LittleEndian(bytes[4..]);
            if (version != singleSequenceVersion && version != splicedVersion)
                throw new FormatException($"Unsupported ZCGP format version {version}.");
            bool spliced = version == splicedVersion;
            int n = BinaryPrimitives.ReadInt32LittleEndian(bytes[12..]);
            if (n < 0) throw new FormatException("ZCGP header does not match the payload length.");

            int offset = binaryHeaderSize, firstWindow = 0, start = 0, chromosome = 0;
            for (int g = 0; g < genes.Length; g++)
            {
                int nextChromosome = spliced ? checked(chromosome + (int)ReadVarint(bytes, ref offset)) : 0;
                if (nextChromosome != chromosome) (start, chromosome) = (0, nextChromosome);
                start += checked((int)ReadVarint(bytes, ref offset));
                int length = checked((int)ReadVarint(bytes, ref offset));
                ulong flags = ReadVarint(bytes, ref offset);
                int introns = spliced ? checked((int)ReadVarint(bytes, ref offset)) : 0;
                if (introns > length) throw new FormatException("ZCGP gene has more introns than span.");
                int residues = (length - introns) / 3 - ((flags & 2) != 0 ? 1 : 0);
                int window = WindowOf(residues);
                int windows = window == 0 ? 0 : 2 * residues / window;
                genes[g] = new Gene(start, start + length, (flags & 1) != 0, checked((int)(flags >> 2)) - 1, firstWindow, windows, residues,
                    chromosome, introns);
                firstWindow += windows;
            }
            if (firstWindow != n) throw new FormatException("ZCGP gene windows do not add up to the window count.");
            if (bytes.Length != offset + 2L * FeatureCount * n) throw new FormatException("ZCGP header does not match the payload length.");
            return (offset, n);
        }

        private static void WriteVarint(List<byte> destination, ulong value)
        {
            while (value >= 0x80)
            {
                destination.Add((byte)(value | 0x80));
                value >>= 7;
            }
            destination.Add((byte)value);
        }

        private static ulong ReadVarint(ReadOnlySpan<byte> source, ref int offset)
        {
            ulong value = 0;
            for (int shift = 0; ; shift += 7)
            {
                if (offset >= source.Length || shift > 63) throw new FormatException("ZCGP gene list is truncated.");
                byte b = source[offset++];
                value |= (ulong)(b & 0x7F) << shift;
                if (b < 0x80) return value;
            }
        }

        // Codons without a final stop: the protein's length.
        private static int Residues(string gene, string geneticCode)
        {
            int codons = gene.Length / 3;
            return codons > 0 && CodonIndex(gene, 3 * (codons - 1)) is int last && last >= 0 && geneticCode[last] == '*' ? codons - 1 : codons;
        }

        /// <summary>The window rule: codon positions 1 and 2 per window for a protein of this length; 0 when it's too short to keep.</summary>
        public static int WindowOf(int residues) =>
            residues >= LongResidues ? WindowBases : residues >= MinResidues ? 2 * residues / SplitWindows : 0;

        // Complete windows of `window` positions-1-and-2 bases over the first `bases`: per window, for k = 0, 1, 2 and each axis,
        // sqrt(2k + 1) / W * sum(a_i P_k(t_i)) with t_i = -1 + 2 (i + 0.5) / W; laid out window by window, 3 * k + axis.
        private static double[] Windows(int bases, int window, Func<int, (double X, double Y, double Z)> axesAt)
        {
            int count = window == 0 ? 0 : bases / window;
            var result = new double[count * FeatureCount];
            for (int w = 0; w < count; w++)
            {
                for (int i = 0; i < window; i++)
                {
                    var (x, y, z) = axesAt(w * window + i);
                    double t = -1 + 2 * (i + 0.5) / window;
                    double p2 = (3 * t * t - 1) / 2;
                    int slot = w * FeatureCount;
                    result[slot] += x;
                    result[slot + 1] += y;
                    result[slot + 2] += z;
                    result[slot + 3] += x * t;
                    result[slot + 4] += y * t;
                    result[slot + 5] += z * t;
                    result[slot + 6] += x * p2;
                    result[slot + 7] += y * p2;
                    result[slot + 8] += z * p2;
                }
                for (int f = 0; f < FeatureCount; f++)
                    result[w * FeatureCount + f] *= Math.Sqrt(2 * (f / 3) + 1) / window;
            }
            return result;
        }

        private static (double X, double Y, double Z) Axes(char @base)
        {
            var (x, y, z) = ZCurveProfile.Axes(@base);
            return (x, y, z);
        }

        private static int CodonIndex(string sequence, int i)
        {
            int first = codonBases.IndexOf(char.ToUpperInvariant(sequence[i]));
            int second = codonBases.IndexOf(char.ToUpperInvariant(sequence[i + 1]));
            int third = codonBases.IndexOf(char.ToUpperInvariant(sequence[i + 2]));
            return first < 0 || second < 0 || third < 0 ? -1 : 16 * first + 4 * second + third;
        }

        private static string ReverseComplement(ReadOnlySpan<char> sequence)
        {
            var result = new char[sequence.Length];
            for (int i = 0; i < sequence.Length; i++)
                result[sequence.Length - 1 - i] = char.ToUpperInvariant(sequence[i]) switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => 'N' };
            return new string(result);
        }
    }
}
