using System.Buffers.Binary;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Channels;
using Npgsql;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A protein query as a search job keeps it (migration 0015): the expected Z-curve windows of
    /// <see cref="ZCurveGeneProfile.QueryFromProtein"/>, the length and the SHA-256 of the letters, never the letters themselves.
    /// </summary>
    /// <param name="GeneticCode">The NCBI translation table the synonymous codons come from: 11, or 4 (TGA is Trp).</param>
    public sealed record ProteinQuery(double[] Axes, int LengthAa, byte[]? Sha256, int GeneticCode)
    {
        public const int StandardCode = 11;
        public const int MycoplasmaCode = 4;

        public int WindowCount => Axes.Length / ZCurveGeneProfile.FeatureCount;

        /// <summary>
        /// Upper-cases the letters, drops white space and one final stop (*), and checks them: letters only, no internal stop, not DNA,
        /// and <see cref="ProteinSearch.MinLengthAa"/> to <see cref="ProteinSearch.MaxLengthAa"/> long. Codon weights are uniform.
        /// </summary>
        public static ProteinQuery FromProtein(string protein, int geneticCode = StandardCode)
        {
            var code = geneticCode switch
            {
                StandardCode => ZCurveGeneProfile.standardCode,
                MycoplasmaCode => ZCurveGeneProfile.mycoplasmaCode,
                _ => throw new ArgumentException($"The genetic code is {StandardCode} (standard) or {MycoplasmaCode} (Mycoplasma and relatives, TGA = Trp)."),
            };
            var letters = new StringBuilder(protein.Length);
            foreach (char c in protein)
                if (!char.IsWhiteSpace(c)) letters.Append(char.ToUpperInvariant(c));
            if (letters.Length > 0 && letters[^1] == '*') letters.Length--;
            var text = letters.ToString();

            for (int i = 0; i < text.Length; i++)
                if (text[i] is < 'A' or > 'Z')
                    throw new ArgumentException(text[i] == '*'
                        ? $"The protein has an internal stop (*) at residue {i + 1:N0}; search one reading frame's protein at a time."
                        : $"A protein is amino-acid letters only; residue {i + 1:N0} is '{text[i]}'.");
            if (text.Length > 0 && text.All(c => c is 'A' or 'C' or 'G' or 'T' or 'U' or 'N'))
                throw new ArgumentException("This looks like DNA, not a protein: use find-similar or localize for DNA.");
            if (text.Length < ProteinSearch.MinLengthAa)
                throw new ArgumentException($"A protein under {ProteinSearch.MinLengthAa} aa ({text.Length:N0} here) can't be searched: genes that short " +
                    "aren't kept.");
            if (text.Length > ProteinSearch.MaxLengthAa)
                throw new ArgumentException($"A protein over {ProteinSearch.MaxLengthAa:N0} aa ({text.Length:N0} here) is longer than any known one.");

            return new ProteinQuery(ZCurveGeneProfile.QueryFromProtein(text, null, code), text.Length,
                SHA256.HashData(Encoding.ASCII.GetBytes(text)), geneticCode);
        }

        /// <summary>The windows as search_run.query_profile stores them: little-endian float64.</summary>
        public byte[] AxesToBinary()
        {
            var bytes = new byte[8 * Axes.Length];
            for (int i = 0; i < Axes.Length; i++) BinaryPrimitives.WriteDoubleLittleEndian(bytes.AsSpan(8 * i), Axes[i]);
            return bytes;
        }

        public static double[] AxesFromBinary(ReadOnlySpan<byte> bytes)
        {
            if (bytes.Length % (8 * ZCurveGeneProfile.FeatureCount) != 0) throw new FormatException("Not a protein query's windows.");
            var axes = new double[bytes.Length / 8];
            for (int i = 0; i < axes.Length; i++) axes[i] = BinaryPrimitives.ReadDoubleLittleEndian(bytes[(8 * i)..]);
            return axes;
        }
    }

    /// <param name="GeneStart">0-based on the record's stored strand.</param>
    /// <param name="GeneEnd">Exclusive; the stop codon is inside, as BV-BRC lists it.</param>
    /// <param name="Family">The BV-BRC PGfam number, or <see cref="ZCurveGeneProfile.NoFamily"/>.</param>
    /// <param name="P">How often a gene of another family scores this high against a protein of this length (<see cref="ProteinScoreNull"/>);
    /// null before the null is built.</param>
    /// <param name="Log10E">log10 of the expected number of genes from other families scoring this high in this search: p × genes scored.</param>
    public sealed record ProteinHit(long SequenceId, int GeneStart, int GeneEnd, bool Reverse, int Family, double Score, PValue? P = null,
        double? Log10E = null);

    /// <param name="Candidates">Genes scored.</param>
    /// <param name="WeakSignal">Always false: no length band from 50 aa found the right family first under 90% of the time
    /// (results/ProteinSearch_CatalogueV2.csv, 100% in every band). Kept because search_run.weak_signal is shared with the other searches.</param>
    public sealed record ProteinSearchResult(IReadOnlyList<ProteinHit> Hits, long Candidates, bool WeakSignal);

    /// <summary>
    /// Which stored genes encode a protein (P7): the query's expected windows (<see cref="ProteinQuery"/>) against the genes of the
    /// records the filter matches (<see cref="GeneProfiles"/>), by <see cref="ZCurveGeneProfile.Score"/> over the common windows. A
    /// protein of 500 aa or more is compared with the genes that long (333-position windows from the start codon); a shorter one, split
    /// into 5 equal windows as the stored short genes are, only with genes within <see cref="LengthBand"/> of its length, read from the
    /// 2-3 length buckets that band touches (S7.S, results/ShortProtein_*.csv). The buckets are streamed from the database for every
    /// query: one reader, decoding and scoring spread over the cores.
    ///
    /// Before S7.S (fixed windows only), right family at rank 1 over the catalogue from the protein alone was 4.3% with one window
    /// (167-332 aa), 84.4% with two and 98.8% with three or more (results/ProteinSearch_Catalogue*.csv).
    /// </summary>
    public static class ProteinSearch
    {
        public const int MinLengthAa = ZCurveGeneProfile.MinResidues;
        public const int MaxLengthAa = 40_000;
        public const double LengthBand = 0.1;
        private const int sharpWindows = 3;

        // Ascending from the worst hit: lower score, then (on a tie) the later record and gene, so the kept top is the same whichever
        // worker scored what.
        private static readonly Comparer<ProteinHit> worseFirst = Comparer<ProteinHit>.Create((a, b) =>
        {
            int byScore = a.Score.CompareTo(b.Score);
            if (byScore != 0) return byScore;
            int byRecord = b.SequenceId.CompareTo(a.SequenceId);
            return byRecord != 0 ? byRecord : b.GeneStart.CompareTo(a.GeneStart);
        });

        public static async Task<ProteinSearchResult> SearchAsync(NpgsqlDataSource dataSource, ProteinQuery query, CatalogueFilter? filter = null,
            int top = 50, CancellationToken cancellationToken = default)
        {
            ArgumentOutOfRangeException.ThrowIfLessThan(top, 1);
            var visible = await Catalogue.IdsAsync(dataSource, filter ?? new CatalogueFilter(), cancellationToken);
            var candidate = CandidateRule(query);
            var rows = Channel.CreateBounded<(long Id, byte[] Payload, byte[] Sha256)>(256);
            // A worker that fails stops the reader, which would otherwise wait on a full channel; its own error is what's thrown.
            using var stop = CancellationTokenSource.CreateLinkedTokenSource(cancellationToken);
            var workers = Enumerable.Range(0, Environment.ProcessorCount).Select(_ => Task.Run(async () =>
            {
                try { return await ScoreAsync(rows.Reader, query.Axes, candidate, top, stop.Token); }
                catch { stop.Cancel(); throw; }
            }, cancellationToken)).ToArray();

            try
            {
                await using var command = dataSource.CreateCommand(
                    "select sequence_id, payload, sha256 from gene_profile where bucket = any($1) and window_count > 0");
                command.Parameters.AddWithValue(GeneProfiles.BucketsFor(query.LengthAa, LengthBand).Select(b => (short)b).ToArray());
                await using var reader = await command.ExecuteReaderAsync(stop.Token);
                while (await reader.ReadAsync(stop.Token))
                {
                    long id = reader.GetInt64(0);
                    if (visible.Contains(id))
                        await rows.Writer.WriteAsync((id, reader.GetFieldValue<byte[]>(1), reader.GetFieldValue<byte[]>(2)), stop.Token);
                }
            }
            catch (OperationCanceledException) when (!cancellationToken.IsCancellationRequested)
            {
                // A worker failed; awaiting the workers below throws its error.
            }
            finally
            {
                rows.Writer.Complete();
            }

            var results = await Task.WhenAll(workers);
            long scored = results.Sum(r => r.Scored);
            var nullModel = await ProteinScoreNull.LoadAsync(dataSource, cancellationToken);
            var hits = results.SelectMany(r => r.Best.UnorderedItems.Select(item => item.Element))
                .OrderByDescending(h => h, worseFirst).Take(top)
                .Select(h => nullModel?.Assess(query.LengthAa, h.Score, scored) is { } e ? h with { P = e.P, Log10E = e.Log10E } : h).ToList();
            return new ProteinSearchResult(hits, scored, WeakSignal: false);
        }

        /// <summary>Which genes a query is compared with: from 500 aa the long genes, below that the genes within <see cref="LengthBand"/> of its length.</summary>
        public static Func<ZCurveGeneProfile.Gene, bool> CandidateRule(ProteinQuery query) => query.LengthAa >= ZCurveGeneProfile.LongResidues
            ? gene => gene.WindowCount >= Math.Min(query.WindowCount, sharpWindows)
            : gene => Math.Abs(gene.Residues - query.LengthAa) <= LengthBand * query.LengthAa;

        private static async Task<(PriorityQueue<ProteinHit, ProteinHit> Best, long Scored)> ScoreAsync(
            ChannelReader<(long Id, byte[] Payload, byte[] Sha256)> rows, double[] axes, Func<ZCurveGeneProfile.Gene, bool> candidate, int top,
            CancellationToken cancellationToken)
        {
            var best = new PriorityQueue<ProteinHit, ProteinHit>(top + 1, worseFirst);
            long scored = 0;
            await foreach (var (id, payload, sha256) in rows.ReadAllAsync(cancellationToken))
            {
                if (!SHA256.HashData(payload).AsSpan().SequenceEqual(sha256))
                    throw new InvalidDataException($"Gene profile of record {id} doesn't match its SHA-256.");
                var profile = ZCurveGeneProfile.FromBinary(payload);
                for (int g = 0; g < profile.Genes.Count; g++)
                {
                    var gene = profile.Genes[g];
                    if (!candidate(gene)) continue;
                    scored++;
                    var hit = new ProteinHit(id, gene.Start, gene.End, gene.Reverse, gene.Family, profile.Score(axes, g));
                    if (best.Count < top) best.Enqueue(hit, hit);
                    else if (worseFirst.Compare(hit, best.Peek()) > 0) best.EnqueueDequeue(hit, hit);
                }
            }
            return (best, scored);
        }
    }
}
