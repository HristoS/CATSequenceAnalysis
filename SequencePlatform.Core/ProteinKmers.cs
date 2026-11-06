using System.Buffers.Binary;
using System.IO.Compression;
using Npgsql;
using NpgsqlTypes;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <summary>
    /// A record's amino-acid k-mer counts (S7.R, migration 0014) over its translated coding sequences, as
    /// <see cref="AminoAcidKmers.Count"/> counts them. Payload v1: the 20^k counts as uint16 little-endian, deflate-compressed; a count
    /// above 65,535 is stored as 65,535. Like <see cref="RecordKmers"/>, a row also records why a record has no counts. Only 4-mers
    /// are stored since migration 0016; <see cref="Shorter"/> gives the 3- and 2-mers.
    /// </summary>
    public static class ProteinKmers
    {
        public static readonly IReadOnlyList<int> Ks = new[] { 4 };
        public const string Done = "done";
        public const string Mismatch = "mismatch";
        public const string Missing = "missing";
        public const string NoCds = "no_cds";
        private const short formatVersion = 1;

        public static byte[] Encode(int[] counts)
        {
            var raw = new byte[2 * counts.Length];
            for (int i = 0; i < counts.Length; i++)
                BinaryPrimitives.WriteUInt16LittleEndian(raw.AsSpan(2 * i), (ushort)Math.Clamp(counts[i], 0, ushort.MaxValue));
            using var output = new MemoryStream();
            using (var deflate = new DeflateStream(output, CompressionLevel.Optimal, leaveOpen: true))
                deflate.Write(raw);
            return output.ToArray();
        }

        /// <summary>
        /// The (k-1)-mer counts of k-mer counts, summed over the last letter: short only of each gene's last (k-1)-mer, which is how
        /// CVTree's 3-mers come from the stored 4-mers (results/CvTree_Derived3mers.csv), and its 2-mers from those 3-mers.
        /// </summary>
        public static int[] Shorter(int[] counts)
        {
            int letters = AminoAcidKmers.Alphabet.Length;
            var shorter = new int[counts.Length / letters];
            for (int w = 0; w < counts.Length; w++) shorter[w / letters] += counts[w];
            return shorter;
        }

        public static int[] Decode(byte[] payload, int k)
        {
            int space = (int)Math.Pow(AminoAcidKmers.Alphabet.Length, k);
            var raw = new byte[2 * space];
            using (var deflate = new DeflateStream(new MemoryStream(payload), CompressionMode.Decompress))
                deflate.ReadExactly(raw);
            var counts = new int[space];
            for (int i = 0; i < space; i++)
                counts[i] = BinaryPrimitives.ReadUInt16LittleEndian(raw.AsSpan(2 * i));
            return counts;
        }

        /// <param name="counts">Required for <see cref="Done"/>, null otherwise.</param>
        public static async Task WriteAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, long sequenceId, int k, string status,
            int? geneCount, int[]? counts)
        {
            if ((status == Done) != (counts is not null))
                throw new ArgumentException($"Counts go with status '{Done}' and only with it.", nameof(counts));
            if (!Ks.Contains(k)) throw new ArgumentOutOfRangeException(nameof(k));

            await using var command = new NpgsqlCommand(
                "insert into protein_kmer (sequence_id, k, format_version, gene_count, payload, status) values ($1, $2, $3, $4, $5, $6) " +
                "on conflict (sequence_id, k) do update set format_version = excluded.format_version, gene_count = excluded.gene_count, " +
                "payload = excluded.payload, status = excluded.status, checked_at = now()",
                connection, transaction);
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue((short)k);
            command.Parameters.AddWithValue(formatVersion);
            command.Parameters.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Integer, Value = (object?)geneCount ?? DBNull.Value });
            command.Parameters.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bytea, Value = counts is null ? DBNull.Value : Encode(counts) });
            command.Parameters.AddWithValue(status);
            await command.ExecuteNonQueryAsync();
        }

        /// <summary>The record's counts for <paramref name="k"/>, or null when it has none stored.</summary>
        public static async Task<int[]?> ReadAsync(NpgsqlDataSource dataSource, long sequenceId, int k)
        {
            await using var command = dataSource.CreateCommand("select payload from protein_kmer where sequence_id = $1 and k = $2 and status = $3");
            command.Parameters.AddWithValue(sequenceId);
            command.Parameters.AddWithValue((short)k);
            command.Parameters.AddWithValue(Done);
            return await command.ExecuteScalarAsync() is byte[] payload ? Decode(payload, k) : null;
        }
    }

    /// <summary>BV-BRC protein families named once (migration 0014): the number of a PGF_ id and its product.</summary>
    public static class ProteinFamilies
    {
        /// <summary>The number of a PGF_ id (PGF_00421071 -> 421071), or <see cref="ZCurveGeneProfile.NoFamily"/>.</summary>
        public static int Parse(string? pgfamId) =>
            pgfamId is not null && pgfamId.StartsWith("PGF_", StringComparison.Ordinal) && int.TryParse(pgfamId.AsSpan(4), out var id) ? id : ZCurveGeneProfile.NoFamily;

        /// <summary>Adds families not stored yet; an existing product is kept, and the first of duplicate ids wins.</summary>
        public static async Task AddAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, IReadOnlyCollection<(int Id, string Product)> families)
        {
            if (families.Count == 0) return;
            // In id order, so concurrent genomes sharing families lock their rows in the same order and can't deadlock.
            var ordered = families.DistinctBy(f => f.Id).OrderBy(f => f.Id).ToList();
            await using var command = new NpgsqlCommand(
                "insert into protein_family (id, product) select * from unnest($1::integer[], $2::text[]) on conflict (id) do nothing",
                connection, transaction);
            command.Parameters.AddWithValue(ordered.Select(f => f.Id).ToArray());
            command.Parameters.AddWithValue(ordered.Select(f => f.Product).ToArray());
            await command.ExecuteNonQueryAsync();
        }

        /// <summary>The PGF_ id as BV-BRC writes it (421071 -> PGF_00421071).</summary>
        public static string Format(int id) => $"PGF_{id:D8}";

        /// <summary>The stored product of each of these families; families not stored are left out.</summary>
        public static async Task<Dictionary<int, string>> ProductsAsync(NpgsqlDataSource dataSource, IReadOnlyCollection<int> ids,
            CancellationToken cancellationToken = default)
        {
            var products = new Dictionary<int, string>();
            if (ids.Count == 0) return products;
            await using var command = dataSource.CreateCommand("select id, product from protein_family where id = any($1) and product is not null");
            command.Parameters.AddWithValue(ids.ToArray());
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken)) products[reader.GetInt32(0)] = reader.GetString(1);
            return products;
        }
    }
}
