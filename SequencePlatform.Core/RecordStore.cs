using System.Security.Cryptography;
using Npgsql;
using NpgsqlTypes;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="HeaderId">The FASTA header id as <see cref="FastaReader"/> parses it; a provider may store a cleaned-up form.</param>
    /// <param name="SequenceSha256">Null, like <paramref name="DinucleotideCounts"/> and <paramref name="HexamerCounts"/>, only for a profile upload: there was no DNA.</param>
    public sealed record ProfiledRecord(string HeaderId, string Description, long LengthBp, double? GcFraction,
        IReadOnlyList<(int Factor, ZCurveWindowedProfile Level)> Levels, byte[]? SequenceSha256, long[]? DinucleotideCounts,
        int[]? HexamerCounts = null);

    /// <param name="RawSha256">Of the file as downloaded or uploaded; the file itself is discarded.</param>
    /// <param name="OwnerId">The uploading user (D3); null for every downloaded source.</param>
    public sealed record SourceRow(string Provider, string? Url, string? Accession, byte[] RawSha256, long RawByteSize, string? License,
        long? OwnerId = null);

    /// <param name="AnnotationsJson">Serialized jsonb; free-form, per provider.</param>
    public sealed record RecordRow(string RecordId, string? Name, string Kind, int? NcbiTaxid, Lineage Lineage,
        string AnnotationsJson, string? ProviderMd5 = null, string DataClass = "public_microbe");

    /// <param name="SequenceIds">The new records' ids, in the order given.</param>
    /// <param name="SkippedLevels">Levels over <see cref="ProfileLevels.MaxInlineBytes"/>, not stored (they need object storage).</param>
    public sealed record StoredSource(long SourceId, IReadOnlyList<long> SequenceIds, IReadOnlyList<(string RecordId, int WindowSize)> SkippedLevels)
    {
        public bool FineLevelSkipped => SkippedLevels.Any(l => l.WindowSize == RecordStore.FineWindowSize);
    }

    /// <summary>
    /// The provider-independent half of ingestion (docs/database-decision.md section 4.1). It profiles every
    /// record of a FASTA file, streaming each one through <see cref="SequenceFeatureReader"/> into
    /// <see cref="ZCurveWindowedProfile.FromStream"/>, so no record is ever held as a whole string. It then stores
    /// the source, its records and their ladder levels in the caller's transaction. The BV-BRC Ingest, the
    /// bundled-corpus loader and the P4 upload use it. Deleting the file stays the caller's job.
    /// </summary>
    public static class RecordStore
    {
        public const string Method = "zcurve-windowed";
        public const string MethodVersion = "1";
        public const int FineWindowSize = ZCurveWindowedProfile.DefaultWindowSize;

        /// <summary>
        /// The platform's ladder: every level rounded to float16, as <see cref="ProfileLevels"/> stores it (ZCWP version 2), so the
        /// levels a query builds from its DNA equal the stored ones of the same sequence and an identical record scores exactly 1.
        /// </summary>
        public static ZCurveProfileLadder Ladder(ZCurveWindowedProfile fine) => new(fine, halfPrecision: true);

        public static IEnumerable<ProfiledRecord> Profile(string fastaPath)
        {
            foreach (var (headerId, description) in FastaReader.ReadHeaders(fastaPath))
                yield return Profile(headerId, description, FastaReader.OpenSequenceReader(fastaPath, headerId));
        }

        /// <summary>The same for FASTA text that can be read only once, such as an upload (P4 S4.2b): one pass, nothing materialized.</summary>
        public static IEnumerable<ProfiledRecord> Profile(TextReader fasta)
        {
            foreach (var (headerId, description, sequence) in FastaReader.ReadRecordStreams(fasta))
                yield return Profile(headerId, description, sequence);
        }

        private static ProfiledRecord Profile(string headerId, string description, TextReader sequence)
        {
            using var features = new SequenceFeatureReader(sequence);
            var fine = ZCurveWindowedProfile.FromStream(features, FineWindowSize);
            var (sequenceSha256, dinucleotideCounts) = features.Finish();
            var (levels, gcFraction) = Levels(fine);
            return new ProfiledRecord(headerId, description.Trim(), features.CharactersRead, gcFraction, levels,
                sequenceSha256, dinucleotideCounts, features.HexamerCounts);
        }

        /// <summary>
        /// A profile upload (D2): its fine level as given, the coarser levels derived from it as for DNA. The length is what the
        /// uploader declares; it must agree with the window count, since FromStream drops only a trailing partial window.
        /// </summary>
        public static ProfiledRecord FromFineProfile(string recordId, ZCurveWindowedProfile fine, long lengthBp)
        {
            if (fine.WindowSize != FineWindowSize)
                throw new ArgumentException($"A profile upload is the {FineWindowSize} bp level; this one has {fine.WindowSize} bp windows.", nameof(fine));
            if (lengthBp < 1 || lengthBp / FineWindowSize != fine.WindowCount)
                throw new ArgumentException($"A length of {lengthBp} bp gives {Math.Max(lengthBp, 0) / FineWindowSize} windows of {FineWindowSize} bp, " +
                    $"but the profile has {fine.WindowCount}.", nameof(lengthBp));
            var (levels, gcFraction) = Levels(fine);
            return new ProfiledRecord(recordId, string.Empty, lengthBp, gcFraction, levels, null, null);
        }

        private static (List<(int, ZCurveWindowedProfile)> Levels, double? GcFraction) Levels(ZCurveWindowedProfile fine)
        {
            var levels = new List<(int, ZCurveWindowedProfile)>();
            if (fine.WindowCount == 0) return (levels, null);
            var ladder = Ladder(fine);
            foreach (var factor in ladder.Factors)
            {
                var level = ladder.Level(factor);
                if (level.WindowCount > 0) levels.Add((factor, level));
            }
            // z = +1 for G/C, -1 for A/T, so the mean z over all windows is 2*GC - 1.
            return (levels, (fine.Windows.Average(w => w.Z) + 1) / 2);
        }

        public static (byte[] Sha256, long ByteSize) HashFile(string path)
        {
            using var stream = File.OpenRead(path);
            return (SHA256.HashData(stream), stream.Length);
        }

        /// <summary>
        /// Inserts the source, then each record with its levels, all in the caller's transaction so a crawl mark
        /// can commit with them. A level too large to store inline is left out and reported in the result.
        /// </summary>
        public static async Task<StoredSource> StoreAsync(NpgsqlConnection connection, NpgsqlTransaction transaction, SourceRow source,
            IReadOnlyList<(ProfiledRecord Profile, RecordRow Row)> records)
        {
            long sourceId;
            await using (var insertSource = new NpgsqlCommand(
                "insert into sequence_source (provider, url, accession, retrieved_at, raw_sha256, raw_byte_size, license, owner_id) " +
                "values ($1, $2, $3, now(), $4, $5, $6, $7) returning id", connection, transaction))
            {
                var p = insertSource.Parameters;
                p.AddWithValue(source.Provider);
                p.AddWithValue(OrNull(source.Url));
                p.AddWithValue(OrNull(source.Accession));
                p.AddWithValue(source.RawSha256);
                p.AddWithValue(source.RawByteSize);
                p.AddWithValue(OrNull(source.License));
                p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bigint, Value = (object?)source.OwnerId ?? DBNull.Value });
                sourceId = (long)(await insertSource.ExecuteScalarAsync())!;
            }

            var sequenceIds = new List<long>();
            var skipped = new List<(string, int)>();
            foreach (var (profile, row) in records)
            {
                long sequenceId;
                await using (var insertRecord = new NpgsqlCommand(
                    "insert into sequence_record (source_id, record_id, name, kind, data_class, length_bp, gc_fraction, ncbi_taxid, " +
                    "tax_superkingdom, tax_realm, tax_kingdom, tax_phylum, tax_class, tax_order, tax_family, tax_genus, tax_species, annotations, " +
                    "sequence_sha256, dinucleotide_counts, provider_md5) " +
                    "values ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21) returning id",
                    connection, transaction))
                {
                    var p = insertRecord.Parameters;
                    var lineage = row.Lineage;
                    p.AddWithValue(sourceId);
                    p.AddWithValue(row.RecordId);
                    p.AddWithValue(row.Name is null ? DBNull.Value : row.Name);
                    p.AddWithValue(row.Kind);
                    p.AddWithValue(row.DataClass);
                    p.AddWithValue(profile.LengthBp);
                    p.AddWithValue(profile.GcFraction.HasValue ? (object)(float)profile.GcFraction.Value : DBNull.Value);
                    p.AddWithValue(row.NcbiTaxid.HasValue ? (object)row.NcbiTaxid.Value : DBNull.Value);
                    foreach (var rank in new[] { lineage.Superkingdom, lineage.Realm, lineage.Kingdom, lineage.Phylum, lineage.Class,
                                 lineage.Order, lineage.Family, lineage.Genus, lineage.Species })
                        p.AddWithValue(OrNull(rank));
                    p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Jsonb, Value = row.AnnotationsJson });
                    p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Bytea, Value = (object?)profile.SequenceSha256 ?? DBNull.Value });
                    p.Add(new NpgsqlParameter { NpgsqlDbType = NpgsqlDbType.Array | NpgsqlDbType.Bigint, Value = (object?)profile.DinucleotideCounts ?? DBNull.Value });
                    p.AddWithValue(OrNull(row.ProviderMd5));
                    sequenceId = (long)(await insertRecord.ExecuteScalarAsync())!;
                }
                sequenceIds.Add(sequenceId);
                if (profile.HexamerCounts is not null)
                    await RecordKmers.WriteAsync(connection, transaction, sequenceId, RecordKmers.Done, profile.HexamerCounts);

                foreach (var (factor, level) in profile.Levels)
                    if (!await ProfileLevels.InsertAsync(connection, transaction, sequenceId, Method, MethodVersion, factor, level))
                        skipped.Add((row.RecordId, level.WindowSize));
            }

            return new StoredSource(sourceId, sequenceIds, skipped);
        }

        private static object OrNull(string? value) => string.IsNullOrEmpty(value) ? DBNull.Value : value;
    }
}
