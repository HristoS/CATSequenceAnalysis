using System.Security.Cryptography;
using Npgsql;
using NpgsqlTypes;
using SequenceAnalyses;

namespace SequencePlatform.Core
{
    /// <param name="HeaderId">The FASTA header id as <see cref="FastaReader"/> parses it; a provider may store a cleaned-up form.</param>
    public sealed record ProfiledRecord(string HeaderId, string Description, long LengthBp, double? GcFraction,
        IReadOnlyList<(int Factor, ZCurveWindowedProfile Level)> Levels, byte[] SequenceSha256, long[] DinucleotideCounts);

    /// <param name="RawSha256">Of the file as downloaded or uploaded; the file itself is discarded.</param>
    public sealed record SourceRow(string Provider, string? Url, string? Accession, byte[] RawSha256, long RawByteSize, string? License);

    /// <param name="AnnotationsJson">Serialized jsonb; free-form, per provider.</param>
    public sealed record RecordRow(string RecordId, string? Name, string Kind, int? NcbiTaxid, Lineage Lineage,
        string AnnotationsJson, string? ProviderMd5 = null, string DataClass = "public_microbe");

    /// <param name="SkippedLevels">Levels over <see cref="ProfileLevels.MaxInlineBytes"/>, not stored (they need object storage).</param>
    public sealed record StoredSource(long SourceId, IReadOnlyList<(string RecordId, int WindowSize)> SkippedLevels)
    {
        public bool FineLevelSkipped => SkippedLevels.Any(l => l.WindowSize == RecordStore.FineWindowSize);
    }

    /// <summary>
    /// The provider-independent half of ingestion (docs/database-decision.md section 4.1). It profiles every
    /// record of a FASTA file, streaming each one through <see cref="SequenceFeatureReader"/> into
    /// <see cref="ZCurveWindowedProfile.FromStream"/>, so no record is ever held as a whole string. It then stores
    /// the source, its records and their ladder levels in the caller's transaction. The BV-BRC Ingest and the
    /// bundled-corpus loader use it, and the P4 FASTA upload is meant to. Deleting the file stays the caller's job.
    /// </summary>
    public static class RecordStore
    {
        public const string Method = "zcurve-windowed";
        public const string MethodVersion = "1";
        public const int FineWindowSize = ZCurveWindowedProfile.DefaultWindowSize;

        public static IEnumerable<ProfiledRecord> Profile(string fastaPath)
        {
            foreach (var (headerId, description) in FastaReader.ReadHeaders(fastaPath))
            {
                using var features = new SequenceFeatureReader(FastaReader.OpenSequenceReader(fastaPath, headerId));
                var fine = ZCurveWindowedProfile.FromStream(features, FineWindowSize);
                var (sequenceSha256, dinucleotideCounts) = features.Finish();

                var levels = new List<(int, ZCurveWindowedProfile)>();
                double? gcFraction = null;
                if (fine.WindowCount > 0)
                {
                    var ladder = new ZCurveProfileLadder(fine);
                    foreach (var factor in ladder.Factors)
                    {
                        var level = ladder.Level(factor);
                        if (level.WindowCount > 0) levels.Add((factor, level));
                    }
                    // z = +1 for G/C, -1 for A/T, so the mean z over all windows is 2*GC - 1.
                    gcFraction = (fine.Windows.Average(w => w.Z) + 1) / 2;
                }

                yield return new ProfiledRecord(headerId, description.Trim(), features.CharactersRead, gcFraction, levels,
                    sequenceSha256, dinucleotideCounts);
            }
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
                "insert into sequence_source (provider, url, accession, retrieved_at, raw_sha256, raw_byte_size, license) " +
                "values ($1, $2, $3, now(), $4, $5, $6) returning id", connection, transaction))
            {
                var p = insertSource.Parameters;
                p.AddWithValue(source.Provider);
                p.AddWithValue(OrNull(source.Url));
                p.AddWithValue(OrNull(source.Accession));
                p.AddWithValue(source.RawSha256);
                p.AddWithValue(source.RawByteSize);
                p.AddWithValue(OrNull(source.License));
                sourceId = (long)(await insertSource.ExecuteScalarAsync())!;
            }

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
                    p.AddWithValue(profile.SequenceSha256);
                    p.AddWithValue(profile.DinucleotideCounts);
                    p.AddWithValue(OrNull(row.ProviderMd5));
                    sequenceId = (long)(await insertRecord.ExecuteScalarAsync())!;
                }

                foreach (var (factor, level) in profile.Levels)
                    if (!await ProfileLevels.InsertAsync(connection, transaction, sequenceId, Method, MethodVersion, factor, level))
                        skipped.Add((row.RecordId, level.WindowSize));
            }

            return new StoredSource(sourceId, skipped);
        }

        private static object OrNull(string? value) => string.IsNullOrEmpty(value) ? DBNull.Value : value;
    }
}
