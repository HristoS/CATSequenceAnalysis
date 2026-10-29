using System.Security.Claims;
using Microsoft.AspNetCore.Http.HttpResults;
using Npgsql;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>
    /// POST /api/sequences (P4 S4.2b, D2): a private upload, as FASTA or as a ZCWP profile. The body is the file itself, with the
    /// metadata in the query string: a multipart form would go through IFormFile, which buffers anything over 64 KB to a temp file.
    ///
    /// FASTA streams in one pass through <see cref="FastaReader.ReadRecordStreams"/> and <see cref="RecordStore.Profile(TextReader)"/>:
    /// every record is profiled as it arrives and never exists as a string, and the body is hashed on the way
    /// (<see cref="HashingReadStream"/>). Only the profile, the hashes and the metadata are stored, as for a download.
    /// </summary>
    internal static class SequenceUpload
    {
        private const int maxNameLength = 500;

        private static readonly HashSet<string> kinds = ["virus", "bacterium", "archaeon", "eukaryote", "plasmid", "organelle", "other"];
        private static readonly HashSet<string> dataClasses = ["restricted", "public_microbe", "public_eukaryote", "human"];

        /// <param name="kind">virus, bacterium, archaeon, eukaryote, plasmid, organelle or other: what every record of the upload is.</param>
        /// <param name="name">Every record's name; by default its FASTA header description, else its id.</param>
        /// <param name="dataClass">restricted (the default), public_microbe, public_eukaryote or human.</param>
        /// <param name="lengthBp">ZCWP only, required: the sequence's length, which must agree with the window count.</param>
        /// <param name="recordId">ZCWP only: the record's id within the upload, by default "profile".</param>
        public static async Task<Results<Created<UploadResponse>, ProblemHttpResult>> UploadAsync(HttpRequest request, NpgsqlDataSource dataSource,
            ClaimsPrincipal user, ILoggerFactory loggers, string? kind, string? name, string? dataClass, int? taxid, long? lengthBp, string? recordId,
            CancellationToken cancellationToken)
        {
            if (kind is null || !kinds.Contains(kind))
                throw new ArgumentException($"kind is required: one of {string.Join(", ", kinds)}.");
            dataClass ??= "restricted";
            if (!dataClasses.Contains(dataClass))
                throw new ArgumentException($"dataClass must be one of {string.Join(", ", dataClasses)}.");
            if (name is { Length: > maxNameLength }) throw new ArgumentException($"name is longer than {maxNameLength} characters.");

            var format = RequestBody.Format(request);
            if (format == BodyFormat.Unsupported) return RequestBody.Unsupported();
            List<ProfiledRecord> records;
            byte[] rawSha256;
            long rawByteSize;
            await using (var body = new HashingReadStream(request.Body))
            {
                records = format == BodyFormat.Fasta
                    ? ProfileFasta(request.HttpContext, body)
                    : [await ReadProfileAsync(body, lengthBp, recordId ?? "profile", cancellationToken)];
                rawSha256 = body.Sha256();
                rawByteSize = body.BytesRead;
            }

            var rows = records.Select(r => (r, new RecordRow(r.HeaderId, name ?? (r.Description.Length > 0 ? r.Description[..Math.Min(r.Description.Length, maxNameLength)] : r.HeaderId),
                kind, taxid, new Lineage(null, null, null, null, null, null, null, null, null), "{}", DataClass: dataClass))).ToList();
            StoredSource stored;
            await using (var connection = await dataSource.OpenConnectionAsync(cancellationToken))
            await using (var transaction = await connection.BeginTransactionAsync(cancellationToken))
            {
                stored = await RecordStore.StoreAsync(connection, transaction,
                    new SourceRow("user upload", null, null, rawSha256, rawByteSize, null, user.UserId()), rows);
                await transaction.CommitAsync(cancellationToken);
            }

            // Hash and lengths only: the sequence never reaches a log (S4.4).
            loggers.CreateLogger(typeof(SequenceUpload).FullName!).LogInformation(
                "User {UserId} uploaded source {SourceId} ({Bytes} bytes, sha256 {Sha256}): {Records} record(s), {Bases} bp",
                user.UserId(), stored.SourceId, rawByteSize, Convert.ToHexStringLower(rawSha256), records.Count, records.Sum(r => r.LengthBp));

            var uploaded = records.Zip(stored.SequenceIds, (r, id) => new UploadedRecordResponse(id, r.HeaderId, r.LengthBp,
                r.Levels.Count == 0 ? 0 : r.Levels[0].Level.WindowCount, $"/api/sequences/{id}")).ToList();
            return TypedResults.Created(uploaded[0].Href, new UploadResponse(uploaded));
        }

        private static List<ProfiledRecord> ProfileFasta(HttpContext context, Stream body)
        {
            using var text = RequestBody.OpenFasta(context, body);
            var records = RecordStore.Profile(text).ToList();

            if (records.Count == 0) throw new FormatException("The upload holds no FASTA record (no '>' header line).");
            if (records.FirstOrDefault(r => r.LengthBp == 0) is { } empty) throw new FormatException($"Record {empty.HeaderId} has no bases.");
            if (records.GroupBy(r => r.HeaderId).FirstOrDefault(g => g.Count() > 1) is { } repeated)
                throw new FormatException($"Record id {repeated.Key} appears more than once.");
            return records;
        }

        private static async Task<ProfiledRecord> ReadProfileAsync(Stream body, long? lengthBp, string recordId, CancellationToken cancellationToken)
        {
            if (lengthBp is null) throw new ArgumentException("lengthBp is required with a profile upload.");
            using var bytes = new MemoryStream();
            await body.CopyToAsync(bytes, cancellationToken);
            return RecordStore.FromFineProfile(recordId, ZCurveWindowedProfile.FromBinary(bytes.ToArray()), lengthBp.Value);
        }
    }
}
