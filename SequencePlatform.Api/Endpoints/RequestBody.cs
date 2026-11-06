using System.Buffers.Binary;
using Microsoft.AspNetCore.Http.Features;
using Microsoft.AspNetCore.Http.HttpResults;
using Microsoft.Net.Http.Headers;
using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    internal enum BodyFormat { Fasta, Profile, Unsupported }

    /// <summary>
    /// Request bodies that carry sequence data, for uploads and searches alike: the raw file, FASTA or ZCWP, never a multipart form
    /// (IFormFile buffers anything over 64 KB to a temp file).
    /// </summary>
    internal static class RequestBody
    {
        public const string FastaContentType = "text/x-fasta";
        public const string ProfileContentType = "application/x-zcwp";

        // The ZCWP header (docs/database-decision.md section 4.3): magic, version at byte 4, reserved, window size, then the window
        // count at byte 12; x, y and z per window follow, 12 bytes as float32 (version 1) or 6 as float16 (version 2).
        private const int profileHeaderSize = 16;
        private const int profileVersionOffset = 4;
        private const int profileCountOffset = 12;
        // Far over the longest known protein (titin, about 35,000 aa), so only a wrong file reaches it.
        private const int maxProteinChars = 1_000_000;

        private static readonly HashSet<string> fastaTypes = new(StringComparer.OrdinalIgnoreCase) { FastaContentType, "application/x-fasta", "text/plain" };

        public static BodyFormat Format(HttpRequest request)
        {
            var mediaType = MediaTypeHeaderValue.TryParse(request.ContentType, out var parsed) ? parsed.MediaType.Value : null;
            if (mediaType is null) return BodyFormat.Unsupported;
            if (fastaTypes.Contains(mediaType)) return BodyFormat.Fasta;
            return string.Equals(mediaType, ProfileContentType, StringComparison.OrdinalIgnoreCase) ? BodyFormat.Profile : BodyFormat.Unsupported;
        }

        public static ProblemHttpResult Unsupported() =>
            TypedResults.Problem($"Send FASTA as {FastaContentType} or a profile as {ProfileContentType}.", statusCode: StatusCodes.Status415UnsupportedMediaType);

        /// <summary>FASTA text over the body. FromStream reads through a TextReader, which is synchronous; the request size limit bounds the read.</summary>
        public static StreamReader OpenFasta(HttpContext context, Stream body)
        {
            context.Features.GetRequiredFeature<IHttpBodyControlFeature>().AllowSynchronousIO = true;
            return new StreamReader(body, leaveOpen: true);
        }

        /// <summary>
        /// A search query: one FASTA record, profiled on both strands from its DNA, which is then dropped (it lives only for this
        /// request, bounded by the size limit); or the 500 bp ZCWP level of the forward strand, optionally followed directly by the
        /// reverse strand's (the Ingest's profile command writes X.zcwp and X.rc.zcwp; send them concatenated). A profile's length is
        /// windows × 500 unless <paramref name="lengthBp"/> gives it exactly.
        /// </summary>
        public static async Task<QueryProfile> ReadQueryAsync(HttpRequest request, BodyFormat format, long? lengthBp, CancellationToken cancellationToken)
        {
            if (format == BodyFormat.Fasta)
            {
                using var text = OpenFasta(request.HttpContext, request.Body);
                using var records = FastaReader.ReadRecordStreams(text).GetEnumerator();
                if (!records.MoveNext()) throw new FormatException("The query holds no FASTA record (no '>' header line).");
                var sequence = records.Current.Sequence.ReadToEnd();
                if (records.MoveNext()) throw new FormatException("A search takes one FASTA record; this query has more.");
                if (sequence.Length == 0) throw new FormatException($"Record {records.Current.Id} has no bases.");
                return QueryProfile.FromSequence(sequence);
            }

            using var bytes = new MemoryStream();
            await request.Body.CopyToAsync(bytes, cancellationToken);
            var strands = SplitProfiles(bytes.ToArray());
            if (strands.Count is 0 or > 2) throw new FormatException("A profile query is one ZCWP profile, or the forward then the reverse strand.");
            var forward = strands[0];
            var reverse = strands.Count == 2 ? strands[1] : null;
            if (strands.Any(s => s.WindowSize != RecordStore.FineWindowSize))
                throw new ArgumentException($"A profile query is the {RecordStore.FineWindowSize} bp level.");
            if (reverse is not null && reverse.WindowCount != forward.WindowCount)
                throw new ArgumentException($"The reverse strand has {reverse.WindowCount} windows, the forward {forward.WindowCount}.");
            if (lengthBp is long length && (length < 1 || length / RecordStore.FineWindowSize != forward.WindowCount))
                throw new ArgumentException($"A length of {length} bp gives {Math.Max(length, 0) / RecordStore.FineWindowSize} windows of " +
                    $"{RecordStore.FineWindowSize} bp, but the profile has {forward.WindowCount}.");
            return QueryProfile.FromFineProfiles(forward, reverse, lengthBp ?? (long)forward.WindowCount * forward.WindowSize, null, null);
        }

        /// <summary>
        /// A protein search's query: one FASTA record, or the bare letters with no header. Returns the letters as sent, white space
        /// included; <see cref="ProteinQuery.FromProtein"/> checks them. Never echoed in an error.
        /// </summary>
        public static async Task<string> ReadProteinAsync(HttpRequest request, CancellationToken cancellationToken)
        {
            using var reader = new StreamReader(request.Body, leaveOpen: true);
            var buffer = new char[maxProteinChars + 1];
            int length = 0, read;
            while (length < buffer.Length && (read = await reader.ReadAsync(buffer.AsMemory(length), cancellationToken)) > 0) length += read;
            if (length > maxProteinChars) throw new ArgumentException($"A protein query is at most {maxProteinChars:N0} characters.");

            var lines = new string(buffer, 0, length).Split('\n');
            int headers = lines.Count(l => l.StartsWith('>'));
            if (headers > 1) throw new FormatException("A protein search takes one FASTA record; this query has more.");
            if (headers == 1 && !lines.First(l => l.Trim().Length > 0).StartsWith('>'))
                throw new FormatException("A FASTA query starts with its '>' header line.");
            return string.Concat(lines.Where(l => !l.StartsWith('>') && !l.StartsWith(';')));
        }

        private static List<ZCurveWindowedProfile> SplitProfiles(byte[] bytes)
        {
            var profiles = new List<ZCurveWindowedProfile>();
            int offset = 0;
            while (offset < bytes.Length)
            {
                if (bytes.Length - offset < profileHeaderSize) throw new FormatException("Not a ZCWP profile.");
                long count = BinaryPrimitives.ReadInt32LittleEndian(bytes.AsSpan(offset + profileCountOffset));
                int bytesPerWindow = BinaryPrimitives.ReadUInt16LittleEndian(bytes.AsSpan(offset + profileVersionOffset)) == 2 ? 6 : 12;
                long size = profileHeaderSize + bytesPerWindow * count;
                if (count < 0 || size > bytes.Length - offset) throw new FormatException("ZCWP header does not match the payload length.");
                profiles.Add(ZCurveWindowedProfile.FromBinary(bytes.AsSpan(offset, (int)size)));
                offset += (int)size;
            }
            return profiles;
        }
    }
}
