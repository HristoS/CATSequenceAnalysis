using System.IO.Compression;
using System.Text;

namespace SequenceAnalyses
{
    /// <summary>
    /// One FASTA record: the header line's id (everything up to the first space) and the rest of
    /// the description, plus the full sequence concatenated across all wrapped lines and
    /// uppercased (the FASTA spec permits lowercase for soft-masked/repeat regions; several of the
    /// real assemblies under data/genomes/ are otherwise all-uppercase already, but this is not
    /// something to assume of arbitrary FASTA input).
    /// </summary>
    public sealed record FastaRecord(string Id, string Description, string Sequence)
    {
        public int Length => Sequence.Length;
    }

    /// <summary>
    /// Minimal FASTA parser (transparently handles ".gz") for real reference-genome files, so
    /// ZCurveWindowedProfile/WindowedProfile can be built from the downloaded RefSeq assemblies under
    /// data/genomes/ instead of only from Generator's synthetic sequences.
    ///
    /// Several of those real assemblies carry IUPAC ambiguity codes beyond ACGT (N, and in three
    /// cases K/M/R/S/W/Y - confirmed by direct inspection, not assumed), which real sequencing and
    /// assembly data legitimately contains. No filtering or validation of the alphabet happens
    /// here: ZCurveProfile.Axes treats an unrecognized character as a harmless zero contribution,
    /// and WindowedProfile restarts its rolling k-mer at one, rather than either throwing,
    /// so passing these characters through unchanged is safe for every profile type in this
    /// project.
    /// </summary>
    public static class FastaReader
    {
        private static Stream OpenRawStream(string path)
        {
            Stream stream = File.OpenRead(path);
            if (path.EndsWith(".gz", StringComparison.OrdinalIgnoreCase))
                stream = new GZipStream(stream, CompressionMode.Decompress);
            return stream;
        }

        private static (string Id, string Description) ParseHeader(string headerLine)
        {
            var header = headerLine.Substring(1); // drop the leading '>'
            var spaceIndex = header.IndexOf(' ');
            return spaceIndex < 0
                ? (header, string.Empty)
                : (header.Substring(0, spaceIndex), header.Substring(spaceIndex + 1));
        }

        /// <summary>
        /// Lazily yields every record in the file, each with its full sequence materialized as one
        /// uppercase string. Fine at bacterial-genome scale (data/genomes/'s largest single record
        /// is ~8.7 Mbp); for sequences too large to comfortably hold as a string - see CLAUDE.md's
        /// scale discussion - use <see cref="OpenSequenceReader"/> instead, which never
        /// materializes a full sequence.
        /// </summary>
        public static IEnumerable<FastaRecord> ReadRecords(string path)
        {
            using var stream = OpenRawStream(path);
            using var reader = new StreamReader(stream);

            string? id = null;
            string description = string.Empty;
            StringBuilder? sequence = null;

            string? line;
            while ((line = reader.ReadLine()) != null)
            {
                if (line.Length == 0) continue;

                if (line[0] == '>')
                {
                    if (id is not null)
                        yield return new FastaRecord(id, description, sequence!.ToString());

                    (id, description) = ParseHeader(line);
                    sequence = new StringBuilder();
                }
                else
                {
                    if (sequence is null)
                        throw new FormatException($"'{path}' has sequence data before its first '>' header line.");
                    sequence.Append(line.ToUpperInvariant());
                }
            }

            if (id is not null)
                yield return new FastaRecord(id, description, sequence!.ToString());
        }

        /// <summary>
        /// The single longest record in the file - the chromosome, in an assembly that also
        /// carries plasmids/contigs (RefSeq bacterial assemblies list the chromosome first, but
        /// this picks the true longest record rather than assuming that ordering). Reads the whole
        /// file once, holding at most two candidate records' sequences in memory at a time rather
        /// than every record's.
        /// </summary>
        public static FastaRecord ReadPrimaryRecord(string path)
        {
            FastaRecord? best = null;
            foreach (var record in ReadRecords(path))
                if (best is null || record.Length > best.Length)
                    best = record;

            return best ?? throw new FormatException($"'{path}' contains no FASTA records.");
        }

        /// <summary>
        /// Every record's header (id and description, as <see cref="ReadRecords"/> parses them) in file
        /// order, without keeping any sequence: one pass that only looks at header lines. Each id can be
        /// handed to <see cref="OpenSequenceReader"/> to stream that record.
        /// </summary>
        public static IEnumerable<(string Id, string Description)> ReadHeaders(string path)
        {
            using var stream = OpenRawStream(path);
            using var reader = new StreamReader(stream);
            string? line;
            while ((line = reader.ReadLine()) != null)
                if (line.Length > 0 && line[0] == '>')
                    yield return ParseHeader(line);
        }

        /// <summary>
        /// Opens a streaming <see cref="TextReader"/> over one record's bases only (header lines
        /// and every other record's body are skipped without being buffered), suitable for
        /// <see cref="ZCurveWindowedProfile.FromStream"/> - never materializes the sequence as a
        /// single string, so this is the entry point for sequences too large for that (see
        /// <see cref="ReadRecords"/>'s doc comment). Requires an explicit <paramref name="recordId"/>
        /// (from <see cref="ReadRecords"/> or <see cref="ReadPrimaryRecord"/>) rather than
        /// defaulting to "the longest record", since determining that would itself require a full
        /// pass over the file - defeating the point of a non-materializing reader.
        /// </summary>
        public static TextReader OpenSequenceReader(string path, string recordId)
        {
            if (recordId is null) throw new ArgumentNullException(nameof(recordId));
            return new FastaSequenceReader(OpenRawStream(path), recordId);
        }

        /// <summary>
        /// Every record of FASTA text that can be read only once (an upload's request body), in order, each as a streaming
        /// reader over its bases - uppercased, line breaks dropped, exactly what <see cref="OpenSequenceReader"/> gives for a
        /// file. Nothing is buffered beyond a fixed-size block, not even a line, so a single-line record of any length streams
        /// too. Read each record's reader before moving on: advancing skips whatever is left of it, and the old reader then
        /// reads nothing.
        /// </summary>
        public static IEnumerable<(string Id, string Description, TextReader Sequence)> ReadRecordStreams(TextReader reader)
        {
            if (reader is null) throw new ArgumentNullException(nameof(reader));
            var walker = new FastaRecordWalker(reader);
            while (walker.NextHeader() is { } header)
            {
                var (id, description) = ParseHeader(header);
                yield return (id, description, new RecordBodyReader(walker, walker.Record));
            }
        }

        // One pass over FASTA text, character by character through a fixed buffer. A '>' at the start of a line opens a
        // record; everything else up to the next one is that record's body.
        private sealed class FastaRecordWalker
        {
            private const int maxHeaderLength = 65536;

            private readonly TextReader inner;
            private readonly char[] buffer = new char[65536];
            private int position;
            private int length;
            private bool atLineStart = true;

            public FastaRecordWalker(TextReader inner) => this.inner = inner;

            /// <summary>Counts the headers read; a body reader reads only while its record is the current one.</summary>
            public int Record { get; private set; }

            /// <returns>The next header line, '>' included; null at the end. Skips what is left of the current body.</returns>
            public string? NextHeader()
            {
                while (Fill())
                {
                    char c = buffer[position++];
                    if (atLineStart && c == '>')
                    {
                        Record++;
                        return ReadHeaderRest();
                    }
                    if (Record == 0 && !char.IsWhiteSpace(c))
                        throw new FormatException("FASTA text has sequence data before its first '>' header line.");
                    atLineStart = c == '\n';
                }
                return null;
            }

            private string ReadHeaderRest()
            {
                var header = new StringBuilder(">");
                while (Fill())
                {
                    char c = buffer[position++];
                    if (c == '\n') break;
                    if (c == '\r') continue;
                    if (header.Length == maxHeaderLength) throw new FormatException($"A FASTA header line is longer than {maxHeaderLength} characters.");
                    header.Append(c);
                }
                atLineStart = true;
                return header.ToString();
            }

            /// <returns>Body characters of the current record copied into <paramref name="destination"/>; 0 at its end.</returns>
            public int ReadBody(Span<char> destination)
            {
                int written = 0;
                while (written < destination.Length && Fill())
                {
                    char c = buffer[position];
                    if (atLineStart && c == '>') break; // the next record's header: left for NextHeader
                    position++;
                    if (c == '\n') atLineStart = true;
                    else if (c != '\r')
                    {
                        atLineStart = false;
                        destination[written++] = char.ToUpperInvariant(c);
                    }
                }
                return written;
            }

            private bool Fill()
            {
                if (position < length) return true;
                length = inner.Read(buffer, 0, buffer.Length);
                position = 0;
                return length > 0;
            }
        }

        private sealed class RecordBodyReader : TextReader
        {
            private readonly FastaRecordWalker walker;
            private readonly int record;

            public RecordBodyReader(FastaRecordWalker walker, int record)
            {
                this.walker = walker;
                this.record = record;
            }

            public override int Read(char[] buffer, int index, int count) => Read(buffer.AsSpan(index, count));

            public override int Read(Span<char> buffer) => walker.Record == record ? walker.ReadBody(buffer) : 0;

            public override int Read()
            {
                Span<char> one = stackalloc char[1];
                return Read(one) == 1 ? one[0] : -1;
            }
        }

        private sealed class FastaSequenceReader : TextReader
        {
            private readonly StreamReader inner;
            private readonly string recordId;
            private bool foundStart;
            private bool finished;
            private string? currentLine;
            private int linePos;

            public FastaSequenceReader(Stream stream, string recordId)
            {
                inner = new StreamReader(stream);
                this.recordId = recordId;
            }

            public override int Read(char[] buffer, int index, int count)
            {
                if (finished) return 0;

                int written = 0;
                while (written < count)
                {
                    if (currentLine is null || linePos >= currentLine.Length)
                    {
                        if (!AdvanceLine())
                        {
                            finished = true;
                            break;
                        }
                        continue;
                    }

                    int toCopy = Math.Min(currentLine.Length - linePos, count - written);
                    currentLine.CopyTo(linePos, buffer, index + written, toCopy);
                    linePos += toCopy;
                    written += toCopy;
                }

                return written;
            }

            // Advances to the next sequence line belonging to the target record. Returns false
            // once that record's sequence is exhausted (next '>' header, or end of file).
            private bool AdvanceLine()
            {
                while (true)
                {
                    var line = inner.ReadLine();
                    if (line is null) return false;
                    if (line.Length == 0) continue;

                    if (line[0] == '>')
                    {
                        if (foundStart) return false; // reached the next record: this one is done

                        var (id, _) = ParseHeader(line);
                        foundStart = id == recordId;
                        continue;
                    }

                    if (!foundStart) continue; // body line of a record we're skipping past

                    currentLine = line.ToUpperInvariant();
                    linePos = 0;
                    return true;
                }
            }

            protected override void Dispose(bool disposing)
            {
                if (disposing) inner.Dispose();
                base.Dispose(disposing);
            }
        }
    }
}
