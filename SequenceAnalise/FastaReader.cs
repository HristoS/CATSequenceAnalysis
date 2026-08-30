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
    /// CATProfile/ZCurveWindowedProfile can be built from the downloaded RefSeq assemblies under
    /// data/genomes/ instead of only from Generator's synthetic sequences.
    ///
    /// Several of those real assemblies carry IUPAC ambiguity codes beyond ACGT (N, and in three
    /// cases K/M/R/S/W/Y - confirmed by direct inspection, not assumed), which real sequencing and
    /// assembly data legitimately contains. No filtering or validation of the alphabet happens
    /// here: both Benchmark.NearMatch/ExactMatch (DNA.cs) and ZCurveProfile.Axes already treat an
    /// unrecognized character as a harmless zero-contribution/full-mismatch rather than throwing,
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

        private sealed class FastaSequenceReader : TextReader
        {
            private readonly StreamReader _inner;
            private readonly string _recordId;
            private bool _foundStart;
            private bool _finished;
            private string? _currentLine;
            private int _linePos;

            public FastaSequenceReader(Stream stream, string recordId)
            {
                _inner = new StreamReader(stream);
                _recordId = recordId;
            }

            public override int Read(char[] buffer, int index, int count)
            {
                if (_finished) return 0;

                int written = 0;
                while (written < count)
                {
                    if (_currentLine is null || _linePos >= _currentLine.Length)
                    {
                        if (!AdvanceLine())
                        {
                            _finished = true;
                            break;
                        }
                        continue;
                    }

                    int toCopy = Math.Min(_currentLine.Length - _linePos, count - written);
                    _currentLine.CopyTo(_linePos, buffer, index + written, toCopy);
                    _linePos += toCopy;
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
                    var line = _inner.ReadLine();
                    if (line is null) return false;
                    if (line.Length == 0) continue;

                    if (line[0] == '>')
                    {
                        if (_foundStart) return false; // reached the next record: this one is done

                        var (id, _) = ParseHeader(line);
                        _foundStart = id == _recordId;
                        continue;
                    }

                    if (!_foundStart) continue; // body line of a record we're skipping past

                    _currentLine = line.ToUpperInvariant();
                    _linePos = 0;
                    return true;
                }
            }

            protected override void Dispose(bool disposing)
            {
                if (disposing) _inner.Dispose();
                base.Dispose(disposing);
            }
        }
    }
}
