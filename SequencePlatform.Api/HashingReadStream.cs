using System.Security.Cryptography;

namespace SequencePlatform.Api
{
    /// <summary>
    /// Passes a request body through while hashing it, so an upload's raw_sha256 and raw_byte_size come from the same single
    /// read that profiles it; the body is never kept.
    /// </summary>
    internal sealed class HashingReadStream(Stream inner) : Stream
    {
        private readonly IncrementalHash sha256 = IncrementalHash.CreateHash(HashAlgorithmName.SHA256);

        public long BytesRead { get; private set; }

        /// <summary>Call once the body has been read to its end.</summary>
        public byte[] Sha256() => sha256.GetHashAndReset();

        public override int Read(byte[] buffer, int offset, int count) => Read(buffer.AsSpan(offset, count));

        public override int Read(Span<byte> buffer)
        {
            int read = inner.Read(buffer);
            Count(buffer[..read]);
            return read;
        }

        public override async ValueTask<int> ReadAsync(Memory<byte> buffer, CancellationToken cancellationToken = default)
        {
            int read = await inner.ReadAsync(buffer, cancellationToken);
            Count(buffer.Span[..read]);
            return read;
        }

        public override Task<int> ReadAsync(byte[] buffer, int offset, int count, CancellationToken cancellationToken) =>
            ReadAsync(buffer.AsMemory(offset, count), cancellationToken).AsTask();

        private void Count(ReadOnlySpan<byte> chunk)
        {
            sha256.AppendData(chunk);
            BytesRead += chunk.Length;
        }

        public override bool CanRead => true;
        public override bool CanSeek => false;
        public override bool CanWrite => false;
        public override long Length => throw new NotSupportedException();
        public override long Position { get => BytesRead; set => throw new NotSupportedException(); }
        public override void Flush() { }
        public override long Seek(long offset, SeekOrigin origin) => throw new NotSupportedException();
        public override void SetLength(long value) => throw new NotSupportedException();
        public override void Write(byte[] buffer, int offset, int count) => throw new NotSupportedException();

        protected override void Dispose(bool disposing)
        {
            if (disposing) sha256.Dispose();
            base.Dispose(disposing);
        }
    }
}
