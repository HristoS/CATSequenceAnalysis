using System.Net;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// GETs started at least <c>minRequestInterval</c> apart, also when callers run in parallel, with 429/5xx and
    /// network errors retried on exponential backoff that honours Retry-After. Shared by the BV-BRC, NCBI and GTDB
    /// clients, none of which we want to load harder than a careful person would.
    /// </summary>
    internal sealed class PoliteHttpClient : IDisposable
    {
        private const int maxAttempts = 5;

        private readonly HttpClient http;
        private readonly string serviceName;
        private readonly TimeSpan minRequestInterval;
        private readonly Lock slotLock = new();
        private DateTime nextSlotUtc = DateTime.MinValue;

        public PoliteHttpClient(string serviceName, string? contactEmail, TimeSpan minRequestInterval)
        {
            this.serviceName = serviceName;
            this.minRequestInterval = minRequestInterval;
            http = new HttpClient { Timeout = TimeSpan.FromMinutes(10) };
            var agent = "SequencePlatform.Ingest/0.1 (research test" + (contactEmail is null ? "" : $"; mailto:{contactEmail}") + ")";
            http.DefaultRequestHeaders.UserAgent.ParseAdd(agent);
        }

        public async Task<HttpResponseMessage> SendAsync(string url, string accept, HttpCompletionOption completion = HttpCompletionOption.ResponseContentRead)
        {
            for (int attempt = 1; ; attempt++)
            {
                var now = DateTime.UtcNow;
                var wait = ReserveSlot(now) - now;
                if (wait > TimeSpan.Zero) await Task.Delay(wait);

                using var request = new HttpRequestMessage(HttpMethod.Get, url);
                request.Headers.Accept.ParseAdd(accept);

                HttpResponseMessage? response = null;
                try
                {
                    response = await http.SendAsync(request, completion);
                    if (response.IsSuccessStatusCode) return response;
                    // Other 4xx won't succeed on retry; a distinct exception type keeps it out of the retry catch below.
                    if (response.StatusCode != HttpStatusCode.TooManyRequests && (int)response.StatusCode < 500)
                        throw new InvalidOperationException($"{serviceName} returned {(int)response.StatusCode} for {url}");
                }
                catch (Exception e) when (e is HttpRequestException or TaskCanceledException && attempt < maxAttempts)
                {
                    Console.WriteLine($"  request failed ({e.Message}), retrying");
                }

                if (attempt >= maxAttempts)
                    throw new HttpRequestException($"{serviceName} request failed after {maxAttempts} attempts: {url}");

                var backoff = response?.Headers.RetryAfter?.Delta ?? TimeSpan.FromSeconds(Math.Pow(2, attempt));
                response?.Dispose();
                await Task.Delay(backoff);
            }
        }

        // Each caller takes its own slot before waiting, so parallel callers queue up rather than all waking at the same moment.
        internal DateTime ReserveSlot(DateTime nowUtc)
        {
            lock (slotLock)
            {
                var slot = nextSlotUtc > nowUtc ? nextSlotUtc : nowUtc;
                nextSlotUtc = slot + minRequestInterval;
                return slot;
            }
        }

        public void Dispose() => http.Dispose();
    }
}
