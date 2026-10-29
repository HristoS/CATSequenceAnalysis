using System.Collections.Concurrent;
using Microsoft.Extensions.Logging;

namespace SequencePlatform.Tests
{
    /// <summary>Keeps every log entry the API writes as text: the formatted message, each structured value, and the exception.</summary>
    internal sealed class LogCapture : ILoggerProvider
    {
        private readonly ConcurrentQueue<string> entries = new();

        public IReadOnlyList<string> Entries => entries.ToList();

        public ILogger CreateLogger(string categoryName) => new Logger(categoryName, entries);

        public void Dispose() { }

        private sealed class Logger(string category, ConcurrentQueue<string> entries) : ILogger
        {
            public IDisposable? BeginScope<TState>(TState state) where TState : notnull
            {
                entries.Enqueue($"{category} scope: {Values(state)}");
                return null;
            }

            public bool IsEnabled(LogLevel logLevel) => true;

            public void Log<TState>(LogLevel logLevel, EventId eventId, TState state, Exception? exception, Func<TState, Exception?, string> formatter) =>
                entries.Enqueue($"{logLevel} {category}: {formatter(state, exception)} | {Values(state)} | {exception}");

            private static string Values(object? state) =>
                state is IEnumerable<KeyValuePair<string, object?>> pairs ? string.Join("; ", pairs.Select(p => $"{p.Key}={p.Value}")) : $"{state}";
        }
    }
}
