using System.Globalization;
using System.Text.Json;
using System.Text.Json.Nodes;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Ingest
{
    /// <summary>
    /// publish-experiments: every experiment an article states (results/experiments.json) published as an 'article' run with its CSVs
    /// (migration 0023), so the experiments page shows each with its numbers and its files for download. The CSVs come from the
    /// registry's folder; nothing is recomputed. A run published again replaces the previous one of the same name and origin.
    /// </summary>
    internal static class ExperimentPublisher
    {
        public static async Task<int> PublishAsync(NpgsqlDataSource dataSource, string registryPath, string? onlyName)
        {
            var registry = JsonNode.Parse(await File.ReadAllTextAsync(registryPath))!.AsObject();
            var folder = Path.GetDirectoryName(Path.GetFullPath(registryPath))!;
            var defaults = registry["defaults"]!.AsObject();
            var entries = registry["experiments"]!.AsArray().Select(e => e!.AsObject())
                .Where(e => onlyName is null || (string)e["name"]! == onlyName).ToList();
            if (entries.Count == 0)
            {
                Console.WriteLine(onlyName is null ? "The registry lists no experiment." : $"No experiment '{onlyName}' in the registry.");
                return 1;
            }

            // Check every entry before publishing any, so a typo doesn't leave the page half updated.
            var problems = entries.SelectMany(e => Problems(e, folder)).ToList();
            if (problems.Count > 0)
            {
                foreach (var problem in problems) Console.WriteLine($"  {problem}");
                return 1;
            }

            foreach (var entry in entries)
            {
                string name = (string)entry["name"]!;
                string? test = (string?)entry["test"], script = (string?)entry["script"];
                string source = test is not null ? $"SequenceAnaliseTests/{test}.cs"
                    : script!.EndsWith(".cs", StringComparison.Ordinal) && !script.Contains('/') ? $"scripts/platform-checks/{script}" : script!;
                string command = (string?)entry["command"] ?? (test is not null
                    ? string.Format(CultureInfo.InvariantCulture, (string)defaults["libraryCommand"]!, test)
                    : string.Format(CultureInfo.InvariantCulture, (string)defaults["scriptCommand"]!, script));
                var summary = new JsonObject();
                foreach (var key in new[] { "section", "rule", "note" })
                    if (entry[key] is { } value) summary[key] = value.DeepClone();
                if (entry["keyNumbers"] is { } numbers) summary["keyNumbers"] = numbers.DeepClone();
                var details = new ExperimentDetails((string?)entry["article"], ExperimentStore.ArticleOrigin, (string?)entry["commit"],
                    (string?)entry["machine"] ?? (string?)defaults["machine"], command, (bool?)entry["runnable"] ?? false);
                var files = entry["files"]!.AsArray().Select(f => (string)f!)
                    .Select(f => new ExperimentFile(f, File.ReadAllBytes(Path.Combine(folder, f)))).ToList();

                await ExperimentStore.PublishAsync(dataSource, name, (string)entry["title"]!, (string)entry["question"]!, (string)entry["verdict"]!,
                    (string)entry["snapshot"]!, source, summary.ToJsonString(), [], details, files);
                Console.WriteLine($"  {name}: {files.Count} file(s), {files.Sum(f => f.Content.Length):N0} bytes");
            }
            Console.WriteLine($"{entries.Count} experiment(s) published from {registryPath}.");
            return 0;
        }

        private static IEnumerable<string> Problems(JsonObject entry, string folder)
        {
            var name = (string?)entry["name"] ?? "(no name)";
            foreach (var key in new[] { "name", "title", "question", "verdict", "snapshot", "files" })
                if (entry[key] is null) yield return $"{name}: no '{key}'";
            if ((string?)entry["verdict"] is { } verdict && verdict is not (ExperimentStore.Pass or ExperimentStore.Fail or ExperimentStore.Context))
                yield return $"{name}: verdict '{verdict}' is not pass, fail or context";
            if ((entry["test"] is null) == (entry["script"] is null)) yield return $"{name}: exactly one of 'test' and 'script'";
            if (entry["keyNumbers"] is JsonArray numbers && numbers.Any(n => n?["value"]?.GetValueKind() != JsonValueKind.Number))
                yield return $"{name}: every key number needs a numeric 'value'";
            foreach (var file in entry["files"]?.AsArray() ?? [])
                if (!File.Exists(Path.Combine(folder, (string)file!))) yield return $"{name}: missing file {file}";
        }
    }
}
