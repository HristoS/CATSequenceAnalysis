using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequencePlatform.Core;

namespace SequencePlatform.Experiments
{
    /// <summary>
    /// The genomes the P9 experiments read, selected exactly as scripts/platform-checks/shared/ExperimentCorpus.cs selects them, so a port
    /// gives its script's numbers on the same database: public records only, the largest bacterial or archaeal record of each source, every
    /// virus and every eukaryote, with dinucleotide counts; no plasmids, no uploads. A record described as partial is marked Excluded. Only
    /// the provider lineage is read: the ports so far need no NCBI or GTDB lineage.
    /// </summary>
    public static class Corpus
    {
        private static readonly Regex partial = new(@"\bpartial\b", RegexOptions.IgnoreCase | RegexOptions.Compiled);

        /// <param name="Fields">BV-BRC's annotations as text (a string, or else its raw JSON), plus the contig's "description".</param>
        public sealed record Genome(long Id, string Group, string Kind, string Name, long LengthBp, double Gc, double[] Signature, Lineage Provider,
            IReadOnlyDictionary<string, string> Fields, bool Excluded)
        {
            /// <summary>A BV-BRC annotation (host_name, isolation_source, ...) or the contig's "description"; null if missing or empty.</summary>
            public string? Field(string name) => Fields.TryGetValue(name, out var value) && value != "" ? value : null;
        }

        public static async Task<List<Genome>> LoadAsync(NpgsqlDataSource dataSource, CancellationToken cancellationToken)
        {
            var genomes = new List<Genome>();
            await using var command = dataSource.CreateCommand(
                "select r.id, r.kind, coalesce(r.name, r.record_id), r.length_bp, r.dinucleotide_counts, r.annotations::text, " +
                "r.tax_superkingdom, r.tax_realm, r.tax_kingdom, r.tax_phylum, r.tax_class, r.tax_order, r.tax_family, r.tax_genus, r.tax_species " +
                "from (select r.*, row_number() over (partition by r.source_id, r.kind in ('bacterium', 'archaeon') order by r.length_bp desc, r.id) as n " +
                "      from sequence_record r join sequence_source s on s.id = r.source_id" +
                "      where s.owner_id is null and r.dinucleotide_counts is not null and r.kind in ('bacterium', 'archaeon', 'virus', 'eukaryote')) r " +
                "where r.kind = 'virus' or r.n = 1 order by r.id");
            await using var reader = await command.ExecuteReaderAsync(cancellationToken);
            while (await reader.ReadAsync(cancellationToken))
            {
                string? Text(int i) => reader.IsDBNull(i) || reader.GetString(i) == "" ? null : reader.GetString(i);
                var kind = reader.GetString(1);
                var counts = reader.GetFieldValue<long[]>(4);
                var fields = Flatten(Text(5));
                double gc = Enumerable.Range(4, 8).Sum(i => (double)counts[i]) / counts.Sum();   // first base C or G: indexes 4..11
                genomes.Add(new Genome(reader.GetInt64(0), SignatureIndex.GroupOf(kind)!, kind, reader.GetString(2), reader.GetInt64(3), gc,
                    SignatureIndex.Signature(counts), new Lineage(Text(6), Text(7), Text(8), Text(9), Text(10), Text(11), Text(12), Text(13), Text(14)),
                    fields, fields.TryGetValue("description", out var description) && partial.IsMatch(description)));
            }
            return genomes;
        }

        private static Dictionary<string, string> Flatten(string? annotations)
        {
            var fields = new Dictionary<string, string>(StringComparer.Ordinal);
            if (annotations is null) return fields;
            using var document = JsonDocument.Parse(annotations);
            if (document.RootElement.TryGetProperty("bvbrc", out var bvbrc))
                foreach (var property in bvbrc.EnumerateObject())
                    fields[property.Name] = property.Value.ValueKind == JsonValueKind.String ? property.Value.GetString()! : property.Value.GetRawText();
            if (document.RootElement.TryGetProperty("contig", out var contig) && contig.TryGetProperty("description", out var description)
                && description.ValueKind == JsonValueKind.String)
                fields["description"] = description.GetString()!;
            return fields;
        }
    }
}
