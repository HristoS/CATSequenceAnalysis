using System.Text.Json;
using System.Text.RegularExpressions;
using Npgsql;
using SequenceAnaliseTests;
using SequencePlatform.Core;

namespace SequencePlatform.Tests
{
    internal sealed record BundledGenome(string FileName, string Name, string Superkingdom, Lineage Lineage);

    /// <summary>
    /// S2.3's local-FASTA adapter: loads the bundled corpus in data/genomes/ (30 bacteria, 31 viruses) through
    /// <see cref="RecordStore"/>, with the taxonomy the analysis tests already use (<see cref="GenomeTaxonomy"/>,
    /// <see cref="ViralTaxonomy"/>, linked from SequenceAnaliseTests). It is meant for P3's validation (S3.7) on the
    /// throwaway test database; the bundled corpus never goes into the real catalogue. Bacteria carry genus,
    /// phylum and species, and viruses carry genus and family: the only taxonomic claims those lists make. No
    /// taxid was recorded at download, so ncbi_taxid stays null.
    /// </summary>
    internal static class BundledCorpus
    {
        public const string Provider = "NCBI RefSeq";
        private const string license = "INSDC: free and unrestricted";

        // The accession at the end of every bundled file name: an assembly (GCF_000005845.2) or a nucleotide record (NC_001806.2, AC_000008.1).
        private static readonly Regex accessionPattern = new(@"(GC[AF]|[A-Z]{2})_\d+\.\d+(?=\.fna\.gz$)");

        public static IReadOnlyList<BundledGenome> Genomes { get; } =
            GenomeTaxonomy.Genomes
                .Select(g => new BundledGenome(g.FileName, g.OrganismName, "Bacteria",
                    new Lineage("Bacteria", null, null, g.Phylum, null, null, null, g.Genus, string.Join(' ', g.OrganismName.Split(' ').Take(2)))))
                .Concat(ViralTaxonomy.Genomes
                    .Select(g => new BundledGenome(g.FileName, g.Name, "Viruses",
                        new Lineage("Viruses", null, null, null, null, null, g.Family, g.Genus, null))))
                .ToList();

        public static string Accession(string fileName) =>
            accessionPattern.Match(Path.GetFileName(fileName)) is { Success: true } match
                ? match.Value
                : throw new FormatException($"No accession at the end of '{fileName}'.");

        // Where data/genomes/README.md and data/genomes/viral/README.md say each file came from.
        private static string Url(string accession) => accession.StartsWith("GC", StringComparison.Ordinal)
            ? $"https://www.ncbi.nlm.nih.gov/datasets/genome/{accession}/"
            : $"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&rettype=fasta&retmode=text";

        /// <returns>The new sequence_source id.</returns>
        public static async Task<long> LoadAsync(NpgsqlDataSource dataSource, BundledGenome genome)
        {
            var path = GenomeCatalog.PathTo(genome.FileName);
            var accession = Accession(genome.FileName);
            var (sha256, byteSize) = RecordStore.HashFile(path);
            var rows = RecordStore.Profile(path)
                .Select(record => (record, new RecordRow(record.HeaderId, genome.Name, Kind(genome, record.Description), null, genome.Lineage,
                    Annotations(genome, record.Description))))
                .ToList();

            await using var connection = await dataSource.OpenConnectionAsync();
            await using var transaction = await connection.BeginTransactionAsync();
            var stored = await RecordStore.StoreAsync(connection, transaction,
                new SourceRow(Provider, Url(accession), accession, sha256, byteSize, license), rows);
            await transaction.CommitAsync();
            return stored.SourceId;
        }

        // RefSeq titles name every bacterial plasmid ("... plasmid pKPN3, complete sequence"); the chromosome never says it.
        private static string Kind(BundledGenome genome, string description) =>
            genome.Superkingdom == "Viruses" ? "virus"
            : description.Contains("plasmid", StringComparison.OrdinalIgnoreCase) ? "plasmid"
            : "bacterium";

        private static string Annotations(BundledGenome genome, string description) =>
            JsonSerializer.Serialize(new Dictionary<string, object?>
            {
                ["bundled"] = new Dictionary<string, object?> { ["file"] = genome.FileName },
                ["contig"] = new Dictionary<string, object?> { ["description"] = description },
            });
    }
}
