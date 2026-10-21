using SequencePlatform.Core;
using SequencePlatform.Ingest;

namespace SequencePlatform.Tests
{
    public class TaxonomyTests
    {
        // Trimmed from a real NCBI Datasets v2 taxonomy report (2026-10-21): an unknown id, HCV's old species-level
        // taxid 11103 (merged into 3052230, Orthohepacivirus hominis) and E. coli 562 (current).
        private const string datasetReport = """
            {"reports": [
              {"query": ["999999999"], "errors": [{"reason": "not a recognized NCBI Taxonomy name", "invalid_identifiers": ["999999999"]}]},
              {"query": ["11103"], "taxonomy": {"tax_id": 3052230, "rank": "SPECIES",
                "current_scientific_name": {"name": "Orthohepacivirus hominis"},
                "classification": {"acellular_root": {"name": "Viruses"}, "realm": {"name": "Riboviria"}, "kingdom": {"name": "Orthornavirae"},
                  "phylum": {"name": "Kitrinoviricota"}, "class": {"name": "Flasuviricetes"}, "order": {"name": "Amarillovirales"},
                  "family": {"name": "Hepaciviridae"}, "genus": {"name": "Orthohepacivirus"}, "species": {"name": "Orthohepacivirus hominis"}}}},
              {"query": ["562"], "taxonomy": {"tax_id": 562, "rank": "SPECIES",
                "current_scientific_name": {"name": "Escherichia coli"},
                "classification": {"domain": {"name": "Bacteria"}, "kingdom": {"name": "Pseudomonadati"}, "phylum": {"name": "Pseudomonadota"},
                  "class": {"name": "Gammaproteobacteria"}, "order": {"name": "Enterobacterales"}, "family": {"name": "Enterobacteriaceae"},
                  "genus": {"name": "Escherichia"}, "species": {"name": "Escherichia coli"}}}}
            ]}
            """;

        [Test]
        public void DatasetReportGivesCurrentMergedAndUnknownTaxa()
        {
            var taxa = TaxonomyLoader.ParseDatasetReport(datasetReport).ToDictionary(t => t.QueriedTaxid);

            Assert.That(taxa[999999999].Status, Is.EqualTo(TaxonomyLoader.Unknown));
            Assert.That(taxa[999999999].CurrentTaxid, Is.Null);

            var hcv = taxa[11103];
            Assert.That((hcv.Status, hcv.CurrentTaxid, hcv.Name, hcv.Rank), Is.EqualTo((TaxonomyLoader.Merged, (int?)3052230, "Orthohepacivirus hominis", "species")));
            Assert.That(hcv.Lineage, Is.EqualTo(new Lineage("Viruses", "Riboviria", "Orthornavirae", "Kitrinoviricota", "Flasuviricetes",
                "Amarillovirales", "Hepaciviridae", "Orthohepacivirus", "Orthohepacivirus hominis")));

            var ecoli = taxa[562];
            Assert.That((ecoli.Status, ecoli.CurrentTaxid), Is.EqualTo((TaxonomyLoader.Current, (int?)562)));
            Assert.That((ecoli.Lineage.Superkingdom, ecoli.Lineage.Realm, ecoli.Lineage.Family), Is.EqualTo(("Bacteria", (string?)null, "Enterobacteriaceae")));
        }

        [Test]
        public void GtdbLinesAndAssemblyAccessionsShareOneKey()
        {
            var taxon = TaxonomyLoader.ParseGtdbLine(
                "GB_GCA_048462655.1\td__Archaea;p__Methanobacteriota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae;g__;s__");

            Assert.That(taxon, Is.EqualTo(new GtdbTaxon("048462655.1", "GB_GCA_048462655.1", "Archaea", "Methanobacteriota",
                "Methanobacteria", "Methanobacteriales", "Methanobacteriaceae", null, null)));
            Assert.That(TaxonomyLoader.AssemblyKey("GCF_048462655.1"), Is.EqualTo(taxon!.AssemblyKey));
            Assert.That(TaxonomyLoader.AssemblyKey("NC_045512.2"), Is.Null);
            Assert.That(TaxonomyLoader.ParseGtdbLine("not a taxonomy line"), Is.Null);
        }

        [Test]
        public async Task LineageDifferenceComparesTheStoredLineageWithNcbis()
        {
            var dataSource = ThrowawayDatabase.DataSource;
            var id = await ThrowawayDatabase.InsertRecordAsync("test-taxonomy", "TAXON_1", 9_646);
            await using (var command = dataSource.CreateCommand(
                "update sequence_record set ncbi_taxid = 11103, tax_superkingdom = 'Viruses', tax_family = 'Flaviviridae', " +
                "tax_genus = 'Hepacivirus' where id = $1"))
            {
                command.Parameters.AddWithValue(id);
                await command.ExecuteNonQueryAsync();
            }
            await TaxonomyLoader.UpsertTaxaAsync(dataSource, TaxonomyLoader.ParseDatasetReport(datasetReport).ToList());

            var rows = new List<(string Rank, string Outcome, string? Provider, string? Ncbi)>();
            await using (var command = dataSource.CreateCommand(
                "select rank, outcome, provider_value, ncbi_value from lineage_difference where sequence_id = $1 order by rank"))
            {
                command.Parameters.AddWithValue(id);
                await using var reader = await command.ExecuteReaderAsync();
                while (await reader.ReadAsync())
                    rows.Add((reader.GetString(0), reader.GetString(1), reader.IsDBNull(2) ? null : reader.GetString(2),
                        reader.IsDBNull(3) ? null : reader.GetString(3)));
            }

            Assert.That(rows.Single(r => r.Rank == "family"), Is.EqualTo(("family", "differs", (string?)"Flaviviridae", (string?)"Hepaciviridae")));
            Assert.That(rows.Single(r => r.Rank == "genus"), Is.EqualTo(("genus", "differs", (string?)"Hepacivirus", (string?)"Orthohepacivirus")));
            Assert.That(rows.Single(r => r.Rank == "order").Outcome, Is.EqualTo("ncbi_only"));
            Assert.That(rows.Any(r => r.Rank == "superkingdom"), Is.False);

            // A second load of the same report updates in place rather than failing on the primary key.
            await TaxonomyLoader.UpsertTaxaAsync(dataSource, TaxonomyLoader.ParseDatasetReport(datasetReport).ToList());
            await using var count = dataSource.CreateCommand("select count(*) from taxon where taxid in (11103, 562, 999999999)");
            Assert.That((long)(await count.ExecuteScalarAsync())!, Is.EqualTo(3));
        }
    }
}
