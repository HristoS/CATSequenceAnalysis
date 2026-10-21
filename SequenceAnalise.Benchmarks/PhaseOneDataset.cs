using System.IO.Compression;
using System.Security.Cryptography;
using System.Text;
using System.Text.Json;
using System.Text.RegularExpressions;
using SequenceAnalyses;

namespace SequenceAnalise.Benchmarks
{
    public sealed record QueryManifest(string QueryId, string Stratum, int LengthBp, string SequenceSha256);

    public sealed record LocusManifest(string LocusId, int StartBp, QueryManifest[] Queries);

    public sealed record ReferenceManifest(
        string FileName,
        string Accession,
        string RecordId,
        int RecordLengthBp,
        string CompressedSha256,
        int SamplingSeed,
        LocusManifest[] Loci);

    public sealed record LocalizationDatasetManifest(
        int SchemaVersion,
        string DatasetId,
        int DatasetSeed,
        int WindowSizeBp,
        int K,
        int LociPerReference,
        string CoordinateConvention,
        int[] PrimaryQueryLengthsBp,
        int[] StressQueryLengthsBp,
        ReferenceManifest[] References);

    public static class PhaseOneDataset
    {
        private const int schemaVersion = 1;
        private const int windowSizeBp = 500;
        private const int lociPerReference = 20;
        private const string datasetId = "phase1-exact-v1";
        private const string coordinateConvention = "0-based-half-open";
        private static readonly int[] primaryQueryLengthsBp = { 2_000, 5_000, 10_000 };
        private static readonly int[] stressQueryLengthsBp = { 500, 1_000 };
        private static readonly Regex accessionPattern = new(@"(?:GCF|GCA)_\d+\.\d+", RegexOptions.Compiled | RegexOptions.CultureInvariant);
        private static readonly JsonSerializerOptions jsonOptions = new() { WriteIndented = true };
        private static readonly string[] referenceFileNames =
        {
            "acinetobacter_baumannii_atcc17978_GCF_013372085.1.fna.gz",
            "bacillus_anthracis_ames_GCF_000007845.1.fna.gz",
            "bacillus_cereus_atcc14579_GCF_006094295.1.fna.gz",
            "borrelia_burgdorferi_b31_GCF_000008685.2.fna.gz",
            "bsubtilis_168_GCF_000009045.1.fna.gz",
            "campylobacter_jejuni_nctc11168_GCF_000009085.1.fna.gz",
            "clostridioides_difficile_630_GCF_000009205.2.fna.gz",
            "deinococcus_radiodurans_r1_GCF_020546685.1.fna.gz",
            "ecoli_k12_mg1655_GCF_000005845.2.fna.gz",
            "ecoli_o157h7_sakai_GCF_000008865.2.fna.gz",
            "enterococcus_faecalis_v583_GCF_000007785.1.fna.gz",
            "haemophilus_influenzae_rd_kw20_GCF_000027305.1.fna.gz",
            "helicobacter_pylori_26695_GCF_000307795.1.fna.gz",
            "klebsiella_pneumoniae_mgh78578_GCF_000016305.1.fna.gz",
            "listeria_monocytogenes_egde_GCF_000196035.1.fna.gz",
            "mycobacterium_smegmatis_mc2155_GCF_000767605.1.fna.gz",
            "mycobacterium_tuberculosis_h37rv_GCF_000195955.2.fna.gz",
            "mycoplasma_genitalium_g37_GCF_000027325.1.fna.gz",
            "mycoplasma_pneumoniae_m129_GCF_910574535.1.fna.gz",
            "pseudomonas_aeruginosa_pao1_GCF_000006765.1.fna.gz",
            "salmonella_typhimurium_lt2_GCF_000006945.2.fna.gz",
            "shigella_flexneri_2a_301_GCF_000006925.2.fna.gz",
            "staphylococcus_aureus_nctc8325_GCF_000013425.1.fna.gz",
            "staphylococcus_epidermidis_atcc12228_GCF_000007645.1.fna.gz",
            "streptococcus_pneumoniae_tigr4_GCF_000273445.1.fna.gz",
            "streptococcus_pyogenes_m1gas_GCF_043231225.1.fna.gz",
            "streptomyces_coelicolor_a3_GCF_008931305.1.fna.gz",
            "vibrio_cholerae_n16961_GCF_900205735.1.fna.gz",
            "yersinia_pestis_co92_GCF_000009065.1.fna.gz",
            "yersinia_pseudotuberculosis_ip32953_GCF_000047365.1.fna.gz",
        };
        private static readonly IReadOnlyList<string> referenceFileNamesView = Array.AsReadOnly(referenceFileNames);

        public const int defaultSeed = 20261001;

        public static IReadOnlyList<string> ReferenceFileNames => referenceFileNamesView;

        public static LocalizationDatasetManifest Generate(
            string genomesDirectory,
            string outputPath,
            int seed = defaultSeed)
        {
            var selectedFiles = referenceFileNames.OrderBy(f => f, StringComparer.Ordinal).ToArray();

            var references = new ReferenceManifest[selectedFiles.Length];
            for (int i = 0; i < selectedFiles.Length; i++)
            {
                var fileName = selectedFiles[i];
                if (Path.GetFileName(fileName) != fileName || !fileName.EndsWith(".fna.gz", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidDataException($"Reference must be a .fna.gz file name, not a path: '{fileName}'.");

                var path = Path.Combine(genomesDirectory, fileName);
                if (!File.Exists(path))
                    throw new FileNotFoundException($"Reference genome not found: {path}", path);

                var record = FastaReader.ReadPrimaryRecord(path);
                var accessionMatch = accessionPattern.Match(fileName);
                if (!accessionMatch.Success)
                    throw new InvalidDataException($"Could not read a versioned assembly accession from '{fileName}'.");

                int samplingSeed = DeriveSamplingSeed(seed, fileName);
                var loci = CreateLoci(record, accessionMatch.Value, samplingSeed);
                using var source = File.OpenRead(path);
                string compressedSha256 = Convert.ToHexString(SHA256.HashData(source)).ToLowerInvariant();

                references[i] = new ReferenceManifest(
                    fileName,
                    accessionMatch.Value,
                    record.Id,
                    record.Length,
                    compressedSha256,
                    samplingSeed,
                    loci);
            }

            var manifest = new LocalizationDatasetManifest(
                schemaVersion,
                datasetId,
                seed,
                windowSizeBp,
                4,
                lociPerReference,
                coordinateConvention,
                primaryQueryLengthsBp,
                stressQueryLengthsBp,
                references);

            ValidateManifest(genomesDirectory, manifest);

            var fullOutputPath = Path.GetFullPath(outputPath);
            var outputDirectory = Path.GetDirectoryName(fullOutputPath)!;
            Directory.CreateDirectory(outputDirectory);
            using var output = new FileStream(fullOutputPath, FileMode.CreateNew, FileAccess.Write, FileShare.None);
            JsonSerializer.Serialize(output, manifest, jsonOptions);
            return manifest;
        }

        public static (int References, int Queries) Validate(string genomesDirectory, string manifestPath)
        {
            using var input = File.OpenRead(manifestPath);
            var manifest = JsonSerializer.Deserialize<LocalizationDatasetManifest>(input, jsonOptions)
                ?? throw new InvalidDataException("Manifest is empty or invalid JSON.");
            ValidateManifest(genomesDirectory, manifest);
            return (manifest.References.Length, manifest.References.Sum(r => r.Loci.Sum(l => l.Queries.Length)));
        }

        private static LocusManifest[] CreateLoci(FastaRecord record, string accession, int seed)
        {
            int maximumQueryLength = primaryQueryLengthsBp.Max();
            int availableStartSlots = (record.Length - maximumQueryLength) / windowSizeBp + 1;
            if (availableStartSlots < 1)
                throw new InvalidDataException($"Record '{record.Id}' is shorter than the maximum query length ({maximumQueryLength} bp).");

            var random = new StableRandom(seed);
            var starts = new List<int>(lociPerReference);
            int minimumSpacingSlots = maximumQueryLength / windowSizeBp;
            for (int i = 0; i < lociPerReference; i++)
            {
                int firstSlot = (int)((long)i * availableStartSlots / lociPerReference);
                int afterLastSlot = (int)((long)(i + 1) * availableStartSlots / lociPerReference);
                int firstCandidate = firstSlot + minimumSpacingSlots / 2;
                int afterLastCandidate = afterLastSlot - minimumSpacingSlots / 2;
                int slotsInStratum = afterLastCandidate - firstCandidate;
                if (slotsInStratum <= minimumSpacingSlots)
                    throw new InvalidDataException($"Record '{record.Id}' is too short to sample {lociPerReference} separated loci.");

                starts.Add((firstCandidate + random.NextInt(slotsInStratum)) * windowSizeBp);
            }

            starts.Sort();
            var allQueryLengths = primaryQueryLengthsBp.Select(length => (Length: length, Stratum: "Primary"))
                .Concat(stressQueryLengthsBp.Select(length => (Length: length, Stratum: "Stress")))
                .ToArray();

            return starts.Select((start, index) => new LocusManifest(
                $"{accession}:locus-{index + 1:D2}",
                start,
                allQueryLengths.Select(query =>
                {
                    string sequence = record.Sequence.Substring(start, query.Length);
                    string digest = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(sequence))).ToLowerInvariant();
                    return new QueryManifest(
                        $"{accession}:locus-{index + 1:D2}:{query.Length}",
                        query.Stratum,
                        query.Length,
                        digest);
                }).ToArray())).ToArray();
        }

        private static void ValidateManifest(string genomesDirectory, LocalizationDatasetManifest manifest)
        {
            if (manifest.SchemaVersion != schemaVersion || manifest.DatasetId != datasetId ||
                manifest.WindowSizeBp != windowSizeBp || manifest.K != 4 ||
                manifest.LociPerReference != lociPerReference || manifest.CoordinateConvention != coordinateConvention)
                throw new InvalidDataException("Manifest protocol metadata does not match phase1-exact-v1.");
            if (!manifest.PrimaryQueryLengthsBp.SequenceEqual(primaryQueryLengthsBp) ||
                !manifest.StressQueryLengthsBp.SequenceEqual(stressQueryLengthsBp))
                throw new InvalidDataException("Manifest query lengths do not match the locked Phase 1 protocol.");
            var expectedFiles = referenceFileNames.OrderBy(f => f, StringComparer.Ordinal);
            if (!manifest.References.Select(r => r.FileName).SequenceEqual(expectedFiles, StringComparer.Ordinal))
                throw new InvalidDataException("Manifest references do not match the curated 30-genome bacterial corpus.");

            foreach (var reference in manifest.References)
            {
                if (Path.GetFileName(reference.FileName) != reference.FileName ||
                    !reference.FileName.EndsWith(".fna.gz", StringComparison.OrdinalIgnoreCase))
                    throw new InvalidDataException($"Unsafe or invalid reference file name '{reference.FileName}'.");

                var path = Path.Combine(genomesDirectory, reference.FileName);
                if (!File.Exists(path))
                    throw new FileNotFoundException($"Reference genome not found: {path}", path);

                using var source = File.OpenRead(path);
                string compressedSha256 = Convert.ToHexString(SHA256.HashData(source)).ToLowerInvariant();
                if (!string.Equals(compressedSha256, reference.CompressedSha256, StringComparison.Ordinal))
                    throw new InvalidDataException($"Compressed FASTA checksum mismatch for '{reference.FileName}'.");

                var record = FastaReader.ReadPrimaryRecord(path);
                if (record.Id != reference.RecordId || record.Length != reference.RecordLengthBp)
                    throw new InvalidDataException($"Primary FASTA record changed for '{reference.FileName}'.");
                if (reference.Loci.Length != lociPerReference)
                    throw new InvalidDataException($"'{reference.FileName}' must have exactly {lociPerReference} loci.");

                var accessionMatch = accessionPattern.Match(reference.FileName);
                if (!accessionMatch.Success || reference.Accession != accessionMatch.Value)
                    throw new InvalidDataException($"Versioned accession does not match '{reference.FileName}'.");
                int expectedSamplingSeed = DeriveSamplingSeed(manifest.DatasetSeed, reference.FileName);
                if (reference.SamplingSeed != expectedSamplingSeed)
                    throw new InvalidDataException($"Sampling seed mismatch for '{reference.FileName}'.");
                var expectedLoci = CreateLoci(record, reference.Accession, expectedSamplingSeed);

                int previousEnd = -1;
                for (int locusIndex = 0; locusIndex < reference.Loci.Length; locusIndex++)
                {
                    var locus = reference.Loci[locusIndex];
                    if (locus.StartBp < 0 || locus.StartBp % windowSizeBp != 0 ||
                        locus.StartBp < previousEnd || locus.StartBp + primaryQueryLengthsBp.Max() > record.Length)
                        throw new InvalidDataException($"Invalid or overlapping locus '{locus.LocusId}'.");
                    if (locus.LocusId != $"{reference.Accession}:locus-{locusIndex + 1:D2}")
                        throw new InvalidDataException($"Unexpected locus ID '{locus.LocusId}'.");
                    if (locus.StartBp != expectedLoci[locusIndex].StartBp)
                        throw new InvalidDataException($"Sampled coordinate mismatch for '{locus.LocusId}'.");

                    var expectedQueries = primaryQueryLengthsBp.Select(length => (Length: length, Stratum: "Primary"))
                        .Concat(stressQueryLengthsBp.Select(length => (Length: length, Stratum: "Stress")))
                        .ToArray();
                    if (locus.Queries.Length != expectedQueries.Length)
                        throw new InvalidDataException($"Locus '{locus.LocusId}' has an unexpected query count.");

                    for (int queryIndex = 0; queryIndex < expectedQueries.Length; queryIndex++)
                    {
                        var query = locus.Queries[queryIndex];
                        var expectedQuery = expectedLoci[locusIndex].Queries[queryIndex];
                        var expected = expectedQueries[queryIndex];
                        if (query.LengthBp != expected.Length || query.Stratum != expected.Stratum ||
                            query.QueryId != $"{locus.LocusId}:{query.LengthBp}" ||
                            query.SequenceSha256 != expectedQuery.SequenceSha256)
                            throw new InvalidDataException($"Unexpected query recipe at locus '{locus.LocusId}'.");

                        string sequence = record.Sequence.Substring(locus.StartBp, query.LengthBp);
                        string digest = Convert.ToHexString(SHA256.HashData(Encoding.UTF8.GetBytes(sequence))).ToLowerInvariant();
                        if (!string.Equals(digest, query.SequenceSha256, StringComparison.Ordinal))
                            throw new InvalidDataException($"Query checksum mismatch for '{query.QueryId}'.");
                    }

                    previousEnd = locus.StartBp + primaryQueryLengthsBp.Max();
                }
            }
        }

        private static int DeriveSamplingSeed(int seed, string fileName)
        {
            unchecked
            {
                uint hash = 2166136261;
                foreach (byte value in Encoding.UTF8.GetBytes($"{seed}|{fileName}"))
                {
                    hash ^= value;
                    hash *= 16777619;
                }
                return (int)hash;
            }
        }

        private struct StableRandom
        {
            private ulong state;

            public StableRandom(int seed) => state = unchecked((uint)seed);

            public int NextInt(int exclusiveMaximum)
            {
                ulong bound = (uint)exclusiveMaximum;
                ulong threshold = unchecked(0UL - bound) % bound;
                ulong value;
                do
                {
                    value = NextUInt64();
                } while (value < threshold);
                return (int)(value % bound);
            }

            private ulong NextUInt64()
            {
                state += 0x9E3779B97F4A7C15UL;
                ulong value = state;
                value = (value ^ (value >> 30)) * 0xBF58476D1CE4E5B9UL;
                value = (value ^ (value >> 27)) * 0x94D049BB133111EBUL;
                return value ^ (value >> 31);
            }
        }
    }
}