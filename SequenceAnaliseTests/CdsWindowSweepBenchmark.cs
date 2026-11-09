using System.Globalization;
using System.Text.RegularExpressions;
using SequenceAnalyses;

namespace SequenceAnaliseTests
{
    using P = ProteinBackTranslationFeasibility;

    /// <summary>
    /// CDS window sweep (2026-11-04, your ideas): a genome as the Z-curve profiles of the segments that code for protein, with smaller
    /// windows inside each CDS, the position of each base inside its window kept through Legendre moments, and the z (GC) axis optionally
    /// left out. Per window and kept axis: means keeps the mean; moments keeps P0, P1 and P2 over the base's position t in [-1, 1], each
    /// scaled by sqrt(2k + 1) so all three carry comparable noise. Distance between two genes: mean over windows of the Euclidean distance.
    /// Search (bacteria, R6f's pool, picks and positions-1-and-2 Prefix profiles at 1,500 bp; windows of 40 / 80 / 166 / 333 / 1,000
    /// positions-1-and-2 bases, about 60 / 120 / 250 / 500 / 1,500 bp of gene; 333 xyz means is R6f). Queries Dna, ExpectedUniform and
    /// ExpectedPerTarget from the protein. Measured: top-1 among the 30-genome pool and the pool-size curve. Extrapolated to the catalogue
    /// with R6f's model, corrected after the first run: a query some competitor already beats keeps the measured share closer, and only
    /// a query no competitor beats is extrapolated from the tail (P.CloserShare's form below the nearest competitor). A cell's extrapolation
    /// counts only if the model fitted on 8-genome pools predicts the measured 30-genome top-1 within 5 points.
    /// Homologs: R6f's 200 E. coli K-12 queries against O157:H7, Shigella, Salmonella and Klebsiella, per cell.
    /// Trees (bacteria and viruses, GeneLevelSignatureBenchmark's genes, first 999 bp of genes at least that long, all positions): windows
    /// of 60 / 120 / 250 / 500 / 999 bp; genome distance = mean best CDS-to-CDS distance, both ways. Karlin and cvTree4 as context.
    /// Rule, fixed before the runs: search passes if a valid cell's ExpectedUniform or ExpectedPerTarget reaches 50% at catalogue scale;
    /// trees pass if a cell beats cvTree4 on bacterial phylum AUC and genus clades (NJ + UPGMA). If both fail, no per-CDS Z-curve row.
    /// Writes CdsSweep_Search.csv, CdsSweep_Homolog.csv and CdsSweep_Trees.csv.
    /// </summary>
    public class CdsWindowSweepBenchmark
    {
        private const int searchLength = 1500;
        private const int treeLength = 999;
        private const int searchPos12Length = searchLength / 3 * 2;
        private const int catalogueGenomes = 15_314;
        private const int tailNeighbours = 100;
        private const int draws = 20;
        private const int fitGenomes = 8;
        private const double modelTolerance = 0.05;
        private static readonly int[] searchWindows = { 40, 80, 166, 333, 1000 };
        private static readonly int[] treeWindows = { 60, 120, 250, 500, 999 };
        private static readonly int[] poolSizes = { 1, 2, 4, 8, 16, 30 };
        private static readonly (string Name, int[] Axes)[] axisSets = { ("xyz", new[] { 0, 1, 2 }), ("xy", new[] { 0, 1 }) };

        [Test]
        public void ProteinSearchAcrossCdsWindows()
        {
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();
            var fileNames = P.Corpus().Where(c => c.Group == "Bacteria").Select(c => c.FileName).ToList();
            int genomeCount = fileNames.Count;
            var pool = new List<(int Genome, string Protein, string Pos12)>();
            var ownTables = new List<(double X, double Y, double Z)[]>();
            var uniformTables = new List<(double X, double Y, double Z)[]>();
            var queries = new List<int>();
            for (int genome = 0; genome < genomeCount; genome++)
            {
                var fileName = fileNames[genome];
                var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;
                var reverse = P.ReverseComplement(forward);
                var code = Path.GetFileName(fileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? P.mycoplasmaCode : P.standardCode;
                var synonyms = P.Synonyms(code);
                ownTables.Add(P.ExpectedAxesTable(synonyms, P.CodonUsage(forward, P.OpenStretches(forward, code, P.codonUsageMinLength))));
                uniformTables.Add(P.ExpectedAxesTable(synonyms, uniformUsage));
                var genes = P.AtgGenes(forward, code, searchLength).Select(g => forward.Substring(g.Start, g.End - g.Start))
                    .Concat(P.AtgGenes(reverse, code, searchLength).Select(g => reverse.Substring(g.Start, g.End - g.Start)))
                    .ToList();
                int first = pool.Count;
                foreach (var gene in genes)
                    pool.Add((genome, P.Translate(gene[..searchLength], code), P.CodonPositions12(gene[..searchLength])));
                // R6d's picker for this length, so the queries are R6e's and R6f's.
                var picker = new Random(P.StableSeed(fileName) ^ searchLength ^ 0x2545f491);
                queries.AddRange(Enumerable.Range(0, genes.Count).OrderBy(_ => picker.Next()).Take(P.regionsPerGenome).OrderBy(i => i).Select(i => first + i));
            }
            var homologCases = HomologCases(uniformUsage);
            var genesOfGenome = new int[genomeCount];
            foreach (var gene in pool)
                genesOfGenome[gene.Genome]++;

            var rows = new List<string>();
            var homologRows = new List<string>();
            foreach (var (axisName, axes) in axisSets)
                foreach (var window in searchWindows)
                    foreach (var moments in new[] { false, true })
                    {
                        var cell = $"{axisName},{window},{window * 3 / 2},{(moments ? "moments" : "means")}";
                        var targets = pool.Select(g => Features(DnaAxes(g.Pos12), searchPos12Length, window, moments, axes)).ToArray();
                        foreach (var variant in new[] { "Dna", "ExpectedUniform", "ExpectedPerTarget" })
                        {
                            var results = new (double Observed, double From8, double Catalogue, double[] Curve)[queries.Count];
                            Parallel.For(0, queries.Count, j =>
                            {
                                int truth = queries[j];
                                int own = pool[truth].Genome;
                                var queryByGenome = new double[genomeCount][][];
                                double[][] Query(int genome) => queryByGenome[variant == "ExpectedPerTarget" ? genome : 0] ??= variant switch
                                {
                                    "Dna" => targets[truth],
                                    "ExpectedUniform" => Features(ExpectedAxes(pool[truth].Protein, uniformTables[own]), searchPos12Length, window, moments, axes),
                                    _ => Features(ExpectedAxes(pool[truth].Protein, ownTables[genome]), searchPos12Length, window, moments, axes),
                                };
                                var trueDistance = Distance(Query(own), targets[truth]);
                                int closerOwn = 0;
                                var closerByGenome = new int[genomeCount];
                                var others = new List<(double Distance, int Genome)>(pool.Count);
                                for (int i = 0; i < pool.Count; i++)
                                {
                                    if (i == truth) continue;
                                    var distance = Distance(Query(pool[i].Genome), targets[i]);
                                    if (pool[i].Genome == own)
                                    {
                                        if (distance < trueDistance) closerOwn++;
                                    }
                                    else
                                    {
                                        others.Add((distance, pool[i].Genome));
                                        if (distance < trueDistance) closerByGenome[pool[i].Genome]++;
                                    }
                                }

                                var rng = new Random(P.StableSeed(pool[truth].Protein) ^ window);
                                var otherGenomes = Enumerable.Range(0, genomeCount).Where(g => g != own).ToArray();
                                var curve = new double[poolSizes.Length];
                                for (int s = 0; s < poolSizes.Length; s++)
                                    for (int d = 0; d < draws; d++)
                                    {
                                        var chosen = otherGenomes.OrderBy(_ => rng.Next()).Take(poolSizes[s] - 1);
                                        if (closerOwn == 0 && chosen.All(g => closerByGenome[g] == 0)) curve[s] += 1.0 / draws;
                                    }
                                var all = others.Select(o => o.Distance).Order().ToArray();
                                var fitSet = otherGenomes.OrderBy(_ => rng.Next()).Take(fitGenomes - 1).ToHashSet();
                                var fit = others.Where(o => fitSet.Contains(o.Genome)).Select(o => o.Distance).Order().ToArray();
                                results[j] = (
                                    closerOwn == 0 && closerByGenome.Sum() == 0 ? 1 : 0,
                                    Survive(closerOwn, fit, trueDistance, all.Length),
                                    Survive(closerOwn, all, trueDistance, catalogueGenomes / (double)genomeCount * pool.Count),
                                    curve);
                            });
                            double observed = results.Average(r => r.Observed), from8 = results.Average(r => r.From8);
                            var valid = Math.Abs(from8 - observed) <= modelTolerance;
                            var curveText = string.Join(",", Enumerable.Range(0, poolSizes.Length).Select(s => Format(results.Average(r => r.Curve[s]))));
                            var line = $"{cell},{variant},{queries.Count},{pool.Count},{Format(observed)},{Format(from8)},{valid.ToString().ToLowerInvariant()},{curveText},{Format(results.Average(r => r.Catalogue))}";
                            rows.Add(line);
                            TestContext.Out.WriteLine(line);
                        }

                        foreach (var line in HomologRanks(homologCases, cell, window, moments, axes))
                            homologRows.Add(line);
                    }

            using (var csv = File.CreateText("CdsSweep_Search.csv"))
            {
                csv.WriteLine("Axes,Pos12Window,ApproxGeneBp,Features,Variant,Queries,Pool,Top1At30,Model8To30,ModelValid," +
                    string.Join(",", poolSizes.Select(s => $"Top1At{s}Genomes")) + ",ModelTop1AtCatalogue");
                foreach (var row in rows)
                    csv.WriteLine(row);
            }
            using (var csv = File.CreateText("CdsSweep_Homolog.csv"))
            {
                csv.WriteLine("Axes,Pos12Window,ApproxGeneBp,Features,Target,Variant,Queries,Candidates,Top1,Top10");
                foreach (var row in homologRows)
                    csv.WriteLine(row);
                foreach (var row in homologRows)
                    TestContext.Out.WriteLine(row);
            }
            var passing = rows.Select(r => r.Split(',')).Where(r => r[4] != "Dna" && r[9] == "true")
                .Select(r => (Cell: string.Join("/", r[0], r[1], r[3], r[4]), Value: double.Parse(r[^1], CultureInfo.InvariantCulture)))
                .OrderByDescending(r => r.Value).FirstOrDefault();
            TestContext.Out.WriteLine($"Search rule: best valid query-time cell {passing.Cell} at catalogue scale {Format(passing.Value)} -> {(passing.Value >= 0.5 ? "pass" : "stop")}");
            Assert.That(rows.Select(r => r.Split(',')).Where(r => r[4] == "Dna").All(r => double.Parse(r[7], CultureInfo.InvariantCulture) >= 0.9));
        }

        /// <summary>
        /// Follow-ups to the sweep (2026-11-04), on the chosen cell: 500 bp-of-gene windows (333 positions-1-and-2 bases) over x, y, z.
        /// Feature sets with the Euclidean distance: means (3 per window), trend (mean and P1, 6), moments (P0-P2, 9). Distance variants on
        /// moments: capped (each window's distance capped at the median per-window distance between random pool genes), whitened (features
        /// multiplied by the inverse square root of their covariance over all pool windows, computed once), uncertainty (each squared
        /// difference over the query feature's variance under the codon choices plus epsilon = 0.1 / window, the iid-null variance's
        /// tenth). Same queries, pool, extrapolation and model check as ProteinSearchAcrossCdsWindows, and the homolog check per option.
        /// Rule fixed before the run: the simplest option, in the order listed, within 2 points of the best on ExpectedUniform's
        /// 30-genome top-1 and on its mean ortholog top-1 over the four targets, with a valid extrapolation. Writes CdsFollowUp_Search.csv.
        /// </summary>
        [Test]
        public void ProteinSearchFollowUps()
        {
            const int window = 333;
            var axes = new[] { 0, 1, 2 };
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();
            var setup = LoadSearchSetup(uniformUsage);
            var homologCases = HomologCases(uniformUsage);
            var standardUniform = P.ExpectedAxesTable(P.Synonyms(P.standardCode), uniformUsage);
            var options = new (string Name, int Orders, string Distance)[]
            {
                ("means", 1, "euclid"), ("trend", 2, "euclid"), ("moments", 3, "euclid"),
                ("moments", 3, "capped"), ("moments", 3, "whitened"), ("moments", 3, "uncertainty"),
            };

            var rows = new List<(string Option, string Variant, double Observed, double From8, bool Valid, double Catalogue, double[] Orthologs)>();
            foreach (var (name, orders, distanceName) in options)
            {
                double[][] Transformed(double[][] features, double[,]? matrix) => matrix is null ? features : features.Select(f => Multiply(matrix, f)).ToArray();
                var raw = setup.Pool.Select(g => FeaturesOfOrder(DnaAxes(g.Pos12), searchPos12Length, window, orders, axes)).ToArray();
                var whitening = distanceName == "whitened" ? Whitening(raw.SelectMany(f => f).ToArray()) : null;
                var targets = raw.Select(f => Transformed(f, whitening)).ToArray();
                double cap = double.PositiveInfinity;
                if (distanceName == "capped")
                {
                    var rng = new Random(20261104);
                    var sample = Enumerable.Range(0, 20_000).Select(_ => targets[rng.Next(targets.Length)].Zip(targets[rng.Next(targets.Length)], WindowDistance)).SelectMany(d => d).Order().ToArray();
                    cap = sample[sample.Length / 2];
                }
                double Measure(SearchQuery query, double[][] target) => distanceName switch
                {
                    "capped" => query.Features.Zip(target, (q, t) => Math.Min(WindowDistance(q, t), cap)).Average(),
                    "uncertainty" => query.Features.Select((q, w) => WeightedWindowDistance(q, target[w], query.Variances?[w], 0.1 / window)).Average(),
                    _ => Distance(query.Features, target),
                };
                SearchQuery Build(string variant, string protein, string? pos12, (double X, double Y, double Z)[]? table)
                {
                    if (variant == "Dna") return new SearchQuery(Transformed(FeaturesOfOrder(DnaAxes(pos12!), searchPos12Length, window, orders, axes), whitening), null);
                    var features = Transformed(FeaturesOfOrder(ExpectedAxes(protein, table!), searchPos12Length, window, orders, axes), whitening);
                    return new SearchQuery(features, FeatureVariances(ExpectedAxes(protein, table!), searchPos12Length, window, orders, axes));
                }

                foreach (var variant in new[] { "Dna", "ExpectedUniform", "ExpectedPerTarget" })
                {
                    var results = new (double Observed, double From8, double Catalogue)[setup.Queries.Count];
                    Parallel.For(0, setup.Queries.Count, j =>
                    {
                        int truth = setup.Queries[j];
                        int own = setup.Pool[truth].Genome;
                        var byGenome = new SearchQuery[setup.GenomeCount];
                        SearchQuery Query(int genome) => byGenome[variant == "ExpectedPerTarget" ? genome : 0] ??= Build(variant, setup.Pool[truth].Protein,
                            setup.Pool[truth].Pos12, variant == "ExpectedPerTarget" ? setup.OwnTables[genome] : setup.UniformTables[own]);
                        var trueDistance = Measure(Query(own), targets[truth]);
                        int closerOwn = 0, closerOther = 0;
                        var others = new List<(double Distance, int Genome)>(setup.Pool.Count);
                        for (int i = 0; i < setup.Pool.Count; i++)
                        {
                            if (i == truth) continue;
                            var distance = Measure(Query(setup.Pool[i].Genome), targets[i]);
                            if (setup.Pool[i].Genome == own) { if (distance < trueDistance) closerOwn++; }
                            else
                            {
                                others.Add((distance, setup.Pool[i].Genome));
                                if (distance < trueDistance) closerOther++;
                            }
                        }
                        var rng = new Random(P.StableSeed(setup.Pool[truth].Protein) ^ window);
                        var fitSet = Enumerable.Range(0, setup.GenomeCount).Where(g => g != own).OrderBy(_ => rng.Next()).Take(fitGenomes - 1).ToHashSet();
                        var all = others.Select(o => o.Distance).Order().ToArray();
                        var fit = others.Where(o => fitSet.Contains(o.Genome)).Select(o => o.Distance).Order().ToArray();
                        results[j] = (closerOwn == 0 && closerOther == 0 ? 1 : 0,
                            Survive(closerOwn, fit, trueDistance, all.Length),
                            Survive(closerOwn, all, trueDistance, catalogueGenomes / (double)setup.GenomeCount * setup.Pool.Count));
                    });

                    var orthologs = homologCases.GroupBy(c => c.Target).Select(group =>
                    {
                        var candidates = group.First().CandidatePos12.Select(c => Transformed(FeaturesOfOrder(DnaAxes(c), searchPos12Length, window, orders, axes), whitening)).ToArray();
                        return group.AsParallel().Average(c =>
                        {
                            var query = Build(variant, c.Protein, c.QueryPos12, variant == "ExpectedPerTarget" ? c.TargetTable : standardUniform);
                            var trueDistance = Measure(query, candidates[c.TrueSlot]);
                            return Enumerable.Range(0, candidates.Length).Any(s => s != c.TrueSlot && Measure(query, candidates[s]) < trueDistance) ? 0.0 : 1.0;
                        });
                    }).ToArray();
                    double observed = results.Average(r => r.Observed), from8 = results.Average(r => r.From8);
                    rows.Add(($"{name}_{distanceName}", variant, observed, from8, Math.Abs(from8 - observed) <= modelTolerance, results.Average(r => r.Catalogue), orthologs));
                    TestContext.Out.WriteLine($"{name}_{distanceName},{variant},{Format(observed)},{Format(from8)},{Format(results.Average(r => r.Catalogue))},{string.Join("/", orthologs.Select(Format))}");
                }
            }

            using (var csv = File.CreateText("CdsFollowUp_Search.csv"))
            {
                csv.WriteLine("Option,NumbersPerWindow,Variant,Top1At30,Model8To30,ModelValid,ModelTop1AtCatalogue," + string.Join(",", P.homologTargets.Select(t => "Ortholog_" + t.Split('_')[0] + "_" + t.Split('_')[1])));
                foreach (var row in rows)
                {
                    var numbers = 3 * options.First(o => $"{o.Name}_{o.Distance}" == row.Option).Orders;
                    csv.WriteLine($"{row.Option},{numbers},{row.Variant},{Format(row.Observed)},{Format(row.From8)},{row.Valid.ToString().ToLowerInvariant()},{Format(row.Catalogue)},{string.Join(",", row.Orthologs.Select(Format))}");
                }
            }
            var uniform = rows.Where(r => r.Variant == "ExpectedUniform").ToList();
            double bestTop1 = uniform.Max(r => r.Observed), bestOrtholog = uniform.Max(r => r.Orthologs.Average());
            var chosen = uniform.FirstOrDefault(r => r.Valid && r.Observed >= bestTop1 - 0.02 && r.Orthologs.Average() >= bestOrtholog - 0.02);
            TestContext.Out.WriteLine($"Rule (best top-1 {Format(bestTop1)}, best mean ortholog {Format(bestOrtholog)}): choose {chosen.Option ?? "none"}");
            Assert.That(rows.Where(r => r.Variant == "Dna").All(r => r.Observed >= 0.9));
        }

        /// <summary>
        /// Storage precision (2026-11-04): the chosen gene profile (500 bp-of-gene windows, moments, x y z, Euclidean) with the stored
        /// target features as 32-bit floats or rounded to 16-bit halves; queries are computed at query time in double. Same pool, queries,
        /// model and homolog check as ProteinSearchFollowUps. Rule fixed before the run: store halves if every ExpectedUniform and
        /// ExpectedPerTarget number (30-genome top-1, catalogue model, each ortholog top-1) is within 0.5 points of float32.
        /// Writes CdsPrecision_Search.csv.
        /// </summary>
        [Test]
        public void HalfPrecisionTargets()
        {
            const int window = 333, orders = 3;
            var axes = new[] { 0, 1, 2 };
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();
            var setup = LoadSearchSetup(uniformUsage);
            var homologCases = HomologCases(uniformUsage);
            var standardUniform = P.ExpectedAxesTable(P.Synonyms(P.standardCode), uniformUsage);
            double[][][] Stored(double[][][] genes, string precision, IReadOnlyList<int> recordOf) => precision switch
            {
                "half" => HalfRounded(genes),
                "byte" => ByteRounded(genes, recordOf),
                _ => genes.Select(g => g.Select(w => w.Select(v => (double)(float)v).ToArray()).ToArray()).ToArray(),
            };

            var lines = new List<string>();
            var numbers = new Dictionary<(string Precision, string Variant), double[]>();
            // byte (added 2026-11-05): 8 bits per value against each record's largest |value| of that feature (ByteRounded).
            foreach (var precision in new[] { "float", "half", "byte" })
            {
                var targets = Stored(setup.Pool.Select(g => FeaturesOfOrder(DnaAxes(g.Pos12), searchPos12Length, window, orders, axes)).ToArray(),
                    precision, setup.Pool.Select(g => g.Genome).ToList());
                foreach (var variant in new[] { "Dna", "ExpectedUniform", "ExpectedPerTarget" })
                {
                    double[][] Build(string protein, string pos12, (double X, double Y, double Z)[] table) => variant == "Dna"
                        ? FeaturesOfOrder(DnaAxes(pos12), searchPos12Length, window, orders, axes)
                        : FeaturesOfOrder(ExpectedAxes(protein, table), searchPos12Length, window, orders, axes);
                    var results = new (double Observed, double Catalogue)[setup.Queries.Count];
                    Parallel.For(0, setup.Queries.Count, j =>
                    {
                        int truth = setup.Queries[j];
                        int own = setup.Pool[truth].Genome;
                        var byGenome = new double[setup.GenomeCount][][];
                        double[][] Query(int genome) => byGenome[variant == "ExpectedPerTarget" ? genome : 0] ??= Build(setup.Pool[truth].Protein,
                            setup.Pool[truth].Pos12, variant == "ExpectedPerTarget" ? setup.OwnTables[genome] : setup.UniformTables[own]);
                        var trueDistance = Distance(Query(own), targets[truth]);
                        int closerOwn = 0, closerOther = 0;
                        var others = new List<double>(setup.Pool.Count);
                        for (int i = 0; i < setup.Pool.Count; i++)
                        {
                            if (i == truth) continue;
                            var distance = Distance(Query(setup.Pool[i].Genome), targets[i]);
                            if (setup.Pool[i].Genome == own) { if (distance < trueDistance) closerOwn++; }
                            else
                            {
                                others.Add(distance);
                                if (distance < trueDistance) closerOther++;
                            }
                        }
                        results[j] = (closerOwn == 0 && closerOther == 0 ? 1 : 0,
                            Survive(closerOwn, others.Order().ToArray(), trueDistance, catalogueGenomes / (double)setup.GenomeCount * setup.Pool.Count));
                    });
                    var orthologs = homologCases.GroupBy(c => c.Target).Select(group =>
                    {
                        var candidatePos12 = group.First().CandidatePos12;
                        var candidates = Stored(candidatePos12.Select(c => FeaturesOfOrder(DnaAxes(c), searchPos12Length, window, orders, axes)).ToArray(),
                            precision, new int[candidatePos12.Length]);
                        return group.AsParallel().Average(c =>
                        {
                            var query = Build(c.Protein, c.QueryPos12, variant == "ExpectedPerTarget" ? c.TargetTable : standardUniform);
                            var trueDistance = Distance(query, candidates[c.TrueSlot]);
                            return Enumerable.Range(0, candidates.Length).Any(s => s != c.TrueSlot && Distance(query, candidates[s]) < trueDistance) ? 0.0 : 1.0;
                        });
                    }).ToArray();
                    var values = new[] { results.Average(r => r.Observed), results.Average(r => r.Catalogue) }.Concat(orthologs).ToArray();
                    numbers[(precision, variant)] = values;
                    var line = $"{precision},{variant},{string.Join(",", values.Select(Format))}";
                    lines.Add(line);
                    TestContext.Out.WriteLine(line);
                }
            }

            using (var csv = File.CreateText("CdsPrecision_Search.csv"))
            {
                csv.WriteLine("Precision,Variant,Top1At30,ModelTop1AtCatalogue," + string.Join(",", P.homologTargets.Select(t => "Ortholog_" + t.Split('_')[0] + "_" + t.Split('_')[1])));
                foreach (var line in lines)
                    csv.WriteLine(line);
            }
            var worst = new[] { "ExpectedUniform", "ExpectedPerTarget" }
                .Max(v => numbers[("float", v)].Zip(numbers[("half", v)], (a, b) => Math.Abs(a - b)).Max());
            TestContext.Out.WriteLine($"Rule: largest difference {Format(worst)} -> {(worst <= 0.005 ? "store halves" : "store floats")}");
            var worstByte = new[] { "ExpectedUniform", "ExpectedPerTarget" }
                .Max(v => numbers[("half", v)].Zip(numbers[("byte", v)], (a, b) => Math.Abs(a - b)).Max());
            TestContext.Out.WriteLine($"Byte rule: largest difference from half {Format(worstByte)} -> {(worstByte <= 0.005 ? "8 bits pass here" : "keep halves")}");
        }

        private static double[][][] HalfRounded(double[][][] genes) =>
            genes.Select(g => g.Select(w => w.Select(v => (double)(Half)v).ToArray()).ToArray()).ToArray();

        // 8-bit storage (2026-11-05): each feature as a signed byte times that record's largest |value| of the feature / 127, the scale
        // kept with the row as a half; what's compared is the value read back.
        private static double[][][] ByteRounded(double[][][] genes, IReadOnlyList<int> recordOf)
        {
            int width = genes.First(g => g.Length > 0)[0].Length;
            var scales = new Dictionary<int, double[]>();
            for (int g = 0; g < genes.Length; g++)
            {
                if (!scales.TryGetValue(recordOf[g], out var scale)) scales[recordOf[g]] = scale = new double[width];
                foreach (var window in genes[g])
                    for (int f = 0; f < width; f++) scale[f] = Math.Max(scale[f], Math.Abs(window[f]));
            }
            foreach (var scale in scales.Values)
                for (int f = 0; f < width; f++) scale[f] = (double)(Half)scale[f];
            return genes.Select((gene, g) => gene.Select(window => window.Select((v, f) =>
            {
                double scale = scales[recordOf[g]][f];
                return scale == 0 ? 0 : Math.Clamp(Math.Round(v / scale * 127), -127, 127) * scale / 127;
            }).ToArray()).ToArray()).ToArray();
        }

        private sealed record SearchQuery(double[][] Features, double[][]? Variances);

        private sealed record SearchSetup(int GenomeCount, List<(int Genome, string Protein, string Pos12)> Pool,
            List<(double X, double Y, double Z)[]> OwnTables, List<(double X, double Y, double Z)[]> UniformTables, List<int> Queries);

        // R6f's pool (genes >= 1,500 bp from their first ATG, both strands, 30 bacteria) and R6d's 20 picks per genome.
        private static SearchSetup LoadSearchSetup(long[] uniformUsage)
        {
            var fileNames = P.Corpus().Where(c => c.Group == "Bacteria").Select(c => c.FileName).ToList();
            var pool = new List<(int Genome, string Protein, string Pos12)>();
            var ownTables = new List<(double X, double Y, double Z)[]>();
            var uniformTables = new List<(double X, double Y, double Z)[]>();
            var queries = new List<int>();
            for (int genome = 0; genome < fileNames.Count; genome++)
            {
                var fileName = fileNames[genome];
                var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileName)).Sequence;
                var reverse = P.ReverseComplement(forward);
                var code = Path.GetFileName(fileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? P.mycoplasmaCode : P.standardCode;
                var synonyms = P.Synonyms(code);
                ownTables.Add(P.ExpectedAxesTable(synonyms, P.CodonUsage(forward, P.OpenStretches(forward, code, P.codonUsageMinLength))));
                uniformTables.Add(P.ExpectedAxesTable(synonyms, uniformUsage));
                var genes = P.AtgGenes(forward, code, searchLength).Select(g => forward.Substring(g.Start, g.End - g.Start))
                    .Concat(P.AtgGenes(reverse, code, searchLength).Select(g => reverse.Substring(g.Start, g.End - g.Start)))
                    .ToList();
                int first = pool.Count;
                foreach (var gene in genes)
                    pool.Add((genome, P.Translate(gene[..searchLength], code), P.CodonPositions12(gene[..searchLength])));
                var picker = new Random(P.StableSeed(fileName) ^ searchLength ^ 0x2545f491);
                queries.AddRange(Enumerable.Range(0, genes.Count).OrderBy(_ => picker.Next()).Take(P.regionsPerGenome).OrderBy(i => i).Select(i => first + i));
            }
            return new SearchSetup(fileNames.Count, pool, ownTables, uniformTables, queries);
        }

        // Variance of each feature under the codon choices: a base's axis value is +-1, so its variance is 1 - mean^2, and a feature
        // sqrt(2k + 1) / W * sum(a_i P_k(t_i)) of independent bases has variance (2k + 1) / W^2 * sum(var_i P_k(t_i)^2).
        private static double[][] FeatureVariances(Func<int, (double X, double Y, double Z)> expectedAt, int length, int window, int orders, int[] axes)
        {
            var result = new double[length / window][];
            for (int w = 0; w < result.Length; w++)
            {
                var variances = new double[axes.Length * orders];
                for (int i = 0; i < window; i++)
                {
                    var (x, y, z) = expectedAt(w * window + i);
                    double t = -1 + 2 * (i + 0.5) / window;
                    for (int k = 0; k < orders; k++)
                    {
                        double weight = k switch { 0 => 1, 1 => t, _ => (3 * t * t - 1) / 2 };
                        for (int a = 0; a < axes.Length; a++)
                        {
                            double mean = axes[a] switch { 0 => x, 1 => y, _ => z };
                            variances[axes.Length * k + a] += (1 - mean * mean) * weight * weight;
                        }
                    }
                }
                for (int k = 0; k < orders; k++)
                    for (int a = 0; a < axes.Length; a++)
                        variances[axes.Length * k + a] *= (2 * k + 1) / ((double)window * window);
                result[w] = variances;
            }
            return result;
        }

        private static double WindowDistance(double[] a, double[] b)
        {
            double squared = 0;
            for (int f = 0; f < a.Length; f++)
                squared += (a[f] - b[f]) * (a[f] - b[f]);
            return Math.Sqrt(squared);
        }

        private static double WeightedWindowDistance(double[] query, double[] target, double[]? variances, double epsilon)
        {
            double squared = 0;
            for (int f = 0; f < query.Length; f++)
                squared += (query[f] - target[f]) * (query[f] - target[f]) / ((variances?[f] ?? 0) + epsilon);
            return Math.Sqrt(squared);
        }

        // ZCA whitening: V diag(1 / sqrt(lambda)) V^T of the covariance of the rows.
        private static double[,] Whitening(double[][] rows)
        {
            int width = rows[0].Length;
            var means = Enumerable.Range(0, width).Select(f => rows.Average(r => r[f])).ToArray();
            var covariance = new double[width, width];
            foreach (var row in rows)
                for (int a = 0; a < width; a++)
                    for (int b = 0; b < width; b++)
                        covariance[a, b] += (row[a] - means[a]) * (row[b] - means[b]) / rows.Length;
            var (values, vectors) = Ordination.SymmetricEigen(covariance);
            var matrix = new double[width, width];
            for (int a = 0; a < width; a++)
                for (int b = 0; b < width; b++)
                    for (int c = 0; c < width; c++)
                        matrix[a, b] += vectors[a, c] * vectors[b, c] / Math.Sqrt(Math.Max(values[c], 1e-12));
            return matrix;
        }

        private static double[] Multiply(double[,] matrix, double[] vector)
        {
            var result = new double[vector.Length];
            for (int a = 0; a < vector.Length; a++)
                for (int b = 0; b < vector.Length; b++)
                    result[a] += matrix[a, b] * vector[b];
            return result;
        }

        // Survival of one query at `competitors` other-genome genes: none if its own genome already beats it; else the measured share
        // closer if any competitor in `sorted` beats it; else extrapolated from the tail below the nearest competitor.
        private static double Survive(int closerOwn, double[] sorted, double trueDistance, double competitors)
        {
            if (closerOwn > 0) return 0;
            int closer = 0;
            while (closer < sorted.Length && sorted[closer] < trueDistance) closer++;
            var share = closer > 0 ? closer / (double)sorted.Length : P.CloserShare(sorted, trueDistance, tailNeighbours);
            return Math.Exp(-competitors * share);
        }

        private sealed record HomologCase(string Target, string Protein, string QueryPos12, int TrueSlot, string[] CandidatePos12, (double X, double Y, double Z)[] TargetTable);

        // R6d's 200 K-12 queries and their orthologs (largest share of amino-acid 5-mers, at least 0.3) among each target's genes >= 1.5 kbp.
        private static List<HomologCase> HomologCases(long[] uniformUsage)
        {
            var querySequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(P.genericUsageFile)).Sequence;
            var queryReverse = P.ReverseComplement(querySequence);
            var picker = new Random(P.StableSeed(P.genericUsageFile) ^ 0x68e31da4);
            var queryGenes = P.AtgGenes(querySequence, P.standardCode, P.homologQueryLength).Select(g => querySequence.Substring(g.Start, g.End - g.Start))
                .Concat(P.AtgGenes(queryReverse, P.standardCode, P.homologQueryLength).Select(g => queryReverse.Substring(g.Start, g.End - g.Start)))
                .OrderBy(_ => picker.Next())
                .Take(P.homologQueryCount)
                .ToList();
            var cases = new List<HomologCase>();
            foreach (var target in P.homologTargets)
            {
                var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(target)).Sequence;
                var reverse = P.ReverseComplement(forward);
                var table = P.ExpectedAxesTable(P.Synonyms(P.standardCode), P.CodonUsage(forward, P.OpenStretches(forward, P.standardCode, P.codonUsageMinLength)));
                var genes = P.AtgGenes(forward, P.standardCode, P.orthologSearchMinLength).Select(g => forward.Substring(g.Start, g.End - g.Start))
                    .Concat(P.AtgGenes(reverse, P.standardCode, P.orthologSearchMinLength).Select(g => reverse.Substring(g.Start, g.End - g.Start)))
                    .ToList();
                var proteins = genes.Select(g => P.Translate(g, P.standardCode)).ToList();
                var kmerIndex = new Dictionary<int, List<int>>();
                for (int gene = 0; gene < proteins.Count; gene++)
                    foreach (var kmer in P.PeptideKmers(proteins[gene]))
                    {
                        if (!kmerIndex.TryGetValue(kmer, out var list)) kmerIndex[kmer] = list = new List<int>();
                        list.Add(gene);
                    }
                var candidates = Enumerable.Range(0, genes.Count).Where(g => genes[g].Length >= P.homologQueryLength).ToList();
                var slotOf = candidates.Select((gene, slot) => (gene, slot)).ToDictionary(p => p.gene, p => p.slot);
                var candidatePos12 = candidates.Select(g => P.CodonPositions12(genes[g][..P.homologQueryLength])).ToArray();
                foreach (var queryGene in queryGenes)
                {
                    var protein = P.Translate(queryGene, P.standardCode);
                    var queryKmers = P.PeptideKmers(protein);
                    var counts = new Dictionary<int, int>();
                    foreach (var kmer in queryKmers)
                        if (kmerIndex.TryGetValue(kmer, out var list))
                            foreach (var gene in list)
                                counts[gene] = counts.GetValueOrDefault(gene) + 1;
                    if (counts.Count == 0) continue;
                    var best = counts.MaxBy(c => c.Value);
                    if (best.Value / (double)queryKmers.Count < P.orthologMinShare || !slotOf.TryGetValue(best.Key, out var trueSlot)) continue;
                    cases.Add(new HomologCase(target, protein[..(P.homologQueryLength / 3)], P.CodonPositions12(queryGene[..P.homologQueryLength]), trueSlot, candidatePos12, table));
                }
            }
            return cases;
        }

        private static IEnumerable<string> HomologRanks(List<HomologCase> cases, string cell, int window, bool moments, int[] axes)
        {
            var uniformTable = P.ExpectedAxesTable(P.Synonyms(P.standardCode), Enumerable.Repeat(1L, 64).ToArray());
            foreach (var group in cases.GroupBy(c => c.Target))
            {
                var candidates = group.First().CandidatePos12.Select(c => Features(DnaAxes(c), searchPos12Length, window, moments, axes)).ToArray();
                foreach (var variant in new[] { "Dna", "ExpectedUniform", "ExpectedPerTarget" })
                {
                    var ranks = group.AsParallel().AsOrdered().Select(c =>
                    {
                        var query = variant switch
                        {
                            "Dna" => Features(DnaAxes(c.QueryPos12), searchPos12Length, window, moments, axes),
                            "ExpectedUniform" => Features(ExpectedAxes(c.Protein, uniformTable), searchPos12Length, window, moments, axes),
                            _ => Features(ExpectedAxes(c.Protein, c.TargetTable), searchPos12Length, window, moments, axes),
                        };
                        var trueDistance = Distance(query, candidates[c.TrueSlot]);
                        return 1 + Enumerable.Range(0, candidates.Length).Count(s => s != c.TrueSlot && Distance(query, candidates[s]) < trueDistance);
                    }).ToList();
                    yield return $"{cell},{group.Key},{variant},{ranks.Count},{candidates.Length},{Format(ranks.Average(r => r == 1 ? 1.0 : 0))},{Format(ranks.Average(r => r <= 10 ? 1.0 : 0))}";
                }
            }
        }

        /// <summary>
        /// Short proteins with windows sized to the gene (2026-11-05, your idea). Today a gene gets windows of 333 positions-1-and-2 bases
        /// (500 bp of gene), so a protein under 333 aa has at most one and isn't searched (4.3% right family at catalogue scale,
        /// results/ProteinSearch_Catalogue_1Window.csv). Instead, split a short gene's positions 1 and 2 into N equal windows (N = 3, 6, 9)
        /// with the same moments (P0-P2 of x, y, z), and compare a query only with genes within 10% of its length, since windows of
        /// different genes then differ in size. Pool: the 30 bundled bacteria, genes called as GeneLevelSignatureBenchmark does (longest
        /// ORFs first, no large overlaps) from 50 aa, every length. Queries: up to 400 genes per length bin, from the protein alone
        /// (uniform codon weights). Cells: fixed333 (today: windows of 333, candidates with at least the query's windows), fixed333-band
        /// (the same within the length band, to separate the band's effect from the split's), split3, split4, split5, split6, split9 (N = 4 and 5
        /// added on your request after the first run). Measured and
        /// extrapolated as ProteinSearchAcrossCdsWindows: exact-gene top-1 at 30 genomes, the 8-to-30 model check, the catalogue model.
        /// Rule fixed before the run: if no split cell reaches 50% at catalogue scale for 167-332 aa with a valid model, keep the 333 aa
        /// minimum; else the smallest N within 2 points of the best. Writes ShortProtein_Windows.csv.
        /// </summary>
        [Test]
        public void ShortProteinsWithGeneSizedWindows()
        {
            var bins = new (int Min, int Max)[] { (50, 99), (100, 166), (167, 332), (333, 499) };
            // split5-half and split5-byte (added 2026-11-05): split5 with the stored targets rounded to 16 bits or to 8 bits (ByteRounded).
            var cells = new[] { "fixed333", "fixed333-band", "split3", "split4", "split5", "split6", "split9", "split5-half", "split5-byte" };
            var rows = ShortProteinSelfSearchRows(50, bins, cells);

            using (var csv = File.CreateText("ShortProtein_Windows.csv"))
            {
                csv.WriteLine("ProteinAa,Cell,Queries,Searchable,MedianCandidates,Top1At30,Model8To30,ModelValid,ModelTop1AtCatalogue");
                foreach (var row in rows)
                    csv.WriteLine(row);
            }
            var splits = rows.Select(r => r.Split(',')).Where(r => r[0] == "167-332" && r[1].StartsWith("split", StringComparison.Ordinal) && !r[1].Contains('-') && r[7] == "true")
                .Select(r => (Cell: r[1], Value: double.Parse(r[^1], CultureInfo.InvariantCulture))).ToList();
            double best = splits.Count == 0 ? 0 : splits.Max(s => s.Value);
            var chosen = splits.FirstOrDefault(s => s.Value >= best - 0.02);
            TestContext.Out.WriteLine(best >= 0.5
                ? $"Rule: {chosen.Cell} at {Format(chosen.Value)} for 167-332 aa at catalogue scale (best {Format(best)}) -> build card"
                : $"Rule: best valid split cell {Format(best)} for 167-332 aa at catalogue scale -> stop, keep the 333 aa minimum");
            Assert.That(rows, Has.Count.EqualTo(cells.Length * bins.Length));
        }

        // Genes are called from minProtein aa, so the pool also holds every gene shorter than the first bin that a band can reach.
        private static List<string> ShortProteinSelfSearchRows(int minProtein, (int Min, int Max)[] bins, string[] cells, int queriesPerBin = 400)
        {
            const double lengthBand = 0.1;
            const int fixedWindow = 333;
            var xyz = new[] { 0, 1, 2 };
            var uniformUsage = Enumerable.Repeat(1L, 64).ToArray();
            var fileNames = P.Corpus().Where(c => c.Group == "Bacteria").Select(c => c.FileName).ToList();
            int genomeCount = fileNames.Count;
            var pool = new List<(int Genome, string Protein, string Pos12)>();
            var uniformTables = new List<(double X, double Y, double Z)[]>();
            for (int genome = 0; genome < genomeCount; genome++)
            {
                var forward = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(fileNames[genome])).Sequence;
                var code = Path.GetFileName(fileNames[genome]).StartsWith("mycoplasma_", StringComparison.Ordinal) ? P.mycoplasmaCode : P.standardCode;
                uniformTables.Add(P.ExpectedAxesTable(P.Synonyms(code), uniformUsage));
                foreach (var gene in GeneLevelSignatureBenchmark.CallGenes(forward, code, 3 * minProtein))
                    pool.Add((genome, P.Translate(gene, code), P.CodonPositions12(gene)));
            }
            var picker = new Random(20261105);
            var queries = new List<(int Bin, int Gene)>();
            for (int b = 0; b < bins.Length; b++)
                queries.AddRange(Enumerable.Range(0, pool.Count).Where(g => pool[g].Protein.Length >= bins[b].Min && pool[g].Protein.Length <= bins[b].Max)
                    .OrderBy(_ => picker.Next()).Take(queriesPerBin).Select(g => (b, g)));
            TestContext.Out.WriteLine($"{pool.Count} genes in {genomeCount} genomes; queries per bin: " +
                string.Join(", ", Enumerable.Range(0, bins.Length).Select(b => $"{bins[b].Min}-{bins[b].Max} aa {queries.Count(q => q.Bin == b)}")));

            int WindowOf(string cell, int pos12Length) => cell.StartsWith("fixed", StringComparison.Ordinal)
                ? fixedWindow : Math.Max(1, pos12Length / int.Parse(cell.Split('-')[0]["split".Length..], CultureInfo.InvariantCulture));
            // Survive needs two competitors for its tail; with fewer, only what was measured counts.
            double Model(int closerOwn, double[] sorted, double trueDistance, double competitors) => sorted.Length >= 2
                ? Survive(closerOwn, sorted, trueDistance, competitors)
                : closerOwn == 0 && (sorted.Length == 0 || sorted[0] >= trueDistance) ? 1 : 0;

            var rows = new List<string>();
            foreach (var cell in cells)
            {
                var targets = pool.Select(g => Features(DnaAxes(g.Pos12), g.Pos12.Length, WindowOf(cell, g.Pos12.Length), true, xyz)).ToArray();
                if (cell.EndsWith("-half", StringComparison.Ordinal)) targets = HalfRounded(targets);
                if (cell.EndsWith("-byte", StringComparison.Ordinal)) targets = ByteRounded(targets, pool.Select(g => g.Genome).ToList());
                bool banded = cell != "fixed333";
                var results = new (bool Searchable, double Observed, double From8, double Catalogue, int Candidates)[queries.Count];
                Parallel.For(0, queries.Count, j =>
                {
                    var (_, truth) = queries[j];
                    var (own, protein, pos12) = pool[truth];
                    var query = Features(ExpectedAxes(protein, uniformTables[own]), pos12.Length, WindowOf(cell, pos12.Length), true, xyz);
                    if (query.Length == 0) return;
                    var trueDistance = Distance(query, targets[truth]);
                    int closerOwn = 0;
                    var others = new List<(double Distance, int Genome)>();
                    for (int i = 0; i < pool.Count; i++)
                    {
                        if (i == truth || targets[i].Length < Math.Min(query.Length, 3)) continue;
                        if (banded && Math.Abs(pool[i].Protein.Length - protein.Length) > Tolerance(cell, protein.Length, lengthBand)) continue;
                        var distance = Distance(query, targets[i]);
                        if (pool[i].Genome == own)
                        {
                            if (distance < trueDistance) closerOwn++;
                        }
                        else others.Add((distance, pool[i].Genome));
                    }
                    var rng = new Random(P.StableSeed(protein) ^ P.StableSeed(cell));
                    var fitSet = Enumerable.Range(0, genomeCount).Where(g => g != own).OrderBy(_ => rng.Next()).Take(fitGenomes - 1).ToHashSet();
                    var all = others.Select(o => o.Distance).Order().ToArray();
                    var fit = others.Where(o => fitSet.Contains(o.Genome)).Select(o => o.Distance).Order().ToArray();
                    results[j] = (true, closerOwn == 0 && (all.Length == 0 || all[0] >= trueDistance) ? 1 : 0,
                        Model(closerOwn, fit, trueDistance, all.Length),
                        Model(closerOwn, all, trueDistance, catalogueGenomes / (double)genomeCount * all.Length),
                        all.Length);
                });
                for (int b = 0; b < bins.Length; b++)
                {
                    var inBin = Enumerable.Range(0, queries.Count).Where(j => queries[j].Bin == b).Select(j => results[j]).ToList();
                    var searchable = inBin.Where(r => r.Searchable).ToList();
                    // An unsearchable query (no window) counts as a miss, as the API refuses it.
                    double observed = inBin.Sum(r => r.Observed) / inBin.Count, from8 = inBin.Sum(r => r.From8) / inBin.Count;
                    double searchableObserved = searchable.Count == 0 ? 0 : searchable.Average(r => r.Observed);
                    double searchableFrom8 = searchable.Count == 0 ? 0 : searchable.Average(r => r.From8);
                    bool valid = searchable.Count > 0 && Math.Abs(searchableFrom8 - searchableObserved) <= modelTolerance;
                    var candidates = searchable.Select(r => r.Candidates).Order().ToList();
                    var line = string.Join(",", $"{bins[b].Min}-{bins[b].Max}", cell, inBin.Count, searchable.Count,
                        candidates.Count == 0 ? 0 : candidates[candidates.Count / 2], Format(observed), Format(from8),
                        valid.ToString().ToLowerInvariant(), Format(inBin.Sum(r => r.Catalogue) / inBin.Count));
                    rows.Add(line);
                    TestContext.Out.WriteLine(line);
                }
            }
            return rows;
        }

        /// <summary>
        /// Homologs with gene-sized windows (2026-11-05). ShortProteinsWithGeneSizedWindows measured exact self-search, where the gene is
        /// exactly as long as the query; an ortholog is often a few residues longer or shorter, which moves every gene-sized window. Queries:
        /// E. coli K-12 genes called from 50 aa (CallGenes), up to 300 per length bin, from the protein alone (uniform weights). In each of
        /// R6f's four targets, the ortholog is the target gene sharing the largest share of amino-acid 5-mers, at least 0.3 (R6d's rule),
        /// and the candidates are all the target's genes. Cells: fixed333 (today: candidates with at least min(query windows, 3)), split3,
        /// split6 and split9 within 10% of the query's length (split4 and split5 added after the first run), and split6 within 20%. Measured: ortholog top-1 and top-10 (an ortholog the
        /// cell can't reach, outside the band or without the windows, is a miss), and the share of orthologs inside the band.
        /// Rule fixed before the run: split6 (10%) replaces the fixed windows for every gene under 500 aa if its mean ortholog top-1 over the
        /// four targets at 333-499 aa is within 5 points of fixed333's; otherwise only under 333 aa, where the fixed windows can't search.
        /// Writes ShortProtein_Homologs.csv.
        /// </summary>
        [Test]
        public void ShortProteinHomologsWithGeneSizedWindows()
        {
            var bins = new (int Min, int Max)[] { (50, 99), (100, 166), (167, 332), (333, 499) };
            var cells = new (string Name, int Split, double Band)[]
            {
                ("fixed333", 0, double.PositiveInfinity), ("split3", 3, 0.1), ("split4", 4, 0.1), ("split5", 5, 0.1), ("split6", 6, 0.1),
                ("split9", 9, 0.1), ("split6-band20", 6, 0.2), ("split5-half", 5, 0.1), ("split5-byte", 5, 0.1),
            };
            var rows = ShortProteinHomologRows(50, bins, cells);
            WriteHomologCsv("ShortProtein_Homologs.csv", rows);
            double fixedTop1 = MeanHomologTop1(rows, "333-499", "fixed333"), splitTop1 = MeanHomologTop1(rows, "333-499", "split6");
            TestContext.Out.WriteLine($"Rule: at 333-499 aa split6 {Format(splitTop1)} vs fixed333 {Format(fixedTop1)} -> " +
                (splitTop1 >= fixedTop1 - 0.05 ? "split every gene under 500 aa" : "split only under 333 aa"));
            Assert.That(rows, Is.Not.Empty);
        }

        private static List<string> ShortProteinHomologRows(int minProtein, (int Min, int Max)[] bins, (string Name, int Split, double Band)[] cells,
            int queriesPerBin = 300)
        {
            const int fixedWindow = 333;
            var xyz = new[] { 0, 1, 2 };
            var uniformTable = P.ExpectedAxesTable(P.Synonyms(P.standardCode), Enumerable.Repeat(1L, 64).ToArray());
            var queryGenes = GeneLevelSignatureBenchmark.CallGenes(FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(P.genericUsageFile)).Sequence,
                P.standardCode, 3 * minProtein);
            var queryProteins = queryGenes.Select(g => P.Translate(g, P.standardCode)).ToList();
            var picker = new Random(20261105);
            var queries = new List<(int Bin, string Protein)>();
            for (int b = 0; b < bins.Length; b++)
                queries.AddRange(queryProteins.Where(p => p.Length >= bins[b].Min && p.Length <= bins[b].Max)
                    .OrderBy(_ => picker.Next()).Take(queriesPerBin).Select(p => (b, p)));
            int WindowOf(int split, int pos12Length) => split == 0 ? fixedWindow : Math.Max(1, pos12Length / split);

            var rows = new List<string>();
            foreach (var target in P.homologTargets)
            {
                var genes = GeneLevelSignatureBenchmark.CallGenes(FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(target)).Sequence, P.standardCode, 3 * minProtein);
                var proteins = genes.Select(g => P.Translate(g, P.standardCode)).ToList();
                var pos12 = genes.Select(P.CodonPositions12).ToList();
                var kmerIndex = new Dictionary<int, List<int>>();
                for (int gene = 0; gene < proteins.Count; gene++)
                    foreach (var kmer in P.PeptideKmers(proteins[gene]))
                    {
                        if (!kmerIndex.TryGetValue(kmer, out var list)) kmerIndex[kmer] = list = new List<int>();
                        list.Add(gene);
                    }
                var orthologs = queries.Select(q =>
                {
                    var queryKmers = P.PeptideKmers(q.Protein);
                    var counts = new Dictionary<int, int>();
                    foreach (var kmer in queryKmers)
                        if (kmerIndex.TryGetValue(kmer, out var list))
                            foreach (var gene in list)
                                counts[gene] = counts.GetValueOrDefault(gene) + 1;
                    if (counts.Count == 0) return -1;
                    var best = counts.MaxBy(c => c.Value);
                    return best.Value / (double)queryKmers.Count >= P.orthologMinShare ? best.Key : -1;
                }).ToArray();

                foreach (var (name, split, band) in cells)
                {
                    var targets = pos12.Select(p => Features(DnaAxes(p), p.Length, WindowOf(split, p.Length), true, xyz)).ToArray();
                    if (name.EndsWith("-half", StringComparison.Ordinal)) targets = HalfRounded(targets);
                    if (name.EndsWith("-byte", StringComparison.Ordinal)) targets = ByteRounded(targets, new int[targets.Length]);
                    var outcomes = new (bool InBand, int Rank)[queries.Count];
                    Parallel.For(0, queries.Count, j =>
                    {
                        int truth = orthologs[j];
                        if (truth < 0) return;
                        var protein = queries[j].Protein;
                        var query = Features(ExpectedAxes(protein, uniformTable), 2 * protein.Length, WindowOf(split, 2 * protein.Length), true, xyz);
                        bool Eligible(int i) => targets[i].Length >= Math.Min(query.Length, 3) && Math.Abs(proteins[i].Length - protein.Length) <= Tolerance(name, protein.Length, band);
                        if (query.Length == 0 || !Eligible(truth))
                        {
                            outcomes[j] = (false, int.MaxValue);
                            return;
                        }
                        var trueDistance = Distance(query, targets[truth]);
                        int rank = 1;
                        for (int i = 0; i < genes.Count; i++)
                            if (i != truth && Eligible(i) && Distance(query, targets[i]) < trueDistance) rank++;
                        outcomes[j] = (true, rank);
                    });
                    for (int b = 0; b < bins.Length; b++)
                    {
                        var cases = Enumerable.Range(0, queries.Count).Where(j => queries[j].Bin == b && orthologs[j] >= 0).Select(j => outcomes[j]).ToList();
                        if (cases.Count == 0) continue;
                        var line = string.Join(",", $"{bins[b].Min}-{bins[b].Max}", name, Path.GetFileName(target).Split('_')[0] + "_" + Path.GetFileName(target).Split('_')[1],
                            cases.Count, Format(cases.Average(c => c.InBand ? 1.0 : 0)), Format(cases.Average(c => c.Rank == 1 ? 1.0 : 0)),
                            Format(cases.Average(c => c.Rank <= 10 ? 1.0 : 0)));
                        rows.Add(line);
                        TestContext.Out.WriteLine(line);
                    }
                }
            }
            return rows;
        }

        private static void WriteHomologCsv(string fileName, List<string> rows)
        {
            using (var csv = File.CreateText(fileName))
            {
                csv.WriteLine("ProteinAa,Cell,Target,Orthologs,OrthologReachable,Top1,Top10");
                foreach (var row in rows)
                    csv.WriteLine(row);
                csv.WriteLine();
                csv.WriteLine("ProteinAa,Cell,MeanTop1OverTargets,MeanTop10OverTargets");
                foreach (var group in rows.Select(r => r.Split(',')).GroupBy(r => (Bin: r[0], Cell: r[1])))
                    csv.WriteLine(string.Join(",", group.Key.Bin, group.Key.Cell,
                        Format(group.Average(r => double.Parse(r[5], CultureInfo.InvariantCulture))),
                        Format(group.Average(r => double.Parse(r[6], CultureInfo.InvariantCulture)))));
            }
        }

        // How many residues a candidate's length may differ from the query's: "-pm5" in the cell name means 5 aa, "-band20" 20% of the
        // query's length; any other cell keeps the caller's relative band (so the committed cells compute what they did).
        private static double Tolerance(string cell, int residues, double band)
        {
            if (Regex.Match(cell, @"-pm(\d+)") is { Success: true } absolute) return int.Parse(absolute.Groups[1].Value, CultureInfo.InvariantCulture);
            if (Regex.Match(cell, @"-band(\d+)") is { Success: true } relative) return int.Parse(relative.Groups[1].Value, CultureInfo.InvariantCulture) / 100.0 * residues;
            return band * residues;
        }

        private static double MeanHomologTop1(List<string> rows, string bin, string cell) => rows.Select(r => r.Split(',')).Where(r => r[0] == bin && r[1] == cell)
            .Average(r => double.Parse(r[5], CultureInfo.InvariantCulture));

        /// <summary>
        /// Proteins under 50 aa (2026-11-09, your idea: two windows instead of leaving them out). The platform stores no gene under 50 aa
        /// (ZCurveGeneProfile.MinResidues; 4.6% of BV-BRC's genes, 12.7% of one plant annotation), a floor set by the S7.S card, never measured.
        /// Both short-protein checks above, with genes called from 20 aa (the ORF caller places longer ORFs first, so the pool of 50 aa and up
        /// is unchanged and short ORFs fill the gaps; an ORF caller over-calls short genes, which makes the pool harder than an annotation).
        /// Bins 20-34 and 35-49 aa, and 50-99 aa as the anchor; cells split2, split3 and split5 within 10% of the query's length.
        /// Rule fixed before the run: a bin is searchable if some valid cell reaches 90% exact-gene top-1 at catalogue scale (the weakSignal
        /// line) and that cell's mean ortholog top-1 is at least 70%; the cell is the smallest N within 2 points of the best on both. If neither
        /// bin passes, the 50 aa floor stays. Writes ShortProtein_Under50.csv and ShortProtein_Under50_Homologs.csv.
        /// </summary>
        [Test]
        public void ProteinsUnder50WithFewWindows()
        {
            var bins = new (int Min, int Max)[] { (20, 34), (35, 49), (50, 99) };
            var selfRows = ShortProteinSelfSearchRows(20, bins, ["split2", "split3", "split5"]);
            using (var csv = File.CreateText("ShortProtein_Under50.csv"))
            {
                csv.WriteLine("ProteinAa,Cell,Queries,Searchable,MedianCandidates,Top1At30,Model8To30,ModelValid,ModelTop1AtCatalogue");
                foreach (var row in selfRows)
                    csv.WriteLine(row);
            }
            var homologRows = ShortProteinHomologRows(20, bins, [("split2", 2, 0.1), ("split3", 3, 0.1), ("split5", 5, 0.1)]);
            WriteHomologCsv("ShortProtein_Under50_Homologs.csv", homologRows);

            foreach (var bin in new[] { "20-34", "35-49", "50-99" })
            {
                var cells = selfRows.Select(r => r.Split(',')).Where(r => r[0] == bin).Select(r => (Cell: r[1], Valid: r[7] == "true",
                    Exact: double.Parse(r[^1], CultureInfo.InvariantCulture), Ortholog: MeanHomologTop1(homologRows, bin, r[1]))).ToList();
                foreach (var c in cells)
                    TestContext.Out.WriteLine($"{bin} {c.Cell}: exact at catalogue scale {Format(c.Exact)} (model valid {c.Valid}), ortholog {Format(c.Ortholog)}");
                var passing = cells.Where(c => c.Valid && c.Exact >= 0.9 && c.Ortholog >= 0.7).ToList();
                if (passing.Count == 0)
                {
                    TestContext.Out.WriteLine($"Rule: {bin} aa -> not searchable");
                    continue;
                }
                double bestExact = passing.Max(c => c.Exact), bestOrtholog = passing.Max(c => c.Ortholog);
                var chosen = passing.FirstOrDefault(c => c.Exact >= bestExact - 0.02 && c.Ortholog >= bestOrtholog - 0.02);
                TestContext.Out.WriteLine(chosen.Cell is null
                    ? $"Rule: {bin} aa -> searchable, but no cell is within 2 points of the best on both; pick by hand"
                    : $"Rule: {bin} aa -> searchable with {chosen.Cell}");
            }
            Assert.That(selfRows, Has.Count.EqualTo(3 * bins.Length));
        }

        /// <summary>
        /// Wider length bands for proteins under 50 aa (2026-11-09, your go). ProteinsUnder50WithFewWindows failed 35-49 aa only on
        /// orthologs (69.5% against 70%), and almost every miss was an ortholog outside the ±10% band. Same pools and five windows, with
        /// 1,000 queries per bin (400 and 300 before) for more orthologs: split5 at ±10% (today), ±20% (split5-band20) and ±5 aa (split5-pm5).
        /// Rule fixed before the run: a band passes if, in both 20-34 and 35-49 aa, it reaches 90% exact-gene top-1 at catalogue scale
        /// (valid model) and 70% mean ortholog top-1; of the passing bands, the one with the fewest median candidates over the two bins. If
        /// none passes, the 50 aa floor stays and the question is closed. Writes ShortProtein_Under50Bands.csv and
        /// ShortProtein_Under50Bands_Homologs.csv.
        /// </summary>
        [Test]
        public void ProteinsUnder50WithWiderBands()
        {
            const int queriesPerBin = 1000;
            var bins = new (int Min, int Max)[] { (20, 34), (35, 49) };
            string[] cells = ["split5", "split5-band20", "split5-pm5"];
            var selfRows = ShortProteinSelfSearchRows(20, bins, cells, queriesPerBin);
            using (var csv = File.CreateText("ShortProtein_Under50Bands.csv"))
            {
                csv.WriteLine("ProteinAa,Cell,Queries,Searchable,MedianCandidates,Top1At30,Model8To30,ModelValid,ModelTop1AtCatalogue");
                foreach (var row in selfRows)
                    csv.WriteLine(row);
            }
            var homologRows = ShortProteinHomologRows(20, bins, [("split5", 5, 0.1), ("split5-band20", 5, 0.2), ("split5-pm5", 5, 0.1)], queriesPerBin);
            WriteHomologCsv("ShortProtein_Under50Bands_Homologs.csv", homologRows);

            var passing = new List<(string Cell, double Candidates)>();
            foreach (var cell in cells)
            {
                var rows = selfRows.Select(r => r.Split(',')).Where(r => r[1] == cell).ToList();
                bool passes = true;
                foreach (var r in rows)
                {
                    double exact = double.Parse(r[^1], CultureInfo.InvariantCulture), ortholog = MeanHomologTop1(homologRows, r[0], cell);
                    TestContext.Out.WriteLine($"{r[0]} {cell}: exact at catalogue scale {Format(exact)} (model valid {r[7]}), ortholog {Format(ortholog)}, median candidates {r[4]}");
                    passes &= r[7] == "true" && exact >= 0.9 && ortholog >= 0.7;
                }
                if (passes) passing.Add((cell, rows.Average(r => double.Parse(r[4], CultureInfo.InvariantCulture))));
            }
            TestContext.Out.WriteLine(passing.Count == 0
                ? "Rule: no band passes in both bins -> the 50 aa floor stays, question closed"
                : $"Rule: passes with {passing.MinBy(p => p.Candidates).Cell} (passing: {string.Join(", ", passing.Select(p => p.Cell))})");
            Assert.That(selfRows, Has.Count.EqualTo(cells.Length * bins.Length));
        }

        [Test]
        public void TreesAcrossCdsWindows()
        {
            using var report = File.CreateText("CdsSweep_Trees.csv");
            report.WriteLine("Corpus, Distance, Rank, Positives, Negatives, Auc, Taxa, NjClades, UpgmaClades, NjNullMean, UpgmaNullMean, NjP, UpgmaP");
            var bacteria = GenomeTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.OrganismName}|{g.Phylum}", g.Genus, Higher: g.Phylum)).ToList();
            var viruses = ViralTaxonomy.Genomes.Select(g => (g.FileName, Label: $"{g.Name}|{g.Family}", g.Genus, Higher: g.Family)).ToList();
            var results = new Dictionary<(string Corpus, string Distance, string Rank), (double Auc, int Clades)>();
            var candidates = new List<string>();
            foreach (var (corpus, genomes, higher) in new[] { ("bacteria", bacteria, "phylum"), ("viruses", viruses, "family") })
            {
                int n = genomes.Count;
                var prefixes = new string[n][];
                var counts = new long[n][];
                var cvTree4 = new double[n][];
                for (int i = 0; i < n; i++)
                {
                    var sequence = FastaReader.ReadPrimaryRecord(GenomeCatalog.PathTo(genomes[i].FileName)).Sequence;
                    var code = Path.GetFileName(genomes[i].FileName).StartsWith("mycoplasma_", StringComparison.Ordinal) ? P.mycoplasmaCode : P.standardCode;
                    var genes = GeneLevelSignatureBenchmark.CallGenes(sequence, code);
                    counts[i] = GenomeSignatures.DinucleotideCounts(sequence);
                    cvTree4[i] = GeneLevelSignatureBenchmark.CvTreeVector(genes.Select(g => P.Translate(g, code)).ToList(), 4);
                    prefixes[i] = genes.Where(g => g.Length >= treeLength).Select(g => g[..treeLength]).ToArray();
                }
                var distances = new List<(string Name, double[,] Matrix)>
                {
                    ("karlin", GeneLevelSignatureBenchmark.Matrix(n, (a, b) => GenomeSignatures.KarlinDelta(GenomeSignatures.KarlinRho(counts[a]), GenomeSignatures.KarlinRho(counts[b])))),
                    ("cvTree4", GeneLevelSignatureBenchmark.Matrix(n, (a, b) => GeneLevelSignatureBenchmark.CosineDistance(cvTree4[a], cvTree4[b]))),
                };
                foreach (var (axisName, axes) in axisSets)
                    foreach (var window in treeWindows)
                        foreach (var moments in new[] { false, true })
                        {
                            var name = $"cds_{axisName}_{window}_{(moments ? "moments" : "means")}";
                            if (!candidates.Contains(name)) candidates.Add(name);
                            var features = prefixes.Select(p => p.Select(gene => Features(DnaAxes(gene), treeLength, window, moments, axes)).ToArray()).ToArray();
                            distances.Add((name, BestHitMatrix(features)));
                            TestContext.Out.WriteLine($"  {corpus} {name} done");
                        }
                GeneLevelSignatureBenchmark.ScoreTrees(corpus, genomes, higher, distances, "CdsSweep", report, results);
            }

            GeneLevelSignatureBenchmark.ReportRule(candidates, results);
            var cvPhylum = results[("bacteria", "cvTree4", "phylum")].Auc;
            var cvClades = results[("bacteria", "cvTree4", "genus")].Clades;
            var beating = candidates.Where(c => results[("bacteria", c, "phylum")].Auc > cvPhylum && results[("bacteria", c, "genus")].Clades > cvClades).ToList();
            TestContext.Out.WriteLine(string.Create(CultureInfo.InvariantCulture,
                $"Tree rule (cvTree4 phylum AUC {cvPhylum:F4}, genus clades {cvClades}): {(beating.Count > 0 ? "pass for " + string.Join(", ", beating) : "stop")}"));
        }

        // One genome pair: each gene's best distance among the other genome's genes, averaged both ways; no genes -> the largest distance seen.
        private static double[,] BestHitMatrix(double[][][][] features)
        {
            int n = features.Length;
            var matrix = new double[n, n];
            Parallel.For(0, n * n, cell =>
            {
                int a = cell / n, b = cell % n;
                if (b <= a) return;
                if (features[a].Length == 0 || features[b].Length == 0)
                {
                    matrix[a, b] = matrix[b, a] = double.NaN;
                    return;
                }
                var bestA = Enumerable.Repeat(double.PositiveInfinity, features[a].Length).ToArray();
                var bestB = Enumerable.Repeat(double.PositiveInfinity, features[b].Length).ToArray();
                for (int i = 0; i < features[a].Length; i++)
                    for (int j = 0; j < features[b].Length; j++)
                    {
                        var distance = Distance(features[a][i], features[b][j]);
                        if (distance < bestA[i]) bestA[i] = distance;
                        if (distance < bestB[j]) bestB[j] = distance;
                    }
                matrix[a, b] = matrix[b, a] = (bestA.Average() + bestB.Average()) / 2;
            });
            double worst = 0;
            for (int a = 0; a < n; a++)
                for (int b = a + 1; b < n; b++)
                    if (double.IsFinite(matrix[a, b])) worst = Math.Max(worst, matrix[a, b]);
            for (int a = 0; a < n; a++)
                for (int b = 0; b < n; b++)
                    if (a != b && !double.IsFinite(matrix[a, b])) matrix[a, b] = worst;
            return matrix;
        }

        private static Func<int, (double X, double Y, double Z)> DnaAxes(string bases) => i =>
        {
            var (x, y, z) = P.BaseAxes(bases[i]);
            return (x, y, z);
        };

        // Each positions-1-and-2 base's expected axes from the protein: codon position 1 then 2 of every residue.
        private static Func<int, (double X, double Y, double Z)> ExpectedAxes(string protein, (double X, double Y, double Z)[] table) =>
            i => table[2 * (protein[i / 2] - 'A') + i % 2];

        // Per complete window over the first `length` bases: for each kept axis, sqrt(2k + 1) * mean(a_i P_k(t_i)) for k = 0 (means) or
        // k = 0, 1, 2 (moments), where t_i is the base's position in [-1, 1]; sqrt(2k + 1) gives every P_k unit mean square. The length
        // is the prefix profiled: searchPos12Length positions-1-and-2 bases for search, treeLength bases for trees.
        private static double[][] Features(Func<int, (double X, double Y, double Z)> axesAt, int length, int window, bool moments, int[] axes) =>
            FeaturesOfOrder(axesAt, length, window, moments ? 3 : 1, axes);

        private static double[][] FeaturesOfOrder(Func<int, (double X, double Y, double Z)> axesAt, int length, int window, int orders, int[] axes)
        {
            var result = new double[length / window][];
            for (int w = 0; w < result.Length; w++)
            {
                var features = new double[axes.Length * orders];
                for (int i = 0; i < window; i++)
                {
                    var (x, y, z) = axesAt(w * window + i);
                    double t = -1 + 2 * (i + 0.5) / window;
                    for (int k = 0; k < orders; k++)
                    {
                        double weight = k switch { 0 => 1, 1 => t, _ => (3 * t * t - 1) / 2 };
                        for (int a = 0; a < axes.Length; a++)
                            features[axes.Length * k + a] += (axes[a] switch { 0 => x, 1 => y, _ => z }) * weight;
                    }
                }
                for (int k = 0; k < orders; k++)
                    for (int a = 0; a < axes.Length; a++)
                        features[axes.Length * k + a] *= Math.Sqrt(2 * k + 1) / window;
                result[w] = features;
            }
            return result;
        }

        // Mean over the common windows of the per-window Euclidean distance: 1 - R6's score.
        private static double Distance(double[][] a, double[][] b)
        {
            int count = Math.Min(a.Length, b.Length);
            double sum = 0;
            for (int w = 0; w < count; w++)
            {
                double squared = 0;
                for (int f = 0; f < a[w].Length; f++)
                    squared += (a[w][f] - b[w][f]) * (a[w][f] - b[w][f]);
                sum += Math.Sqrt(squared);
            }
            return sum / count;
        }

        private static string Format(double value) => value.ToString("F4", CultureInfo.InvariantCulture);
    }
}
