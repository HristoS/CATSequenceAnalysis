using SequenceAnalyses;
using SequencePlatform.Core;

namespace SequencePlatform.Api
{
    /// <summary>A search accepted for processing. Poll <c>Href</c> until its status is Done or Failed.</summary>
    public sealed record SearchAcceptedResponse(Guid Id, SearchStatus Status, string Href);

    /// <summary>What the search was asked: its catalogue filters, top K and kind-specific parameters.</summary>
    /// <param name="LengthRatio">Find-similar: candidates are records of length [ratio × query, query / ratio].</param>
    /// <param name="WithinSequenceId">Localize inside this one record only.</param>
    public sealed record SearchParametersResponse(string? Kind, string? Rank, string? Taxon, long? MinLengthBp, long? MaxLengthBp, string? Provider,
        string? Name, int Top, double? LengthRatio, long? WithinSequenceId);

    /// <summary>The shape of a finished search's result.</summary>
    /// <param name="WindowSize">Find-similar: the level every hit is scored at; localize: the coarse scan level.</param>
    /// <param name="BandMinBp">Find-similar: the shortest candidate length.</param>
    /// <param name="BandMaxBp">Find-similar: the longest candidate length; absent when unbounded.</param>
    /// <param name="RankedTotal">Find-similar: candidates ranked before the top K was taken.</param>
    /// <param name="Candidates">Localize: records scanned.</param>
    /// <param name="WeakSignal">The query is short enough that the ranking reflects composition more than order (find-similar under
    /// about 8 kbp, localize under 8 kbp): treat it as a hint, and check it with an aligner.</param>
    /// <param name="NeighbourRanked">Find-similar: the top hits are ordered by shared neighbours, then score.</param>
    public sealed record SearchSummaryResponse(int? WindowSize, long? BandMinBp, long? BandMaxBp, int? RankedTotal, int? Candidates,
        bool? WeakSignal, bool? NeighbourRanked);

    /// <summary>One hit. Which fields are present depends on the list and the kind of search.</summary>
    /// <param name="Position">1-based order within its list; for the ranked list, the rank.</param>
    /// <param name="Sequence">The record; absent if it has been deleted since.</param>
    /// <param name="Strand">Reverse: the query's reverse complement is what matches the record's stored strand.</param>
    /// <param name="CoarseScore">Localize: the score at the coarse level that brought the record into the refined top.</param>
    /// <param name="OffsetBp">Localize: where the fragment starts on the record's stored strand, to about one 500 bp window.</param>
    /// <param name="Log10P">log10 of the chance an unrelated record of this shape scores as high.</param>
    /// <param name="PKind">Empirical, or Bound: the null sample has too few scores this high, so p is only an upper bound.</param>
    /// <param name="Log10E">log10 of the expected number of unrelated candidates scoring as high.</param>
    /// <param name="SharedNeighbours">Find-similar: how many of the query's top 20 are also among this hit's own top 20.</param>
    /// <param name="Rotated">Find-similar: the best score came with the query or the record laid out from its replication origin.</param>
    public sealed record SearchHitResponse(int Position, long SequenceId, SequenceSummaryResponse? Sequence, Strand? Strand, double? Score,
        double? CoarseScore, long? OffsetBp, double? Log10P, PValueKind? PKind, double? Log10E, int? SharedNeighbours, bool Rotated);

    /// <summary>One page of one hit list.</summary>
    /// <param name="List">Ranked; Identical (records whose sequence equals the query's on either strand, by SHA-256); or NotComparable
    /// (candidates with no stored level at the search's window size).</param>
    /// <param name="NextAfter">Pass as <c>after</c>, with the same <c>list</c>, for the next page; absent on the last page.</param>
    public sealed record SearchHitPageResponse(HitList List, IReadOnlyList<SearchHitResponse> Items, int? NextAfter);

    /// <summary>
    /// A search and, once done, a page of its results. Results are kept 30 days (D7), or until DELETE; after that the search is
    /// Expired and keeps only its parameters, the query's length and hash, and the summary.
    /// </summary>
    public sealed record SearchRunResponse(Guid Id, SearchKind Kind, SearchStatus Status, string? Error, long QueryLengthBp,
        DateTimeOffset CreatedAt, DateTimeOffset? StartedAt, DateTimeOffset? CompletedAt, DateTimeOffset ExpiresAt,
        SearchParametersResponse Parameters, SearchSummaryResponse? Summary, SearchHitPageResponse? Hits)
    {
        internal static SearchRunResponse From(SearchRun run, SearchHitPageResponse? hits)
        {
            var r = run.Request;
            var f = r.Filter;
            return new SearchRunResponse(run.PublicId, r.Kind, run.Status, run.Error, run.QueryLengthBp, run.CreatedAt, run.StartedAt,
                run.CompletedAt, run.ExpiresAt,
                new SearchParametersResponse(f.Kind, f.Rank, f.Taxon, f.MinLengthBp, f.MaxLengthBp, f.Provider, f.NameContains, r.Top,
                    r.LengthRatio, r.WithinSequenceId),
                run.CompletedAt is null || run.Status == SearchStatus.Failed ? null : new SearchSummaryResponse(run.WindowSize, run.Lengths?.Min,
                    run.Lengths?.Max - 1, run.RankedTotal, run.Candidates, run.WeakSignal, run.NeighbourRanked),
                hits);
        }
    }
}
