using System;

namespace SequenceAnalyses
{
    public static class NeedlemanWunsch
    {
        private const int matchScore = 2;
        private const int mismatchScore = -1;
        private const int gapScore = -2;

        //http://www.codeproject.com/Tips/638377/Needleman-Wunsch-Algorithm-in-Csharp
        public static double Calculate(string refSeq, string alignSeq)
        {
            // Keep the shorter sequence as the "column" dimension so the two rolling
            // rows below cost O(min(len)) instead of O(max(len)); the DP score itself
            // is symmetric in ref/align, but swapping transposes the matrix, so the
            // left/up tie-break priority below must be flipped to match to keep the
            // exact same chosen alignment (and match count) as the un-swapped case.
            bool swapped = refSeq.Length > alignSeq.Length;
            if (swapped)
            {
                (refSeq, alignSeq) = (alignSeq, refSeq);
            }

            int refSeqCnt = refSeq.Length + 1;
            int alignSeqCnt = alignSeq.Length + 1;

            int[] prevScore = new int[refSeqCnt];
            int[] currScore = new int[refSeqCnt];
            int[] prevMatches = new int[refSeqCnt];
            int[] currMatches = new int[refSeqCnt];

            for (int j = 0; j < refSeqCnt; j++)
            {
                prevScore[j] = j * gapScore;
            }

            for (int i = 1; i < alignSeqCnt; i++)
            {
                currScore[0] = i * gapScore;
                currMatches[0] = 0;

                char alignChar = alignSeq[i - 1];

                for (int j = 1; j < refSeqCnt; j++)
                {
                    bool isMatch = refSeq[j - 1] == alignChar;
                    int scoreDiag = prevScore[j - 1] + (isMatch ? matchScore : mismatchScore);
                    int scoreLeft = currScore[j - 1] + gapScore;
                    int scoreUp = prevScore[j] + gapScore;

                    // Tie-break priority mirrors what a backward traceback over the
                    // (possibly transposed) full matrix would pick: diag first, then
                    // whichever of left/up corresponds to the original "left" move.
                    bool diagWins = swapped
                        ? scoreDiag >= scoreUp && scoreDiag >= scoreLeft
                        : scoreDiag >= scoreLeft && scoreDiag >= scoreUp;

                    if (diagWins)
                    {
                        currScore[j] = scoreDiag;
                        currMatches[j] = prevMatches[j - 1] + (isMatch ? 1 : 0);
                    }
                    else if (swapped ? scoreUp >= scoreLeft : scoreLeft >= scoreUp)
                    {
                        currScore[j] = swapped ? scoreUp : scoreLeft;
                        currMatches[j] = swapped ? prevMatches[j] : currMatches[j - 1];
                    }
                    else
                    {
                        currScore[j] = swapped ? scoreLeft : scoreUp;
                        currMatches[j] = swapped ? currMatches[j - 1] : prevMatches[j];
                    }
                }

                (prevScore, currScore) = (currScore, prevScore);
                (prevMatches, currMatches) = (currMatches, prevMatches);
            }

            int exactMatch = prevMatches[refSeqCnt - 1];

            return exactMatch / (double)Math.Min(refSeq.Length, alignSeq.Length);
        }
    }
}
