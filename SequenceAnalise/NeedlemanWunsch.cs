using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

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
            int refSeqCnt = refSeq.Length + 1;
            int alineSeqCnt = alignSeq.Length + 1;

            int[,] scoringMatrix = new int[alineSeqCnt, refSeqCnt];

            //Initialization Step - filled with 0 for the first row and the first column of matrix
            for (int i = 0; i < alineSeqCnt; i++)
            {
                scoringMatrix[i, 0] = i * gapScore;
            }

            for (int j = 0; j < refSeqCnt; j++)
            {
                scoringMatrix[0, j] = j * gapScore;
            }

            //Matrix Fill Step
            for (int i = 1; i < alineSeqCnt; i++)
            {
                for (int j = 1; j < refSeqCnt; j++)
                {
                    int scroeDiag = 0;
                    if (refSeq.Substring(j - 1, 1) == alignSeq.Substring(i - 1, 1))
                        scroeDiag = scoringMatrix[i - 1, j - 1] + matchScore;
                    else
                        scroeDiag = scoringMatrix[i - 1, j - 1] + mismatchScore;

                    int scroeLeft = scoringMatrix[i, j - 1] + gapScore;
                    int scroeUp = scoringMatrix[i - 1, j] + gapScore;

                    int maxScore = Math.Max(Math.Max(scroeDiag, scroeLeft), scroeUp);

                    scoringMatrix[i, j] = maxScore;
                }
            }

            //Traceback Step
            char[] alineSeqArray = alignSeq.ToCharArray();
            char[] refSeqArray = refSeq.ToCharArray();

            string AlignmentA = string.Empty;
            string AlignmentB = string.Empty;
            int m = alineSeqCnt - 1;
            int n = refSeqCnt - 1;
            while (m > 0 && n > 0)
            {
                int scroeDiag = 0;

                //Remembering that the scoring scheme is +2 for a match, -1 for a mismatch and -2 for a gap
                if (alineSeqArray[m - 1] == refSeqArray[n - 1])
                    scroeDiag = 2;
                else
                    scroeDiag = -1;

                if (m > 0 && n > 0 && scoringMatrix[m, n] == scoringMatrix[m - 1, n - 1] + scroeDiag)
                {
                    AlignmentA = refSeqArray[n - 1] + AlignmentA;
                    AlignmentB = alineSeqArray[m - 1] + AlignmentB;
                    m = m - 1;
                    n = n - 1;
                }
                else if (n > 0 && scoringMatrix[m, n] == scoringMatrix[m, n - 1] - 2)
                {
                    AlignmentA = refSeqArray[n - 1] + AlignmentA;
                    AlignmentB = "-" + AlignmentB;
                    n = n - 1;
                }
                else //if (m > 0 && scoringMatrix[m, n] == scoringMatrix[m - 1, n] + -2)
                {
                    AlignmentA = "-" + AlignmentA;
                    AlignmentB = alineSeqArray[m - 1] + AlignmentB;
                    m = m - 1;
                }
            }

            int exactMatch = 0;
            int gaps = 0;
            for (int i = 0; i < AlignmentA.Length; i++)
            {
                exactMatch += AlignmentA[i] == AlignmentB[i] ? 1 : 0;
                gaps += AlignmentA[i] == '-' || AlignmentB[i] == '-' ? 1 : 0;
            }

            //return (exactMatch, gaps);
            return exactMatch / (double)Math.Min(refSeq.Length, alignSeq.Length);
            //int[,] matrix = InitScoringMatrix(alignSeq.Length, refSeq.Length);

            //FillScoringMatrix(refSeq, alignSeq, matrix);
            //print(refSeq, alignSeq, matrix);
            //string copyRefSeq = refSeq;
            //string copyAlignSeq = alignSeq;
            //return Traceback(ref copyRefSeq, ref copyAlignSeq, matrix);
        }
    }
}