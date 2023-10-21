using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace SequenceAnalise
{
    public static class PatternSearchingKMP
    {
        public static List<int> KMPSearch(string pat, string txt)
        {
            List<int> result = new List<int>();

            int M = pat.Length;
            int N = txt.Length;

            // Create lps[] that will hold the longest
            // prefix suffix values for pattern
            int[] lps = new int[M];

            // Index for pat[]
            int j = 0;

            // Preprocess the pattern (calculate lps[]
            // array)
            computeLPSArray(pat, M, lps);

            int i = 0;
            while ((N - i) >= (M - j))
            {
                if (pat[j] == txt[i])
                {
                    j++;
                    i++;
                }
                if (j == M)
                {
                    //Console.Write("Found pattern "
                    //              + "at index " + (i - j));
                    result.Add(i - j);
                    j = lps[j - 1];
                }

                // Mismatch after j matches
                else if (i < N && pat[j] != txt[i])
                {
                    // Do not match lps[0..lps[j-1]] characters,
                    // they will match anyway
                    if (j != 0)
                        j = lps[j - 1];
                    else
                        i = i + 1;
                }
            }

            return result;
        }

        private static void computeLPSArray(string pat, int M, int[] lps)
        {
            // Length of the previous longest prefix suffix
            int len = 0;
            int i = 1;
            lps[0] = 0;

            // The loop calculates lps[i] for i = 1 to M-1
            while (i < M)
            {
                if (pat[i] == pat[len])
                {
                    len++;
                    lps[i] = len;
                    i++;
                }
                else // (pat[i] != pat[len])
                {
                    // This is tricky. Consider the example.
                    // AAACAAAA and i = 7. The idea is similar
                    // to search step.
                    if (len != 0)
                    {
                        len = lps[len - 1];

                        // Also, note that we do not increment
                        // i here
                    }
                    else // len = 0
                    {
                        lps[i] = len;
                        i++;
                    }
                }
            }
        }
    }
}