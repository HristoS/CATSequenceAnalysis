using System;
using System.Collections.Generic;
using System.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading.Tasks;

namespace SequenceAnalyses
{
    public class Generator
    {
        private static RandomNumberGenerator rngCsp = RandomNumberGenerator.Create();

        public static string GetRandomDNA(int length)
        {
            byte[] dna = new byte[length];
            rngCsp.GetBytes(dna);
            StringBuilder sb = new StringBuilder();
            foreach (byte basis in dna)
            {
                Base b = (Base)(basis % Enum.GetValues(typeof(Base)).Length);
                sb.Append(b.ToString());
            }

            return sb.ToString();
        }
    }
}