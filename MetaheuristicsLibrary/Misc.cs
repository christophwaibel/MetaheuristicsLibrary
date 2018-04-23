using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaheuristicsLibrary.Misc
{
    /// <summary>
    /// Additional random distributions. 
    /// </summary>
    public class RandomDistributions : Random
    {

        public RandomDistributions(int rndSeed)
            : base(rndSeed)
        {        }

        /// <summary>
        /// Normal distributed random number.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="stdDev">Standard deviation of the distribution.</param>
        /// <returns>Normal distributed random number.</returns>
        public double NextGaussian(double mean, double stdDev)
        {
            //Random rand = new Random(); //reuse this if you are generating many
            double u1 = base.NextDouble(); //these are uniform(0,1) random doubles
            double u2 = base.NextDouble();
            double randStdNormal = Math.Sqrt(-2.0 * Math.Log(u1)) *
                         Math.Sin(2.0 * Math.PI * u2); //random normal(0,1)
            double randNormal =
                         mean + stdDev * randStdNormal; //random normal(mean,stdDev^2)

            return randNormal;

        }

        /// <summary>
        /// Normal distributed random number, normalized between 0 and 1. Assuming range: -5 to +5.
        /// </summary>
        /// <param name="mean">Mean of the distribution.</param>
        /// <param name="stdDev">Standard deviation of the distribution.</param>
        /// <returns>Normal distributed random number, normalized between 0 and 1.</returns>
        public double NextGaussianNorm(double mean, double stdDev)
        {
            double gauss = this.NextGaussian(mean, stdDev);
            if (gauss < -5) gauss = -5;
            else if (gauss > 5) gauss = 5;
            return (gauss + 5) / 10;
        }

    }


    /// <summary>
    /// Vector Operations.
    /// </summary>
    public static class Vector
    {
        public static double Norm(double[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; i++)
            {
                sum += Math.Pow(x[i], 2);
            }

            return Math.Sqrt(sum);

        }

        /// <summary>
        /// Computes the centroid of points
        /// </summary>
        /// <param name="X">Input points</param>
        /// <returns>Centroid</returns>
        public static double[] Centroid(double[][] X)
        {
            double[] centroid = new double[X[0].Length];
            for (int i = 0; i < X[0].Length; i++)
            {
                double sum = 0;
                for (int j = 0; j < X.Length; j++)
                {
                    sum += X[j][i];
                }
                centroid[i] = sum / X.Length;
            }
            return centroid;
        }

    }


    /// <summary>
    /// Miscellaneous functions
    /// </summary>
    public static class Misc
    {

        /// <summary>
        /// Converts decimal to binary number.
        /// </summary>
        /// <param name="n">Amount of binaries.</param>
        /// <returns>Binary encoded decimal number with n binaries.</returns>
        public static string Dec2Bin(int n)
        {
            if (n < 2) return n.ToString();

            var divisor = n / 2;
            var remainder = n % 2;

            return Dec2Bin(divisor) + remainder;
        }

        /// <summary>
        /// Convert an integer to a bitstring with specified length.
        /// Source: http://stackoverflow.com/questions/1838963/easy-and-fast-way-to-convert-an-int-to-binary
        /// </summary>
        /// <param name="value">Integer to be converted</param>
        /// <param name="len">Bitstring length</param>
        /// <returns>Bitstring</returns>
        public static string Dec2Bin(int value, int len)
        {
            return (len > 1 ? Dec2Bin(value >> 1, len - 1) : null) + "01"[value & 1];
        }

        /// <summary>
        /// Convert a bitstring to an integer.
        /// </summary>
        /// <param name="binary"></param>
        /// <returns></returns>
        public static int Bin2Dec(string binary)
        {
            int result = 0;
            foreach (char c in binary)
            {
                if (c == '1' || c == '0')
                    result = result * 2 + (c - '0');
                else
                {
                    // throw some sort of exception here
                }
            }
            return result;
        }


        /// <summary>
        /// Tells the required Bitstring length of an array of integer numbers.
        /// </summary>
        /// <param name="value">Array of double values. If they are indicated as integer, they will be round to integer.</param>
        /// <param name="intx">Indicator, wether variable is integer (true), or not.</param>
        /// <param name="lb">Array of lower bounds.</param>
        /// <param name="ub">Array of upper bounds</param>
        /// <returns>Required length of bitsting.</returns>
        public static int BinReqLength(bool[] intx, double[] lb, double[] ub)
        {
            List<int> diff = new List<int>();
            for (int i = 0; i < intx.Length; i++)
            {
                if (intx[i])
                {
                    diff.Add(Convert.ToInt32(ub[i] - lb[i]));
                }
            }
            diff.Sort(); //biggest number is last element
            int BitStringLength = Dec2Bin(diff[diff.Count - 1]).Length;

            return BitStringLength;
        }

    }


}
