using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaheuristicsLibrary.TestFunctions
{
    /// <summary>
    ///  Repository of benchmark functions for multi objective optimization.</summary>
    ///     <remarks> Types:
    ///     <para/>- MO2: 2 objectives,
    ///     <para/>- MO3: 3 obejctives,
    ///     <para/>- MOM: arbitrary number of objectives, 
    ///     <para/>- MOC: functions with constraints.
    ///     <para/>(Web-)Sources:
    ///     <para/>http://www.tik.ee.ethz.ch/sop/pisa/?page=selvar.php
    ///     <para/>https://en.wikipedia.org/wiki/Test_functions_for_optimization 
    ///     <para/>Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results.
    ///     <para/>Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.</remarks>
    public static class MO
    {

        #region 2 Objectives, no constraints
        /* 2 Objectives, NO constraints*/
        /// <summary>
        /// FONSECA-FLEMING Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [-4, 4], for i = 1, …, n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_FonsecaFleming(double[] x)
        {
            int n = x.Length;

            double sum1 = 0;
            double sum2 = 0;
            for (int i = 0; i < n; i++)
            {
                double term = x[i] - (1 / Math.Sqrt(n));
                sum1 -= Math.Pow(term, 2);
                term = x[i] + (1 / Math.Sqrt(n));
                sum2 -= Math.Pow(term, 2);
            }

            double f1 = 1 - Math.Exp(sum1);
            double f2 = 1 - Math.Exp(sum2);

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// KURSAWE Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞ (usually 3).
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> In H. P. Schwefel and R. Männer, editors, Parallel Problem Solving from Nature. 1st Workshop, PPSN I, volume 496 of Lecture Notes in Computer Science, pages 193-197, Berlin, Germany, oct 1991. Springer-Verlag. 
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [-5, 5], for i = 1, …, n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_Kursawe(double[] x)
        {
            int n = x.Length;

            double sum1 = 0;
            for (int i = 1; i <= (n - 1); i++)
            {
                sum1 -= (10 * Math.Exp(-0.2 * Math.Sqrt((Math.Pow(x[i - 1], 2) + Math.Pow(x[i], 2)))));
            }

            double sum2 = 0;
            for (int i = 1; i <= n; i++)
            {
                sum2 += (Math.Pow(Math.Abs(x[i - 1]), 0.8) + 5 * (Math.Pow(Math.Sin(x[i - 1]), 3)));
            }

            double[] fx = new double[2] { sum1, sum2 };
            return fx;
        }

        /// <summary>
        /// SCHAFFER N.1 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: 1.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// </summary>
        /// <param name="x">Decision variable vector. x ∈ [-A, A]. Values of A form 10 to 10^{5} have been used successfully. Higher values of A increase the difficulty of the problem</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_Schaffer1(double x)
        {
            double f1 = Math.Pow(x, 2);
            double f2 = Math.Pow(x - 2, 2);

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// SCHAFFER N.2 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: 1.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// </summary>
        /// <param name="x">Decision variable vector. x ∈ [-5, 10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_Schaffer2(double x)
        {
            double f1 = 0;
            if (x <= 1) f1 = -x;
            else if (1 < x && x <= 3) f1 = x - 2;
            else if (3 < x && x <= 4) f1 = 4 - x;
            else if (x > 4) f1 = x - 4;

            double f2 = Math.Pow(x - 5, 2);

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// POLONI Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: 2.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [-π, π], i = 1, 2.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_Poloni(double[] x)
        {
            double a1, a2, b1, b2;
            a1 = 0.5 * Math.Sin(1) - 2 * Math.Cos(1) + Math.Sin(2) - 1.5 * Math.Cos(2);
            a2 = 1.5 * Math.Sin(1) - Math.Cos(1) + 2 * Math.Sin(2) - 0.5 * Math.Cos(2);
            b1 = 0.5 * Math.Sin(x[0]) - 2 * Math.Cos(x[0]) + Math.Sin(x[1]) - 1.5 * Math.Cos(x[1]);
            b2 = 1.5 * Math.Sin(x[0]) - Math.Cos(x[0]) + 2 * Math.Sin(x[1]) - 0.5 * Math.Cos(x[1]);

            double f1 = 1 + Math.Pow(a1 - b1, 2) + Math.Pow(a2 - b2, 2);
            double f2 = Math.Pow(x[0] + 3, 2) + Math.Pow(x[1] + 1, 2);

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// ZITZLER-DEB-THIELE N.1 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞ (usually 30).
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [0, 1], i = 1, ..., n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_ZDT1(double[] x)
        {
            int m = x.Length;
            double f1 = x[0];

            double sum = 0;
            for (int i = 1; i < m; i++) { sum += x[i]; }
            double g = 1 + (9 / (m - 1)) * sum;

            double h = 1 - Math.Sqrt(f1 / g);

            double f2 = g * h;

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// ZITZLER-DEB-THIELE N.2 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞ (usually 30).
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [0, 1], i = 1, ..., n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_ZDT2(double[] x)
        {
            int m = x.Length;
            double f1 = x[0];

            double sum = 0;
            for (int i = 1; i < m; i++) { sum += x[i]; }
            double g = 1 + (9 / (m - 1)) * sum;

            double h = 1 - Math.Pow(f1 / g, 2);

            double f2 = g * h;

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// ZITZLER-DEB-THIELE N.3 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞ (usually 30).
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [0, 1], i = 1, ..., n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_ZDT3(double[] x)
        {
            int m = x.Length;
            double f1 = x[0];

            double sum = 0;
            for (int i = 1; i < m; i++) { sum += x[i]; }
            double g = 1 + (9 / (m - 1)) * sum;

            double h = 1 - Math.Sqrt(f1 / g) - (f1 / g) * Math.Sin(10 * Math.PI * f1);

            double f2 = g * h;

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        /// <summary>
        /// ZITZLER-DEB-THIELE N.4 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞ (usually 10).
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results.
        /// </summary>
        /// <param name="x">Decision variable vector. x1 ∈ [0, 1], xi ∈ [-5, 5], i = 2, ..., n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_ZDT4(double[] x)
        {
            int m = x.Length;
            double f1 = x[0];

            double sum = 0;
            for (int i = 1; i < m; i++) { sum += Math.Pow(x[i], 2) - 10 * Math.Cos(4 * Math.PI * x[i]); }
            double g = 1 + 10 * (m - 1) + sum;

            double h = 1 - Math.Sqrt(f1 / g);

            double f2 = g * h;

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        // ZITZLER-DEB-THIELE N.5
        //    # of real variables = 0
        //    # of bin variables = 11
        //    # of bits for binvar1 = 30
        //    # of bits for binvar2-11 = 5
        //    # of objectives = 2
        //    # of constraints = 0
        //
        //  Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results
        /*
        #ifdef zdt5
        void test_problem (double *xreal, double *xbin, int **gene, double *obj, double *constr)
        {
            int i, j;
            int u[11];
            int v[11];
            double f1, f2, g, h;
            for (i=0; i<11; i++)
            {
                u[i] = 0;
            }
            for (j=0; j<30; j++)
            {
                if (gene[0][j] == 1)
                {
                    u[0]++;
                }
            }
            for (i=1; i<11; i++)
            {
                for (j=0; j<4; j++)
                {
                    if (gene[i][j] == 1)
                    {
                        u[i]++;
                    }
                }
            }
            f1 = 1.0 + u[0];
            for (i=1; i<11; i++)
            {
                if (u[i] < 5)
                {
                    v[i] = 2 + u[i];
                }
                else
                {
                    v[i] = 1;
                }
            }
            g = 0;
            for (i=1; i<11; i++)
            {
                g += v[i];
            }
            h = 1.0/f1;
            f2 = g*h;
            obj[0] = f1;
            obj[1] = f2;
            return;
        }
        #endif*/

        /// <summary>
        /// ZITZLER-DEB-THIELE N.6 Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: ∞ (usually 10).
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Zizler, Deb, Thiele (2000). Comparison of Multiobjective Evolutionary Algorithms: Empirical Results.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [0, 1], i = 1, ..., n.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO2_ZDT6(double[] x)
        {
            int m = x.Length;
            double f1 = 1 - Math.Exp(-4 * x[0]) * Math.Pow(Math.Sin(6 * Math.PI * x[0]), 6);

            double sum = 0;
            for (int i = 1; i < m; i++)
            {
                sum += x[i];
            }
            double g = 1 + 9 * Math.Pow(sum / (m - 1), 0.25);

            double h = 1 - Math.Pow(f1 / g, 2);

            double f2 = g * h;

            if (double.IsNaN(f2)) f2 = double.MaxValue;
            if (double.IsNaN(f1)) f1 = double.MaxValue;

            double[] fx = new double[2] { f1, f2 };
            return fx;
        }

        #endregion





        #region 3 Objectives, no constraints
        /* 3 Objectives, NO constraints */

        /// <summary>
        /// VIENNET Function. 
        /// <para/> Objectives: 3.
        /// <para/> Real variables: 2.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// </summary>
        /// <param name="x">Decision variable vector. xi ∈ [-3, 3], i = 1, 2.</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MO3_Viennet(double[] x)
        {
            double f1 = 0.5 * (Math.Pow(x[0], 2) + Math.Pow(x[1], 2)) + Math.Sin(Math.Pow(x[0], 2) + Math.Pow(x[1], 2));
            double f2 = (Math.Pow(3 * x[0] - 2 * x[1] + 4, 2) / 8) + (Math.Pow(x[0] - x[1] + 1, 2) / 27) + 15;
            double f3 = (1 / (Math.Pow(x[0], 2) + Math.Pow(x[1], 2) + 1)) - 1.1 * Math.Exp(-(Math.Pow(x[0], 2) + Math.Pow(x[1], 2)));

            double[] fx = new double[3] { f1, f2, f3 };
            return fx;
        }

        #endregion





        #region Many Objectives, no constraints
        /* Many Objectives, NO constraints*/
        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.1 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 5 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ1(double[] x, int M)
        {
            int i = 0, j = 0;
            int n = x.Length;   // length of decision variables
            int k = n - M + 1;  // M is number of objectives

            double g = 0;
            for (i = n - k + 1; i <= n; i++)
            {
                g += Math.Pow(x[i - 1] - 0.5, 2) - Math.Cos(20 * Math.PI * (x[i - 1] - 0.5));
            }
            g = 100 * (k + g);

            double[] fx = new double[M];
            for (i = 1; i <= M; i++)
            {
                double f = 0.5 * (1 + g);
                for (j = M - i; j >= 1; j--)
                {
                    f *= x[j - 1];
                }
                if (i > 1)
                {
                    f *= 1 - x[(M - i + 1) - 1];
                }

                fx[i - 1] = f;
            }

            return fx;
        }

        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.2 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 10 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ2(double[] x, int M)
        {
            int i = 0, j = 0;
            int n = x.Length;
            int k = n - M + 1;

            double g = 0;
            for (i = n - k + 1; i <= n; i++)
            {
                g += Math.Pow(x[i - 1] - 0.5, 2);
            }

            double[] fx = new double[M];
            for (i = 1; i <= M; i++)
            {
                double f = (1 + g);
                for (j = M - i; j >= 1; j--)
                {
                    f *= Math.Cos(x[j - 1] * Math.PI / 2);
                }
                if (i > 1)
                {
                    f *= Math.Sin(x[(M - i + 1) - 1] * Math.PI / 2);
                }

                fx[i - 1] = f;
            }

            return fx;
        }

        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.3 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 10 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ3(double[] x, int M)
        {
            int i = 0, j = 0;
            int n = x.Length;
            int k = n - M + 1;

            double g = 0;
            for (i = n - k + 1; i <= n; i++)
            {
                g += Math.Pow(x[i - 1] - 0.5, 2) - Math.Cos(20.0 * Math.PI * (x[i - 1] - 0.5));
            }
            g = 100 * (k + g);

            double[] fx = new double[M];
            for (i = 1; i <= M; i++)
            {
                double f = (1 + g);
                for (j = M - i; j >= 1; j--)
                {
                    f *= Math.Cos(x[j - 1] * Math.PI / 2.0);
                }
                if (i > 1)
                {
                    f *= Math.Sin(x[(M - i + 1) - 1] * Math.PI / 2.0);
                }

                fx[i - 1] = f;
            }

            return fx;
        }

        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.4 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 10 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ4(double[] x, int M, double alpha = 100)
        {
            int i = 0, j = 0;
            int n = x.Length;
            int k = n - M + 1;

            double g = 0;
            for (i = n - k + 1; i <= n; i++)
            {
                g += Math.Pow(x[i - 1] - 0.5, 2);
            }

            double[] fx = new double[M];
            for (i = 1; i <= M; i++)
            {
                double f = (1 + g);
                for (j = M - i; j >= 1; j--)
                {
                    f *= Math.Cos(Math.Pow(x[j - 1], alpha) * Math.PI / 2);
                }
                if (i > 1)
                {
                    f *= Math.Sin(Math.Pow(x[(M - i + 1) - 1], alpha) * Math.PI / 2);
                }

                fx[i - 1] = f;
            }

            return fx;
        }

        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.5 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 10 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ5(double[] x, int M)
        {
            int i = 0, j = 0;
            int n = x.Length;
            int k = n - M + 1;
            double[] theta = new double[M];
            double t = 0;
            double g = 0;

            for (i = n - k + 1; i <= n; i++)
            {
                g += Math.Pow(x[i - 1] - 0.5, 2);
            }

            t = Math.PI / (4 * (1 + g));
            theta[0] = x[0] * Math.PI / 2;
            for (i = 2; i <= M - 1; i++)
            {
                theta[i - 1] = t * (1 + 2 * g * x[i - 1]);
            }

            double[] fx = new double[M];
            for (i = 1; i <= M; i++)
            {
                double f = (1 + g);
                for (j = M - i; j >= 1; j--)
                {
                    f *= Math.Cos(theta[j - 1]);
                }
                if (i > 1)
                {
                    f *= Math.Sin(theta[(M - i + 1) - 1]);
                }

                fx[i - 1] = f;
            }

            return fx;
        }

        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.6 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 10 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ6(double[] x, int M)
        {
            int i = 0, j = 0;
            int n = x.Length;
            int k = n - M + 1;
            double[] theta = new double[M];
            double t = 0;
            double g = 0;

            for (i = n - k + 1; i <= n; i++)
            {
                g += Math.Pow(x[i - 1], 0.1);
            }

            t = Math.PI / (4 * (1 + g));
            theta[0] = x[0] * Math.PI / 2;
            for (i = 2; i <= M - 1; i++)
            {
                theta[i - 1] = t * (1 + 2 * g * x[i - 1]);
            }

            double[] fx = new double[M];
            for (i = 1; i <= M; i++)
            {
                double f = (1 + g);
                for (j = M - i; j >= 1; j--)
                {
                    f *= Math.Cos(theta[j - 1]);
                }
                if (i > 1)
                {
                    f *= Math.Sin(theta[(M - i + 1) - 1]);
                }

                fx[i - 1] = f;
            }

            return fx;
        }

        /// <summary>
        /// DEB-THIELE-LAUMANNS-ZITZLER N.7 Function. 
        /// <para/> Objectives: M ϵ [2,10].
        /// <para/> Real variables: ∞. k = 20 suggested. (k = n - M + 1), where n is number of variables.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Source: Deb, Thiele, Laumanns, Zitzler (2001). Scalable Test Problems for Evolutionary Multi-Objective Optimization.
        /// </summary>
        /// <param name="x">Decision variable vector. xi, i = 1, ..., n.</param>
        /// <param name="M">Number of Objectives. M ϵ [2,10].</param>
        /// <returns>Objective function values f(x).</returns>
        public static double[] MOM_DTLZ7(double[] x, int M)
        {
            int i = 0, j = 0;
            int n = x.Length;
            int k = n - M + 1;
            double g = 0;
            double h = 0;

            for (i = n - k + 1; i <= n; i++)
            {
                g += x[i - 1];
            }
            g = 1 + 9 * g / k;

            double[] fx = new double[M];
            for (i = 1; i <= M - 1; i++)
            {
                fx[i - 1] = x[i - 1];
            }

            for (j = 1; j <= M - 1; j++)
            {
                h += x[j - 1] / (1 + g) * (1 + Math.Sin(3 * Math.PI * x[j - 1]));
            }
            h = M - h;
            fx[M - 1] = (1 + g) * h;

            return fx;
        }

        #endregion




        #region With Constrainst
        /* WITH CONSTRAINTS 
        negative value means constraint is violated.
        0 is ok, if <= or >=
        if only < or >, then 0 is also indicating violation*/

        /// <summary>
        /// BINH-KORN Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: 2.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 2.
        /// </summary>
        /// <param name="x">Decision variable vector. x1 ϵ [0, 5], x2 ϵ [0, 3].</param>
        /// <returns>[0][]: Objective function values f(x).
        /// <para/>[1][]: Constraints violation. Negative value means constraint is violated. 0-value accepted, if smaller/greater-equals constraint, otherwise 0 also indicates violation.</returns>
        public static double[][] MO2C_BinhKorn(double[] x)
        {
            double x1 = x[0];
            double x2 = x[1];

            double f1 = 4 * Math.Pow(x1, 2) + 4 * Math.Pow(x2, 2);
            double f2 = Math.Pow(x1 - 5, 2) + Math.Pow(x2 - 5, 2);

            double g1 = 1 - (Math.Pow(x1 - 5, 2) + Math.Pow(x2 - 5, 2)) / 25;  // <=25 
            double g2 = (Math.Pow(x1 - 8, 2) + Math.Pow(x2 + 3, 2)) / 7.7 - 1; // >= 7.7

            double[] fx = new double[2] { f1, f2 };
            double[] g = new double[2] { g1, g2 };
            double[][] result = new double[2][];
            result[0] = fx;
            result[1] = g;
            return result;
        }

        /// <summary>
        /// OSYCZKA-KUNDU Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: 6.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 6.
        /// </summary>
        /// <param name="x">Decision variable vector. x1, x2, x6 ϵ [0, 10], x3, x5 ϵ [1, 5], x4 ϵ [0, 6].</param>
        /// <returns>[0][]: Objective function values f(x).
        /// <para/>[1][]: Constraints violation. Negative value means constraint is violated. 0-value accepted, if smaller/greater-equals constraint, otherwise 0 also indicates violation.</returns>
        public static double[][] MO2C_OsyczkaKundu(double[] x)
        {
            double f1 = -(25 * Math.Pow(x[0] - 2, 2) + Math.Pow(x[1] - 2, 2) + Math.Pow(x[2] - 1, 2) + Math.Pow(x[3] - 4, 2) + Math.Pow(x[4] - 1, 2));
            double f2 = 0;
            for (int i = 0; i < 6; i++) { f2 += Math.Pow(x[i], 2); }

            double[] g = new double[6];
            g[0] = (x[0] + x[1]) / 2 - 1;                       // >=2
            g[1] = 1 - (x[0] - x[1]) / 6;                       // <=6
            g[2] = 1 - x[1] / 2 + x[0] / 2;                     // <=2
            g[3] = 1 - x[0] / 2 + 3 * x[1] / 2;                 // <= 2
            g[4] = 1 - Math.Pow(x[2] - 3, 2) / 4 - x[3] / 4;    // <= 4
            g[5] = Math.Pow(x[4] - 3, 2) / 4 + x[5] / 4 - 1;    // >= 4

            double[] fx = new double[2] { f1, f2 };
            double[][] result = new double[2][];
            result[0] = fx;
            result[1] = g;
            return result;
        }

        /// <summary>
        /// CHAKONG-HAIMES Function. 
        /// <para/> Objectives: 2.
        /// <para/> Real variables: 2.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 2.
        /// </summary>
        /// <param name="x">Decision variable vector. x1, x2 ϵ [-20, 20].</param>
        /// <returns>[0][]: Objective function values f(x).
        /// <para/>[1][]: Constraints violation. Negative value means constraint is violated. 0-value accepted, if smaller/greater-equals constraint, otherwise 0 also indicates violation.</returns>
        public static double[][] MO2C_ChakongHaimes(double[] x)
        {
            double f1 = 2 + Math.Pow(x[0] - 2, 2) + Math.Pow(x[1] - 1, 2);
            double f2 = 9 * x[0] - Math.Pow(x[1] - 1, 2);
            double g1 = 1 - (Math.Pow(x[0], 2) + Math.Pow(x[1], 2)) / 225;
            double g2 = 3 * x[1] / 10 - x[0] / 10 - 1;

            double[][] result = new double[2][];
            result[0] = new double[] { f1, f2 };
            result[1] = new double[] { g1, g2 };
            return result;
        }

        #endregion
    }
}
