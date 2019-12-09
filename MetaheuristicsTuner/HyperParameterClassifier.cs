using System;
using ILOG.CPLEX;
using ILOG.Concert;


namespace MetaheuristicsTuner
{
    internal static class HyperParameterClassifier
    {
        /// <summary>
        /// Getting maximum congruency for naming parameter sets A or B. Linear Programming optimization model.
        /// Eq. (4) and Table II in DOI:10.1109/CEC.2019.8790261
        /// The matrices need to be manually entered here in the script (cTopA, cTopB, cBtmA, cBtmB)
        /// </summary>
        internal static void IdentifyTwoBestParameterSets(string caseSolver)
        {
            Cplex cpl = new Cplex();

            int N = 7;
            //string caseSolver = "ES"; //"SGA" (default), "ES", "PSO", "FIPS" 

            int[] cTopA;
            int[] cTopB;
            int[] cBtmA;
            int[] cBtmB;
            switch (caseSolver)
            {
                default:
                    cTopA = new int[] { 11, 1, 1, 13, 0, 0, 2 };
                    cTopB = new int[] { 0, 5, 3, 0, 5, 5, 5 };
                    cBtmA = new int[] { 2, 12, 12, 0, 14, 14, 9 };
                    cBtmB = new int[] { 4, 0, 1, 4, 0, 0, 0 };
                    break;
                case "ES":
                    cTopA = new int[] { 4, 6, 9, 3, 11, 11, 5 };
                    cTopB = new int[] { 2, 3, 1, 3, 1, 3, 3 };
                    cBtmA = new int[] { 4, 5, 2, 8, 0, 0, 6 };
                    cBtmB = new int[] { 3, 1, 3, 2, 5, 3, 2 };
                    break;
                case "PSO":
                    cTopA = new int[] { 1, 3, 2, 3, 3, 7, 5 };
                    cTopB = new int[] { 12, 11, 3, 1, 7, 8, 5 };
                    cBtmA = new int[] { 7, 5, 6, 5, 5, 1, 3 };
                    cBtmB = new int[] { 0, 1, 9, 11, 5, 4, 7 };
                    break;
                case "FIPS":
                    cTopA = new int[] { 6, 6, 7, 3, 5, 0, 8 };
                    cTopB = new int[] { 4, 6, 6, 9, 5, 9, 1 };
                    cBtmA = new int[] { 4, 4, 3, 7, 5, 10, 2 };
                    cBtmB = new int[] { 6, 4, 4, 1, 5, 1, 9 };
                    break;
            }


            INumVar[] xTopA = new INumVar[N];
            INumVar[] xBtmB = new INumVar[N];
            INumVar[] xTopB = new INumVar[N];
            INumVar[] xBtmA = new INumVar[N];
            ILinearNumExpr value = cpl.LinearNumExpr();
            for (int n = 0; n < N; n++)
            {
                xTopA[n] = cpl.BoolVar();
                xBtmB[n] = cpl.BoolVar();
                xTopB[n] = cpl.BoolVar();
                xBtmA[n] = cpl.BoolVar();
                cpl.AddEq(cpl.Sum(xTopB[n], xTopA[n]), 1);
                cpl.AddEq(cpl.Sum(xBtmA[n], xBtmB[n]), 1);
                cpl.AddEq(cpl.Sum(xTopA[n], xBtmA[n]), 1);
                cpl.AddEq(cpl.Sum(xTopB[n], xBtmB[n]), 1);
                value.AddTerm(xTopA[n], cTopA[n]);
                value.AddTerm(xTopB[n], cTopB[n]);
                value.AddTerm(xBtmA[n], cBtmA[n]);
                value.AddTerm(xBtmB[n], cBtmB[n]);
            }

            cpl.AddMaximize(value);
            cpl.Solve();

            Console.WriteLine("Parameter Grouping for Solver: {0}", caseSolver);
            for (int n = 0; n < N; n++)
            {
                Console.WriteLine("n: {0}", n);
                Console.WriteLine("xtopA: ;{0};, _____, xTopB: ;{1};", cpl.GetValue(xTopA[n]), cpl.GetValue(xTopB[n]));
                Console.WriteLine("xbtmB: ;{0};, _____, xBtmA: ;{1};", cpl.GetValue(xBtmB[n]), cpl.GetValue(xBtmA[n]));
            }
            Console.WriteLine("cost: {0}", cpl.GetObjValue());
            Console.ReadKey();




        }

    }
}
