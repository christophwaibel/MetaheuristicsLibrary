using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MetaheuristicsLibrary;
using MetaheuristicsLibrary.TestFunctions;
using MetaheuristicsLibrary.SolversMO;
using MetaheuristicsLibrary.SolversSO;

namespace Tester
{
    class Program
    {

        /// <summary>
        /// Test function
        /// </summary>
        /// <param name="x"></param>
        /// <returns></returns>
        static double[] testfnc(double[] x)
        {
            Func<double[], int, double[]> FncIn = MO.MOM_DTLZ3;
            int m = 3;
            double[] val;
            //double alpha = 100;
            val = FncIn(x, m);
            return val;
        }


        static void mMain(string[] args)
        {
            int nVar = 5;                                //number of decision variables                      
            int mObj = 3;                                //number of objectives
            double[] lb = new double[nVar];          //lower bound values for each variable
            double[] ub = new double[nVar];          //upper bound values for each variable
            for (int i = 0; i < nVar; i++)              //here I'm assigning lower bound 0 and upper bound 1 to every variable
            {
                lb[i] = 0;
                ub[i] = 10;
            }
            bool[] intX = new bool[nVar];
            intX[0] = true;
            intX[1] = true;
            intX[2] = true;
            intX[3] = true;
            intX[4] = true;
            //intX[5] = false;
            //intX[6] = false;
            //intX[7] = false;
            //intX[8] = false;
            //intX[9] = false;
            //intX[10] = false;
            //intX[11] = false;

            Func<double[], double[]> testfunc = testfnc;                                            //evaluation function
            int randomSeed = ((int)DateTime.Now.Ticks & 0x0000FFFF);                               // random seed for random number generator


            //SPEA2
            List<double[]> x0pop = new List<double[]>();
            //x0pop.Add(new double[12] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });        //initial solution 1
            //x0pop.Add(new double[12] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });        //initial solution 2
            x0pop.Add(new double[5] { 0, 0, 0, 0, 0 });        //initial solution 1
            x0pop.Add(new double[5] { 1, 1, 1, 1, 1 });        //initial solution 2

            SPEA2 spea2 = new SPEA2(mObj, nVar, lb, ub, testfunc, randomSeed, x0pop, intX);

            //change solver parameters
            spea2.ga_genmax = 300;                        // max amount of generations (which are the solver iterations)
            spea2.ga_nPop = 20;                           //population size: number of solution candidates per iteration
            spea2.ga_PCrossover = 1.0;
            spea2.ga_PMutation = 0.0;


            //initialize solver
            Console.WriteLine("Initializing Algorithm:... {0}", spea2.GetSolverName());
            Console.WriteLine("On Problem:... {0}", "MO.MOM_DTLZ4");
            spea2.initialize();
            Console.ReadKey();
            //foreach (double[] current_best_objvalues in spea2.objvalsBestPop_MO)
            //{
            //    Console.WriteLine("{0}, {1}, {2}", current_best_objvalues[0], current_best_objvalues[1], current_best_objvalues[2]);
            //}

            ////show decision variables for each solution
            //int j = 0;
            //foreach (double[] current_best_x in spea2.xPopulationArchive)
            //{
            //    Console.WriteLine(" ");
            //    Console.Write("individual {0}, x:  ", j);
            //    for (int n = 0; n < spea2.nVar; n++)
            //    {
            //        Console.Write("{0}, ", current_best_x[n]);
            //    }
            //    j++;
            //}



            //for each generation
            //for (int t = 0; t < 100; t++)
            //{
            //solve one generation
            Console.WriteLine("***************************");
            //Console.WriteLine("Generation: {0}", t);
            spea2.Solve();

            //looping through pareto front of the current generation
            //show three objective function values for each solution
            foreach (double[] current_best_objvalues in spea2.objvalsBestPop_MO)
            {
                Console.WriteLine("{0}, {1}, {2}", current_best_objvalues[0], current_best_objvalues[1], current_best_objvalues[2]);
            }

            //show decision variables for each solution
            int j = 0;
            foreach (double[] current_best_x in spea2.xPopulationArchive)
            {
                Console.WriteLine(" ");
                Console.Write("individual {0}, x:  ", j);
                for (int n = 0; n < spea2.nVar; n++)
                {
                    Console.Write("{0}, ", current_best_x[n]);
                }
                j++;
            }

            //Console.WriteLine(" ");
            //Console.ReadKey();
            //}
            Console.WriteLine(" ");
            Console.ReadKey();

        }



        static void emilieMain(string[] args)
        {
            Random rnd = new Random();

            List<double[]> x = new List<double[]>();            //these are your solutions
            List<double[]> fx = new List<double[]>();           //these are the respective objective values

            for (int i = 0; i < 100; i++)
            {
                double[] xi = new double[10];
                for (int u = 0; u < xi.Length; u++)
                {
                    xi[u] = rnd.NextDouble();
                }
                x.Add(xi);

                double[] fxi = new double[2];
                for (int u = 0; u < fxi.Length; u++)
                {
                    fxi[u] = rnd.NextDouble() * 10000;
                }
                fx.Add(fxi);
            }

            List<List<int>> Fronts = new List<List<int>>();         //fronts with indices of solutions per front. refers to unsorted x.
            int?[] Rank = new int?[x.Count];                        //Ranks of sorted x. Rank = 0 is Pareto Front.
            ParetoFront.S_MO_NonDominatedFront(ref x, ref fx, ref Rank, ref Fronts);
            for (int i = 0; i < Fronts[0].Count; i++)
            {
                //Console.WriteLine("fx1: {0}, fx2: {1}, fx3: {2}", fx[i][0], fx[i][1], fx[i][2]);      //3 objectives
                Console.WriteLine("fx1: {0}, fx2: {1}", fx[i][0], fx[i][1]);                            //2 objectives
            }


        }



        /// <summary>
        /// Testing a SO solver from MetaheuristicsLibrary.SovlersSO
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
        {
            int dvar = 2;
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            bool[] xint = new bool[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 100;
                xint[i] = false;
            }



            //Hillclimber hc = new Hillclimber(lb, ub, 100, Testfunctions.L_Ackley, 1, 0.1);
            //hc.solve();
            //Console.WriteLine(hc.get_fxoptimum());
            //Console.ReadLine();


            Dictionary<string, object> settings = new Dictionary<string, object>();
            settings.Add("maxgen", 100);
            settings.Add("popsize", 50);
            settings.Add("pcross", 0.7);
            settings.Add("pmut", 0.3);
            settings.Add("d", 0.1);
            settings.Add("r", 0.1);
            settings.Add("k", 6);
            //double[][] x0 = new double[50][];
            //double[] fx0 = new double[50];
            //for (int i = 0; i < x0.Length; i++)
            //{
            //    x0[i] = new double[dvar];
            //    for (int u = 0; u < x0[i].Length; u++)
            //    {
            //        x0[i][u] = i * u;
            //    }
            //    fx0[i] = i;
            //}
            //SimpleGA ga = new SimpleGA(lb, ub, xint, 100, Testfunctions.L_Ackley, 1, settings, x0, fx0);
            SimpleGA ga = new SimpleGA(lb, ub, xint, 1000, Testfunctions.L_Ackley, 1, settings);
            //SimpleGA ga = new SimpleGA(lb, ub, xint, 100, Testfunctions.stupid, 1, settings);
            ga.solve();
            Console.WriteLine(ga.get_fxoptimum().ToString());
            Console.ReadLine();
        }

        /// <summary>
        /// example for FrOG Single Objective. Not working? Check project: FrOG.vs
        /// </summary>
        /// <param name="args"></param>
        static void mmMain(string[] args)
        {
            int dvar = 5;
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5;
                ub[i] = 5;
            }

            Func<List<decimal>, double> evalFunc = x =>
            {
                double[] dbl = new double[dvar];
                for (int i = 0; i < dvar; i++)
                {
                    dbl[i] = Convert.ToDouble(x[i]);
                }
                return Testfunctions.L_Ackley(dbl);
            };


            List<Variable> thomasVariable = new List<Variable>();
            for (int i = 0; i < dvar; i++)
            {
                Variable var = new Variable(Convert.ToDecimal(lb[i]), Convert.ToDecimal(ub[i]), false);
                thomasVariable.Add(var);
            }

            HillclimberFROG hcFROG = new HillclimberFROG();
            bool run = hcFROG.RunSolver(thomasVariable, evalFunc, "don't know", "don't know");
            Console.WriteLine(hcFROG.GetErrorMessage());
            if (run)
            {
                Console.WriteLine("cost: {0}", hcFROG.get_fxoptimum());
            }
            Console.ReadKey();

        }


    }

    public class Testfunctions
    {
        /// <summary>
        /// ACKLEY Function. 
        /// <para/> Objectives: 1.
        /// <para/> Real variables: ∞.
        /// <para/> Binary variables: 0.
        /// <para/> Constraints: 0.
        /// <para/> Type: Many local minima.
        /// <para/> Description: The Ackley function is widely used for testing optimization algorithms. In its two-dimensional form, it is characterized by a nearly flat outer region, and a large hole at the centre. The function poses a risk for optimization algorithms, particularly hillclimbing algorithms, to be trapped in one of its many local minima.
        /// </summary>
        /// <param name="x">Decision variable vector. Usually: xi ∈ [-32.768, 32.768], or [-5, 5], for all i = 1, …, d.</param>
        /// <param name="a">(optional) parameter. Default 20.</param>
        /// <param name="b">(optional) parameter. Default 0.2.</param>
        /// <param name="c">(optional) parameter. Default 2PI.</param>
        /// <returns>Objective function value f(x). Global minimum f(x*) = 0, at x* = (0, ..., 0).</returns>
        public static double L_Ackley(double[] x)
        {
            double a = 20;
            double b = 0.2;
            double c = (2 * Math.PI);
            double term1, term2, xi, sum1 = 0, sum2 = 0;
            int d = x.Length;

            for (int i = 0; i < d; i++)
            {
                xi = x[i];
                sum1 += (Math.Pow(xi, 2));
                sum2 += (Math.Cos(c * xi));
            }

            term1 = -a * Math.Exp(-b * Math.Sqrt(sum1 / d));
            term2 = -Math.Exp(sum2 / d);

            return (term1 + term2 + a + Math.Exp(1));
        }


        public static double stupid(double[] x)
        {
            double sum = 0;
            for (int i = 0; i < x.Length; i++)
            {
                sum += x[i];
            }
            return sum;
        }
    }
}
