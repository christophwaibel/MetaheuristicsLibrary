using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MetaheuristicsLibrary;
using MetaheuristicsTuner.Testfunctions;
using MetaheuristicsLibrary.SolversMO;
using MetaheuristicsLibrary.SolversSO;
using MetaheuristicsLibrary.Misc;
using System.IO;

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

        /// <summary>
        /// multi objective SPEA-2 for Emilie
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
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
            //intX[0] = true;
            //intX[1] = true;
            //intX[2] = true;
            //intX[3] = true;
            //intX[4] = true;
            //intX[5] = false;
            //intX[6] = false;
            //intX[7] = false;
            //intX[8] = false;
            //intX[9] = false;
            //intX[10] = false;
            //intX[11] = false;

            Func<double[], double[]> testfunc = testfnc;                                            //evaluation function
            //int randomSeed = ((int)DateTime.Now.Ticks & 0x0000FFFF);                               // random seed for random number generator
            int randomSeed = 1;

            //SPEA2
            List<double[]> x0pop = new List<double[]>();
            //x0pop.Add(new double[12] { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 });        //initial solution 1
            //x0pop.Add(new double[12] { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 });        //initial solution 2
            //x0pop.Add(new double[5] { 0, 0, 0, 0, 0 });        //initial solution 1
            //x0pop.Add(new double[5] { 1, 1, 1, 1, 1 });        //initial solution 2

            SPEA2 spea2 = new SPEA2(mObj, nVar, lb, ub, testfunc, randomSeed, x0pop, intX);

            //change solver parameters
            spea2.ga_genmax = 3000;                        // max amount of generations (which are the solver iterations)
            spea2.ga_nPop = 20;                           //population size: number of solution candidates per iteration
            spea2.ga_PCrossover = 0.8;
            spea2.ga_PMutation = 0.2;


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
            for (int t = 0; t < 100; t++)
            {
            //solve one generation
            Console.WriteLine("***************************");
            Console.WriteLine("Generation: {0}", t);
            spea2.Solve(1);

            //looping through pareto front of the current generation
            //show three objective function values for each solution
            foreach (double[] current_best_objvalues in spea2.objvalsBestPop_MO)
            {
                Console.WriteLine("{0}, {1}, {2}", current_best_objvalues[0], current_best_objvalues[1], current_best_objvalues[2]);
            }

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

            Console.WriteLine(" ");
            Console.ReadKey();
            }
            Console.WriteLine(" ");
            Console.ReadKey();

        }


        /// <summary>
        /// pareto front example
        /// </summary>
        /// <param name="args"></param>
        static void pfMain(string[] args)
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

        static void rndMain(string[] args)
        {
            List<string> log = new List<string>();
            RandomDistributions rnd = new RandomDistributions(0);
            for (int i = 0; i < 5000; i++)
            {
                log.Add(Convert.ToString(rnd.NextGaussian(0, 0.3)));
            }
            Console.ReadKey();


            string fileName = @"C:\_CHRIS\rnd.txt";
            using (FileStream fs = new FileStream(fileName, FileMode.Append, FileAccess.Write))
            using (StreamWriter sw = new StreamWriter(fs))
            {
                for (int i = 0; i < log.Count; i++)
                {
                    sw.WriteLine(log[i]);
                }
            }

        }

        /// <summary>
        /// Testing a SO solver from MetaheuristicsLibrary.SovlersSO
        /// </summary>
        /// <param name="args"></param>
        static void sosoMain(string[] args)
        {
           
            int seeds = 20;
            double[] optis = new double[seeds];
            int dvar = 13;
            int evalcount = (dvar+1) * 100;
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            bool[] xint = new bool[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -30;
                //lb[i] = -3;
                ub[i] = 1;
                xint[i] = true;
            }

            Func<double[], double> testfunc = SO.V_Rosenbrock;


            //Hillclimber hc = new Hillclimber(lb, ub, xint, 100, testfunc, 1, 0.1);
            //hc.solve();
            //Console.WriteLine(hc.get_fxoptimum());
            //Console.ReadLine();

            //simpleGAsettings.Add("popsize", simplegapop);            //x[0]  50
            //simpleGAsettings.Add("maxgen", Convert.ToInt32(Math.Floor(_maxfuncs / Convert.ToDouble(simplegapop))));            //no x. its a function of popsize and max evals.  50
            //simpleGAsettings.Add("k", hpo_x[1]);                //x[1]  6
            //simpleGAsettings.Add("pcross", hpo_x[2]);           //x[2]  0.7
            //simpleGAsettings.Add("pmut", hpo_x[3]);             //x[3]  0.3
            //simpleGAsettings.Add("d", hpo_x[4]);                //x[4]  0.1
            //simpleGAsettings.Add("r", hpo_x[5]);                //x[5]  0.1
            //int SGAelite = Convert.ToInt16(hpo_x[6] * (int)Math.Round((double)simplegapop / 2, 0)); //from percentage to integer
            //simpleGAsettings.Add("elite", SGAelite);                //x[6]  0 - (popsize/2). 



            //Dictionary<string, object> settings = new Dictionary<string, object>();
            //int simplegapop = Convert.ToInt32(14);   //x[0]
            //simplegapop += simplegapop % 2;
            //settings.Add("maxgen", Convert.ToInt32(Math.Floor(evalcount / Convert.ToDouble(simplegapop))));
            //settings.Add("popsize", simplegapop);
            //settings.Add("k", 11);               //x[1]  6
            //settings.Add("pcross", 1);          //x[2]  0.7
            //settings.Add("pmut", 0.2);          //x[3]  0.3
            //settings.Add("d", 0.01);             //x[4]  0.1
            //settings.Add("r", 0.2);             //x[5]  0.1
            ////int SGAelite = Convert.ToInt16(0 * (int)Math.Round((double)simplegapop / 2, 0)); // x[6].  from percentage to integer
            //int SGAelite = 1;
            //settings.Add("elite", SGAelite);
            //SimpleGA[] ga = new SimpleGA[seeds];
            //for (int i = 0; i < seeds; i++)
            //{
            //    ga[i] = new SimpleGA(lb, ub, xint, evalcount, testfunc, i, settings);
            //    ga[i].solve();
            //    optis[i] = ga[i].get_fxoptimum();
            //}
            //Console.WriteLine("ga average: {0}", optis.Average());
            //Console.ReadLine();




            //Dictionary<string, object> settingsB = new Dictionary<string, object>();
            //int simplegapopB = Convert.ToInt32(6);   //x[0]
            //simplegapopB += simplegapopB % 2;
            //settingsB.Add("maxgen", Convert.ToInt32(Math.Floor(evalcount / Convert.ToDouble(simplegapopB))));
            //settingsB.Add("popsize", simplegapopB);
            //settingsB.Add("k", 36.7352);               //x[1]  6
            //settingsB.Add("pcross", 0.87624);          //x[2]  0.7
            //settingsB.Add("pmut", 0.8189);          //x[3]  0.3
            //settingsB.Add("d", 0.63072);             //x[4]  0.1
            //settingsB.Add("r", 1.75361);             //x[5]  0.1
            //int SGAeliteB = Convert.ToInt16(0.43565 * (int)Math.Round((double)simplegapopB / 2, 0)); // x[6].  from percentage to integer
            //settingsB.Add("elite", SGAeliteB);
            //SimpleGA[] ga_B = new SimpleGA[seeds];
            //for (int i = 0; i < seeds; i++)
            //{
            //    ga_B[i] = new SimpleGA(lb, ub, xint, evalcount, testfunc, i, settingsB);
            //    ga_B[i].solve();
            //    optis[i] = ga_B[i].get_fxoptimum();
            //}
            //Console.WriteLine("ga B average: {0}", optis.Average());
            //Console.ReadLine();



            Dictionary<string, object> settingsES = new Dictionary<string, object>();
            settingsES.Add("popsize", 3);          // ∈ {2,...,200}
            settingsES.Add("lambda", 1);            // ∈ {1,...,200}
            settingsES.Add("roh", 2);               // ∈ {1,...,popsize}  . in hyperoptimization, express as percentage of lambda
            settingsES.Add("x0sampling", 0);        // ∈ {0,1}  0=uniform, 1=gaussian
            settingsES.Add("stepsize0", 10);       // ∈ [0.01, 10]
            settingsES.Add("stepsize", 1.8);        // ∈ [0.01, 10]
            settingsES.Add("tauc", 0.01);              // ∈ [0.01, 50]
            settingsES.Add("selmode", 1);
            //settingsES.Add("pmut_int", 0.1);        // ∈ [0.01, 0.99] 
            SimpleES[] es = new SimpleES[seeds];
            for (int i = 0; i < seeds; i++)
            {
                es[i] = new SimpleES(lb, ub, xint, evalcount, testfunc, i, settingsES);
                es[i].solve();
                optis[i] = es[i].get_fxoptimum();
            }
            Console.WriteLine("es average: {0}", optis.Average());
            Console.ReadKey();



            //Dictionary<string, object> settingsPSO = new Dictionary<string, object>();
            //settingsPSO.Add("popsize", 29);        // popsize                           ∈ {4,..., 100}
            //settingsPSO.Add("chi", 0.32903);            // constriction coefficient         ∈ [0.001, 1]
            //settingsPSO.Add("phi", 10.87919);              // attraction to best particle     ∈ [0.01, 50]
            //settingsPSO.Add("v0max", 11.24591);          // max velocity at initialisation. fraction of domain. ∈ [0.01, 10]
            //settingsPSO.Add("x0samplingmode", 0);   // 0 = uniform, 1 = gaussian        ∈ [0.01, 10]
            //settingsPSO.Add("pxupdatemode", 0);     //0 = update after population. 1 = update after each evaluation ∈ [0.01, 10]
            //settingsPSO.Add("s0", 1);             //initial step size in case of gaussian x0 ∈ [0.01, 10]
            //settingsPSO.Add("psomode", 0);      //0 = fipso, 1 = inertia, 2 = constriction
            //settingsPSO.Add("phi1", 1.05);      //attraction own best
            //settingsPSO.Add("phi2", 2.05);      //attraction global best
            //settingsPSO.Add("intmut", 0.5);     //mutation probability for integer variables
            //settingsPSO.Add("mutstdev", 0.3);
            //PSO[] pso = new PSO[seeds];
            //for (int i = 0; i < seeds; i++)
            //{
            //    pso[i] = new PSO(lb, ub, xint, evalcount, testfunc, i, settingsPSO);
            //    pso[i].solve();
            //    optis[i] = pso[i].get_fxoptimum();
            //}
            //Console.WriteLine("pso average: {0}", optis.Average());
            //Console.WriteLine("pso min: {0}", optis.Min());
            //Console.WriteLine("pso max: {0}", optis.Max());

            //Console.ReadKey();

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
                return SO.L_Ackley(dbl);
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


}
