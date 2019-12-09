using System;
using System.Collections.Generic;
using System.Threading.Tasks;
using MetaheuristicsLibrary.SingleObjective;
using System.IO;


namespace MetaheuristicsTuner
{
    internal static class HyperParameterOptimization
    {
        /// <summary>
        /// Tuning hyper-parameters of a solver
        /// </summary>
        /// <param name="solver">Options: SGA, ES, FIPSO, PSO</param>
        /// <param name="testfuncdim">dimensionality n</param>
        /// <param name="evalbudgetmultilp">evaluation budget per test function</param>
        internal static void TuneSolver(string solver, int testfuncdim, int evalbudgetmultilp)
        {
            //string solver = "PSO";               //choosing solver to be tuned. string: "SGA", "ES", "FIPSO", "SA"
            int maxfunc = 5000;                  //func calls of the meta-optimizer. 5000
                                                 //int testfuncdim = 20;                //tuned to this n
            int rerunsMeta = 10;                 // 10
            int rerunsTestFuncs = 30;            // 30


            HyperParameterFunctions hf = new HyperParameterFunctions(testfuncdim, rerunsTestFuncs, solver, evalbudgetmultilp);   //problem dim; reruns of testfuncs; solver to be tuned
            int[] seeds_in = new int[rerunsMeta]; //reruns of meta-optimizer
            for (int i = 0; i < seeds_in.Length; i++) seeds_in[i] = i;

            int dvar;
            double[] ub, lb;
            bool[] xint;
            hf.getXandBounds(out dvar, out lb, out ub, out xint, solver);


            Console.WriteLine();
            Console.WriteLine(@"///////////////////////////////////////");
            Console.WriteLine("Hyperparameter optimization of...: {0}", solver);
            Console.WriteLine(@"tuned to n = {0}", testfuncdim);
            Console.WriteLine(@"Meta-optimizers: SGA and ES. Each with {0} re-runs", seeds_in.Length);
            Console.WriteLine(@"///////////////////////////////////////");
            Console.WriteLine();
            Console.WriteLine(@"Please enter an existing output path for results. E.g.: c:\n13\");
            string basepath = Console.ReadLine();
            Console.WriteLine();
            Console.WriteLine(@"How many threads would you like to assign to this? Please enter an integer between 1 and 20.");
            int maxthreads = Convert.ToInt16(Console.ReadLine());
            Console.WriteLine();
            Console.WriteLine("Thanks! Starting... This might take a couple of days.");





            //Func<double[], double>[] manyhyperfuncs = new Func<double[], double>[1];
            //manyhyperfuncs[0] = hf.HyperFunc_B_Perm0db;
            //manyhyperfuncs[1] = hf.HyperFunc_B_RotHypEll;
            //manyhyperfuncs[2] = hf.HyperFunc_L_Ackley;
            //manyhyperfuncs[3] = hf.HyperFunc_O_PermDB;
            //manyhyperfuncs[4] = hf.HyperFunc_L_Schwefel_Edge;
            //manyhyperfuncs[0] = hf.HyperFunc_L_Levy_Edge;

            Func<double[], double>[] manyhyperfuncs = new Func<double[], double>[20];
            manyhyperfuncs[0] = hf.HyperFunc_B_Perm0db;
            manyhyperfuncs[1] = hf.HyperFunc_B_RotHypEll;
            manyhyperfuncs[2] = hf.HyperFunc_B_Sphere;
            manyhyperfuncs[3] = hf.HyperFunc_B_SumSquares;
            manyhyperfuncs[4] = hf.HyperFunc_B_Trid;
            manyhyperfuncs[5] = hf.HyperFunc_L_Ackley;
            manyhyperfuncs[6] = hf.HyperFunc_L_Griewank;
            manyhyperfuncs[7] = hf.HyperFunc_L_Levy;
            manyhyperfuncs[8] = hf.HyperFunc_L_Rastrigin;
            manyhyperfuncs[9] = hf.HyperFunc_L_Schwefel;
            manyhyperfuncs[10] = hf.HyperFunc_O_PermDB;
            manyhyperfuncs[11] = hf.HyperFunc_O_StyblinskiTang;
            manyhyperfuncs[12] = hf.HyperFunc_P_Zhakarov;
            manyhyperfuncs[13] = hf.HyperFunc_V_DixonPrice;
            manyhyperfuncs[14] = hf.HyperFunc_V_Rosenbrock;
            manyhyperfuncs[15] = hf.HyperFunc_L_Levy_Edge;
            manyhyperfuncs[16] = hf.HyperFunc_L_Schwefel_Edge;
            manyhyperfuncs[17] = hf.HyperFunc_O_StyblinskiTang_Edge;
            manyhyperfuncs[18] = hf.HyperFunc_V_DixonPrice_Edge;
            manyhyperfuncs[19] = hf.HyperFunc_V_Rosenbrock_Edge;




            //manyhyperfuncs[15] = hf.HyperFunc_B_RotHypEll_Edge;
            //manyhyperfuncs[16] = hf.HyperFunc_B_Sphere_Edge;
            //manyhyperfuncs[17] = hf.HyperFunc_B_SumSquares_Edge;
            //manyhyperfuncs[18] = hf.HyperFunc_L_Ackley_Edge;
            //manyhyperfuncs[19] = hf.HyperFunc_L_Griewank_Edge;
            //manyhyperfuncs[21] = hf.HyperFunc_L_Rastrigin_Edge;
            //manyhyperfuncs[24] = hf.HyperFunc_P_Zhakarov_Edge;




            //meta optimizer SGA
            Dictionary<string, object> settingsSGA = new Dictionary<string, object>();
            settingsSGA.Add("maxgen", 10000);
            settingsSGA.Add("popsize", 20);
            settingsSGA.Add("elite", 1);
            settingsSGA.Add("pcross", 1);
            settingsSGA.Add("pmut", 0.3);
            settingsSGA.Add("d", 1.1);
            settingsSGA.Add("r", 0.1);
            settingsSGA.Add("k", 6);


            //meta optimizer FIPSO
            Dictionary<string, object> settingsPSO = new Dictionary<string, object>();

            //meta optimizer ES
            Dictionary<string, object> settingsES = new Dictionary<string, object>();
            settingsES.Add("popsize", 20);          // ∈ {2,...,200}
            settingsES.Add("lambda", 10);            // ∈ {1,...,200}
            settingsES.Add("roh", 2);               // ∈ {1,...,popsize}  . in hyperoptimization, express as percentage of lambda
            settingsES.Add("x0sampling", 1);        // ∈ {0,1}  0=uniform, 1=gaussian
            settingsES.Add("stepsize0", 5);       // ∈ [0.01, 10]
            settingsES.Add("stepsize", 2.2);        // ∈ [0.01, 10]
            settingsES.Add("tauc", 1);              // ∈ [0.01, 50]




            Parallel.For(0, manyhyperfuncs.Length, new ParallelOptions { MaxDegreeOfParallelism = maxthreads }, hh =>
            {
            //for(int hh=0; hh<manyhyperfuncs.Length; hh++){
            Console.WriteLine("started hh {0} of {1}", hh, manyhyperfuncs.Length - 1);
                List<string> log_SA = new List<string>();
                List<string> log_ES = new List<string>();
                List<string> log_PSO = new List<string>();
                List<string> log_SGA = new List<string>();
                string str;
                for (int iseeds = 0; iseeds < seeds_in.Length; iseeds++)
                {
                /////////////////////////////////////////            SGA          /////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////
                GeneticAlgorithm ga = new GeneticAlgorithm(lb, ub, xint, maxfunc, manyhyperfuncs[hh], iseeds, settingsSGA);
                    ga.solve();
                    ga.get_fxoptimum();
                    str = "Last call;" + Math.Round(ga.get_fxoptimum(), 4);
                    for (int xx = 0; xx < dvar; xx++)
                    {
                        str += ";" + Math.Round(ga.get_Xoptimum()[xx], 5);
                    }
                    str += ";SGA";
                    log_SGA.Add(str);
                    Console.WriteLine("Metaoptimizer SGA DONE, seed {0} and hh {1}", iseeds, hh);


                /////////////////////////////////////////            ES           /////////////////////////////////////////
                ///////////////////////////////////////////////////////////////////////////////////////////////////////////
                EvolutionStrategy es = new EvolutionStrategy(lb, ub, xint, maxfunc, manyhyperfuncs[hh], iseeds, settingsSGA);
                    es.solve();
                    es.get_fxoptimum();
                    str = "Last call;" + Math.Round(es.get_fxoptimum(), 4);
                    for (int xx = 0; xx < dvar; xx++)
                    {
                        str += ";" + Math.Round(es.get_Xoptimum()[xx], 5);
                    }
                    str += ";ES";
                    log_ES.Add(str);
                    Console.WriteLine("Metaoptimizer ES DONE, seed {0} and hh {1}", iseeds, hh);




                ///////////////////////////////////////////            FIPSO        /////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////////////////////////
                //PSO pso = new PSO(lb, ub, xint, maxfunc, manyhyperfuncs[hh], iseeds, settingsPSO);
                //pso.solve();
                //pso.get_fxoptimum();
                //str = "Last call;" + Math.Round(pso.get_fxoptimum(), 4);
                //for (int xx = 0; xx < dvar; xx++)
                //{
                //    str += ";" + Math.Round(pso.get_Xoptimum()[xx], 5);
                //}
                //str += ";PSO";
                //log_PSO.Add(str);
                //Console.WriteLine("Metaoptimizer PSO DONE, seed {0} and hh {1}", iseeds, hh);
            }


            /////////////////////////////////////////            Writing          ///////////////////////////////////
            ///////////////////////////////////////////////////////////////////////////////////////////////////////////
            string fileName = basepath + solver + @"_hf_" + hh + "_AllSeeds.txt";
                using (FileStream fs = new FileStream(fileName, FileMode.Append, FileAccess.Write))
                using (StreamWriter sw = new StreamWriter(fs))
                {
                    for (int i = 0; i < log_SA.Count; i++)
                        sw.WriteLine(log_SA[i]);
                    for (int i = 0; i < log_ES.Count; i++)
                        sw.WriteLine(log_ES[i]);
                    for (int i = 0; i < log_SGA.Count; i++)
                        sw.WriteLine(log_SGA[i]);
                    for (int i = 0; i < log_PSO.Count; i++)
                        sw.WriteLine(log_PSO[i]);
                }



                Console.WriteLine("Done hh {0} of {1}", hh, manyhyperfuncs.Length);
            });
            //}

            //hf.printLogs(basepath, 15);

            Console.WriteLine();
            Console.WriteLine(@"///////////////////////////////////////");
            Console.WriteLine("Done with everything, tuning the {0}", solver);
            Console.WriteLine(@"///////////////////////////////////////");
            Console.ReadKey();

        }

    }
}
