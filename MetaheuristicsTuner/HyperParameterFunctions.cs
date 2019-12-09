using System;
using System.Collections.Generic;
using System.Linq;
using System.IO;
using MetaheuristicsLibrary.TestFunctions;


namespace MetaheuristicsTuner
{
    internal class HyperParameterFunctions
    {
        int testfuncDim;
        int runsperTestfunc;
        string solver;
        int evalbudget;
        List<string> logs = new List<string>();


        public void printLogs(string basepath, int hh)
        {
            string fileName = basepath + solver + @"_hf_" + hh + "_AllEvals.txt";
            using (FileStream fs = new FileStream(fileName, FileMode.Append, FileAccess.Write))
            using (StreamWriter sw = new StreamWriter(fs))
            {
                foreach (string l in logs)
                {
                    sw.WriteLine(l);
                }
            }
        }


        internal HyperParameterFunctions(int dim, int reruns, string solvername, int evaluationBudgetMultiplier)
        {
            testfuncDim = dim;
            runsperTestfunc = reruns;
            solver = solvername;
            evalbudget = evaluationBudgetMultiplier;

        }


        /// <summary>
        /// Returns Hyper Parameter (HP) information of a solver that can be tuned
        /// </summary>
        /// <param name="dvar">Number of HP of a solver</param>
        /// <param name="lb">lower bound of each HP</param>
        /// <param name="ub">upper bound of each HP</param>
        /// <param name="xint">indicator whether HP is an integer. False, if it is a continuous parameter</param>
        /// <param name="solver">Solver name. Options: SGA, ES, FIPSO, PSO</param>
        public void getXandBounds(out int dvar, out double[] lb, out double[] ub, out bool[] xint, string solver)
        {
            switch (solver)
            {
                default:    //SGA
                    dvar = 7;
                    lb = new double[dvar];
                    ub = new double[dvar];
                    xint = new bool[dvar];
                    lb[0] = 5;          //simpleGAsettings.Add("popsize", Convert.ToInt32(x[0]));            //x[0]  50
                    ub[0] = 200;
                    xint[0] = false;
                    lb[1] = 1;          //simpleGAsettings.Add("k", Convert.ToInt32(x[1]));//x[1]  6
                    ub[1] = 100;
                    xint[1] = false;
                    lb[2] = 0.01;       //simpleGAsettings.Add("pcross", x[2]);           //x[2]  0.7
                    ub[2] = 1.0;
                    xint[2] = false;
                    lb[3] = 0;       //simpleGAsettings.Add("pmut", x[3]);             //x[3]  0.3
                    ub[3] = 1.0;
                    xint[3] = false;
                    lb[4] = 0.01;       //simpleGAsettings.Add("d", x[4]);                //x[4]  0.1
                    ub[4] = 2.0;
                    xint[4] = false;
                    lb[5] = 0.01;       //simpleGAsettings.Add("r", x[5]);                //x[5]  0.1
                    ub[5] = 2.0;
                    xint[5] = false;
                    lb[6] = 0.0;       //simpleGAsettings.Add("elite", x[6]);                //x[6]  
                    ub[6] = 1.0;
                    xint[6] = false;
                    break;
                case "ES":
                    dvar = 8;
                    lb = new double[dvar];
                    ub = new double[dvar];
                    xint = new bool[dvar];
                    lb[0] = 2;      // "popsize" ∈ {2,...,200}
                    ub[0] = 200;
                    xint[0] = false;
                    lb[1] = 1;      // "lambda" (offspring) ∈ {1,...,200}
                    ub[1] = 200;
                    xint[1] = false;
                    lb[2] = 0;      // "roh" ∈ {1,...,popsize} in %
                    ub[2] = 1.0;
                    xint[2] = false;
                    lb[3] = 0;      // "x0sampling" ∈ {0,1}  0=uniform, 1=gaussian
                    ub[3] = 1.0;
                    xint[3] = false;
                    lb[4] = 0.01;   // "stepsize0" ∈ [0.01, 10]
                    ub[4] = 10;
                    xint[4] = false;
                    lb[5] = 0.01;   // "stepsie" ∈ [0.01, 10]
                    ub[5] = 10;
                    xint[5] = false;
                    lb[6] = 0.01;   // "tauc" ∈ [0.01, 50]
                    ub[6] = 50;
                    xint[6] = false;
                    lb[7] = 0;      // "selmode" ∈ {0,1}
                    ub[7] = 1;
                    xint[7] = false;
                    break;
                case "FIPSO":
                    dvar = 5;
                    lb = new double[dvar];
                    ub = new double[dvar];
                    xint = new bool[dvar];
                    lb[0] = 4;
                    ub[0] = 200;
                    xint[0] = false;        // "popsize" ∈ {4,...,200}
                    lb[1] = 0.001;
                    ub[1] = 1;
                    xint[1] = false;        // "chi" constriction coefficient ∈ [0.001,1]
                    lb[2] = 0;
                    ub[2] = 50;
                    xint[2] = false;        // "phi" attraction to best particles ∈ [0,10]
                    lb[3] = 0;
                    ub[3] = 20;
                    xint[3] = false;        // "v0max" initial velocity multiplicator ∈ [0,20]
                    lb[4] = 0;
                    ub[4] = 1;
                    xint[4] = true;         //0 = update after population. 1 = update after each evaluation 
                    //lb[5] = 0;  
                    //ub[5] = 2;
                    //xint[5] = true;         //0 = fipso, 1 = inertia, 2 = constriction 
                    //lb[6] = 0;
                    //ub[6] = 5;      
                    //xint[6] = false;        //attraction own best ∈ [0,5]
                    //lb[7] = 0;
                    //ub[7] = 5;
                    //xint[7] = false;        //attraction global best ∈ [0,5]
                    break;
                case "PSO":
                    dvar = 7;
                    lb = new double[dvar];
                    ub = new double[dvar];
                    xint = new bool[dvar];
                    lb[0] = 4;
                    ub[0] = 200;
                    xint[0] = false;        // "popsize" ∈ {4,...,200}
                    lb[1] = 0.001;
                    ub[1] = 1;
                    xint[1] = false;        // "chi" constriction coefficient ∈ [0.001,1]
                    //lb[2] = 0;
                    //ub[2] = 50;
                    //xint[2] = false;        // "phi" attraction to best particles ∈ [0,50]
                    lb[2] = 0;
                    ub[2] = 20;
                    xint[2] = false;        // "v0max" initial velocity multiplicator ∈ [0,20]
                    lb[3] = 0;
                    ub[3] = 1;
                    xint[3] = true;         //0 = update after population. 1 = update after each evaluation 
                    lb[4] = 1;
                    ub[4] = 2;
                    xint[4] = true;         //0 = fipso, 1 = inertia, 2 = constriction. fipso deactivated here 
                    lb[5] = 0;
                    ub[5] = 5;
                    xint[5] = false;        //attraction own best ∈ [0,5]
                    lb[6] = 0;
                    ub[6] = 5;
                    xint[6] = false;        //attraction global best ∈ [0,5]
                    break;

                    //case "SA":

                    //    break;
            }

        }


        /// <summary>
        /// Computes the average optimal cost values found by a parametrized solver for a certain test function.
        /// This one is for ES, which has 8 parameters for hpo_x:
        /// [0] population size, [1] lambda, [2] roh, [3] x0sampling, [4] stepsize0, [5] stepsize, [6] tauc, [7] selection mode
        /// </summary>
        /// <param name="hpo_x">hyper parameters</param>
        /// <param name="testfunc">test function</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of the test function variables</param>
        /// <param name="ub">upper bound of the test function variables</param>
        /// <returns>Average minimal cost value that this parametrized solver can find for a certain test function</returns>
        private double meanES(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> settingsES = new Dictionary<string, object>();
            int simpleESpop = Convert.ToInt16(hpo_x[0]);
            int simpleESLambda = Convert.ToInt16(hpo_x[1]);
            int simpleESroh = Convert.ToInt16(hpo_x[2] * simpleESLambda) + 1; //from percentage to integer
            settingsES.Add("popsize", simpleESpop);          // ∈ {2,...,200}
            settingsES.Add("lambda", simpleESLambda);            // ∈ {1,...,200}
            settingsES.Add("roh", simpleESroh);               // ∈ {1,...,popsize}  . in hyperoptimization, express as percentage of lambda
            settingsES.Add("x0sampling", Convert.ToInt16(hpo_x[3]));        // ∈ {0,1}  0=uniform, 1=gaussian
            settingsES.Add("stepsize0", hpo_x[4]);       // ∈ [0.01, 10]
            settingsES.Add("stepsize", hpo_x[5]);        // ∈ [0.01, 10]
            settingsES.Add("tauc", hpo_x[6]);              // ∈ [0.01, 10]
            settingsES.Add("selmode", hpo_x[7]);        // 0 = random, 1=roulette wheel
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] simpleES = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];

            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                simpleES[i] = new MetaheuristicsLibrary.SingleObjective.EvolutionStrategy(lb, ub, xint, _maxfuncs, testfunc, i, settingsES);
                simpleES[i].solve();
                mins[i] = simpleES[i].get_fxoptimum();
            }
            return mins.Average();
        }


        /// <summary>
        /// Computes the average optimal cost values found by a parametrized solver for a certain test function.
        /// This one is for SGA, which has 7 parameters for hpo_x:
        /// [0] population size, [1] k, [2] pcross, [3] pmut, [4] d, [5] r, [6] elite
        /// </summary>
        /// <param name="hpo_x">hyper parameters</param>
        /// <param name="testfunc">test function</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of the test function variables</param>
        /// <param name="ub">upper bound of the test function variables</param>
        /// <returns>Average minimal cost value that this parametrized solver can find for a certain test function</returns>
        private double meanSGA(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> simpleGAsettings = new Dictionary<string, object>();
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] simpleGAs = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];
            int simplegapop = Convert.ToInt32(hpo_x[0]);
            simplegapop += simplegapop % 2;
            simpleGAsettings.Add("popsize", simplegapop);            //x[0]  50
            simpleGAsettings.Add("maxgen", Convert.ToInt32(Math.Floor(_maxfuncs / Convert.ToDouble(simplegapop))));            //no x. its a function of popsize and max evals.  50
            simpleGAsettings.Add("k", hpo_x[1]);                //x[1]  6
            simpleGAsettings.Add("pcross", hpo_x[2]);           //x[2]  0.7
            simpleGAsettings.Add("pmut", hpo_x[3]);             //x[3]  0.3
            simpleGAsettings.Add("d", hpo_x[4]);                //x[4]  0.1
            simpleGAsettings.Add("r", hpo_x[5]);                //x[5]  0.1
            int SGAelite = Convert.ToInt16(hpo_x[6] * (int)Math.Round((double)simplegapop / 2, 0)); //from percentage to integer
            simpleGAsettings.Add("elite", SGAelite);                //x[6]  0 - (popsize/2). 

            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                simpleGAs[i] = new MetaheuristicsLibrary.SingleObjective.GeneticAlgorithm(lb, ub, xint, _maxfuncs, testfunc, i, simpleGAsettings);
                simpleGAs[i].solve();
                mins[i] = simpleGAs[i].get_fxoptimum();
            }

            return mins.Average();
        }


        /// <summary>
        /// Computes the average optimal cost values found by a parametrized solver for a certain test function.
        /// This one is for FIPSO, which has 5 parameters for hpo_x:
        /// [0] population size, [1] chi, [2] phi, [3] v0max, [4] pxupdatemode
        /// </summary>
        /// <param name="hpo_x">hyper parameters</param>
        /// <param name="testfunc">test function</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of the test function variables</param>
        /// <param name="ub">upper bound of the test function variables</param>
        /// <returns>Average minimal cost value that this parametrized solver can find for a certain test function</returns>
        private double meanFIPSO(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> PSOsettings = new Dictionary<string, object>();
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] PSO = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];
            int fipsopop = Convert.ToInt16(hpo_x[0]);
            PSOsettings.Add("popsize", fipsopop);             //x[0] ∈ {4,..., 100}
            PSOsettings.Add("chi", hpo_x[1]);                 //x[1] ∈ [0.001, 1], constriction coefficient
            PSOsettings.Add("phi", hpo_x[2]);                 //x[2] ∈ [0, 50], attraction to best particle
            PSOsettings.Add("v0max", hpo_x[3]);               //x[3] ∈ [0, 10], initial velocity
            PSOsettings.Add("pxupdatemode", Convert.ToInt16(hpo_x[4])); //x[4] 0 = update after population. 1 = update after each evaluation 
            PSOsettings.Add("psomode", 0);
            //PSOsettings.Add("psomode", Convert.ToInt16(hpo_x[5])); //x[5] 0 = fipso, 1 = inertia, 2 = constriction 
            //PSOsettings.Add("phi1", hpo_x[6]);                //x[6] ∈ [0,5], attraction own best 
            //PSOsettings.Add("phi2", hpo_x[7]);                //x[7] ∈ [0,5], attraction global best 
            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                PSO[i] = new MetaheuristicsLibrary.SingleObjective.ParticleSwarmOptimization(lb, ub, xint, _maxfuncs, testfunc, i, PSOsettings);
                PSO[i].solve();
                mins[i] = PSO[i].get_fxoptimum();
            }

            return mins.Average();
        }


        /// <summary>
        /// Computes the average optimal cost values found by a parametrized solver for a certain test function.
        /// This one is for PSO, which has 7 parameters for hpo_x:
        /// [0] population size, [1] chi, [2] v0max, [3] pxupdatemode, [4] psomode, [5] phi1, [6] phi2
        /// </summary>
        /// <param name="hpo_x">hyper parameters</param>
        /// <param name="testfunc">test function</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of the test function variables</param>
        /// <param name="ub">upper bound of the test function variables</param>
        /// <returns>Average minimal cost value that this parametrized solver can find for a certain test function</returns>
        private double meanPSO(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> PSOsettings = new Dictionary<string, object>();
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] PSO = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];
            int fipsopop = Convert.ToInt16(hpo_x[0]);
            PSOsettings.Add("popsize", fipsopop);             //x[0] ∈ {4,..., 100}
            PSOsettings.Add("chi", hpo_x[1]);                 //x[1] ∈ [0.001, 1], constriction coefficient
            //PSOsettings.Add("phi", hpo_x[2]);                 //x[2] ∈ [0, 50], attraction to best particle
            PSOsettings.Add("v0max", hpo_x[2]);               //x[3] ∈ [0, 10], initial velocity
            PSOsettings.Add("pxupdatemode", Convert.ToInt16(hpo_x[3])); //x[4] 0 = update after population. 1 = update after each evaluation 
            PSOsettings.Add("psomode", Convert.ToInt16(hpo_x[4])); //x[5] 0 = fipso, 1 = inertia, 2 = constriction 
            PSOsettings.Add("phi1", hpo_x[5]);                //x[6] ∈ [0,5], attraction own best 
            PSOsettings.Add("phi2", hpo_x[6]);                //x[7] ∈ [0,5], attraction global best 
            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                PSO[i] = new MetaheuristicsLibrary.SingleObjective.ParticleSwarmOptimization(lb, ub, xint, _maxfuncs, testfunc, i, PSOsettings);
                PSO[i].solve();
                mins[i] = PSO[i].get_fxoptimum();
            }

            return mins.Average();
        }


        /// <summary>
        /// Returns the average minimal cost value found by a solver for a certain test function.
        /// The solver is set globally when instantiating a HyperParameterFunction object.
        /// </summary>
        /// <param name="x">hyper parameters for the solver</param>
        /// <param name="testfunc">test function to be tested</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of the test function variables</param>
        /// <param name="ub">upper bound of the test function variables</param>
        /// <returns>Average minimal cost value found in all runs,</returns>
        private double switchMean(double[] x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            switch (solver)
            {
                case "ES":
                    return meanES(x, testfunc, _maxfuncs, lb, ub);
                case "SGA":
                    return meanSGA(x, testfunc, _maxfuncs, lb, ub);
                case "FIPSO":
                    return meanFIPSO(x, testfunc, _maxfuncs, lb, ub);
                case "PSO":
                    return meanPSO(x, testfunc, _maxfuncs, lb, ub);
                default:
                    return 0.0;
            }
        }


        /// <summary>
        /// Computes a distribution of minimal cost values found in a number of runs of a parametrized solver.
        /// Solver is set globally when instantiating HyperParameterFunctions object. Number of runs as well.
        /// </summary>
        /// <param name="hpo_x">hyper parameters to be tested</param>
        /// <param name="testfunc">test function to be tested</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of test function variables</param>
        /// <param name="ub">upper bound of test function variables</param>
        /// <returns>A vector of minimal cost values found.</returns>
        private double[] distES(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> settingsES = new Dictionary<string, object>();
            int simpleESpop = Convert.ToInt16(hpo_x[0]);
            int simpleESLambda = Convert.ToInt16(hpo_x[1]);
            int simpleESroh = Convert.ToInt16(hpo_x[2] * simpleESLambda) + 1; //from percentage to integer
            settingsES.Add("popsize", simpleESpop);          // ∈ {2,...,200}
            settingsES.Add("lambda", simpleESLambda);            // ∈ {1,...,200}
            settingsES.Add("roh", simpleESroh);               // ∈ {1,...,popsize}  . in hyperoptimization, express as percentage of lambda
            settingsES.Add("x0sampling", Convert.ToInt16(hpo_x[3]));        // ∈ {0,1}  0=uniform, 1=gaussian
            settingsES.Add("stepsize0", hpo_x[4]);       // ∈ [0.01, 10]
            settingsES.Add("stepsize", hpo_x[5]);        // ∈ [0.01, 10]
            settingsES.Add("tauc", hpo_x[6]);              // ∈ [0.01, 10]
            settingsES.Add("selmode", hpo_x[7]);        // 0 = random, 1=roulette wheel
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] simpleES = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];

            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                simpleES[i] = new MetaheuristicsLibrary.SingleObjective.EvolutionStrategy(lb, ub, xint, _maxfuncs, testfunc, i, settingsES);
                simpleES[i].solve();
                mins[i] = simpleES[i].get_fxoptimum();
            }
            return mins;
        }


        /// <summary>
        /// Computes a distribution of minimal cost values found in a number of runs of a parametrized solver.
        /// Solver is set globally when instantiating HyperParameterFunctions object. Number of runs as well.
        /// </summary>
        /// <param name="hpo_x">hyper parameters to be tested</param>
        /// <param name="testfunc">test function to be tested</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of test function variables</param>
        /// <param name="ub">upper bound of test function variables</param>
        /// <returns>A vector of minimal cost values found.</returns>
        private double[] distSGA(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> simpleGAsettings = new Dictionary<string, object>();
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] simpleGAs = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];
            int simplegapop = Convert.ToInt32(hpo_x[0]);
            simplegapop += simplegapop % 2;
            simpleGAsettings.Add("popsize", simplegapop);            //x[0]  50
            simpleGAsettings.Add("maxgen", Convert.ToInt32(Math.Floor(_maxfuncs / Convert.ToDouble(simplegapop))));            //no x. its a function of popsize and max evals.  50
            simpleGAsettings.Add("k", hpo_x[1]);                //x[1]  6
            simpleGAsettings.Add("pcross", hpo_x[2]);           //x[2]  0.7
            simpleGAsettings.Add("pmut", hpo_x[3]);             //x[3]  0.3
            simpleGAsettings.Add("d", hpo_x[4]);                //x[4]  0.1
            simpleGAsettings.Add("r", hpo_x[5]);                //x[5]  0.1
            int SGAelite = Convert.ToInt16(hpo_x[6] * (int)Math.Round((double)simplegapop / 2, 0)); //from percentage to integer
            simpleGAsettings.Add("elite", SGAelite);                //x[6]  0 - (popsize/2). 

            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                simpleGAs[i] = new MetaheuristicsLibrary.SingleObjective.GeneticAlgorithm(lb, ub, xint, _maxfuncs, testfunc, i, simpleGAsettings);
                simpleGAs[i].solve();
                mins[i] = simpleGAs[i].get_fxoptimum();
            }

            return mins;
        }


        /// <summary>
        /// Computes a distribution of minimal cost values found in a number of runs of a parametrized solver.
        /// Solver is set globally when instantiating HyperParameterFunctions object. Number of runs as well.
        /// </summary>
        /// <param name="hpo_x">hyper parameters to be tested</param>
        /// <param name="testfunc">test function to be tested</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of test function variables</param>
        /// <param name="ub">upper bound of test function variables</param>
        /// <returns>A vector of minimal cost values found.</returns>
        private double[] distFIPSO(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> PSOsettings = new Dictionary<string, object>();
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] PSO = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];
            int fipsopop = Convert.ToInt16(hpo_x[0]);
            PSOsettings.Add("popsize", fipsopop);             //x[0] ∈ {4,..., 100}
            PSOsettings.Add("chi", hpo_x[1]);                 //x[1] ∈ [0.001, 1], constriction coefficient
            PSOsettings.Add("phi", hpo_x[2]);                 //x[2] ∈ [0, 50], attraction to best particle
            PSOsettings.Add("v0max", hpo_x[3]);               //x[3] ∈ [0, 10], initial velocity
            PSOsettings.Add("pxupdatemode", Convert.ToInt16(hpo_x[4])); //x[4] 0 = update after population. 1 = update after each evaluation 
            PSOsettings.Add("psomode", 0);
            //PSOsettings.Add("psomode", Convert.ToInt16(hpo_x[5])); //x[5] 0 = fipso, 1 = inertia, 2 = constriction 
            //PSOsettings.Add("phi1", hpo_x[6]);                //x[6] ∈ [0,5], attraction own best 
            //PSOsettings.Add("phi2", hpo_x[7]);                //x[7] ∈ [0,5], attraction global best 
            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                PSO[i] = new MetaheuristicsLibrary.SingleObjective.ParticleSwarmOptimization(lb, ub, xint, _maxfuncs, testfunc, i, PSOsettings);
                PSO[i].solve();
                mins[i] = PSO[i].get_fxoptimum();
            }

            return mins;
        }


        /// <summary>
        /// Computes a distribution of minimal cost values found in a number of runs of a parametrized solver.
        /// Solver is set globally when instantiating HyperParameterFunctions object. Number of runs as well.
        /// </summary>
        /// <param name="hpo_x">hyper parameters to be tested</param>
        /// <param name="testfunc">test function to be tested</param>
        /// <param name="_maxfuncs">evaluation budget</param>
        /// <param name="lb">lower bound of test function variables</param>
        /// <param name="ub">upper bound of test function variables</param>
        /// <returns>A vector of minimal cost values found.</returns>
        private double[] distPSO(double[] hpo_x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            int dvar = testfuncDim;
            int runs = runsperTestfunc;

            bool[] xint = new bool[dvar];
            Dictionary<string, object> PSOsettings = new Dictionary<string, object>();
            MetaheuristicsLibrary.SingleObjective.SingleObjective[] PSO = new MetaheuristicsLibrary.SingleObjective.SingleObjective[runs];
            int fipsopop = Convert.ToInt16(hpo_x[0]);
            PSOsettings.Add("popsize", fipsopop);             //x[0] ∈ {4,..., 100}
            PSOsettings.Add("chi", hpo_x[1]);                 //x[1] ∈ [0.001, 1], constriction coefficient
            //PSOsettings.Add("phi", hpo_x[2]);                 //x[2] ∈ [0, 50], attraction to best particle
            PSOsettings.Add("v0max", hpo_x[2]);               //x[3] ∈ [0, 10], initial velocity
            PSOsettings.Add("pxupdatemode", Convert.ToInt16(hpo_x[3])); //x[4] 0 = update after population. 1 = update after each evaluation 
            PSOsettings.Add("psomode", Convert.ToInt16(hpo_x[4])); //x[5] 0 = fipso, 1 = inertia, 2 = constriction 
            PSOsettings.Add("phi1", hpo_x[5]);                //x[6] ∈ [0,5], attraction own best 
            PSOsettings.Add("phi2", hpo_x[6]);                //x[7] ∈ [0,5], attraction global best 
            double[] mins = new double[runs];
            for (int i = 0; i < runs; i++)
            {
                PSO[i] = new MetaheuristicsLibrary.SingleObjective.ParticleSwarmOptimization(lb, ub, xint, _maxfuncs, testfunc, i, PSOsettings);
                PSO[i].solve();
                mins[i] = PSO[i].get_fxoptimum();
            }

            return mins;
        }


        /// <summary>
        /// Computes minimal cost values found by a number of runs for a parametrized solver. Solver is set globally. Runs as well
        /// </summary>
        /// <param name="x">hyper parameters to be tested</param>
        /// <param name="testfunc">test function to be tested</param>
        /// <param name="_maxfuncs">evaluation budget for the test function</param>
        /// <param name="lb">lower bound of test function variables</param>
        /// <param name="ub">upper bound of test function variable</param>
        /// <returns>A vector of minimal cost values found for all runs of a parametrized solver</returns>
        private double[] switchDist(double[] x, Func<double[], double> testfunc, int _maxfuncs, double[] lb, double[] ub)
        {
            switch (solver)
            {
                case "ES":
                    return distES(x, testfunc, _maxfuncs, lb, ub);
                case "SGA":
                    return distSGA(x, testfunc, _maxfuncs, lb, ub);
                case "FIPSO":
                    return distFIPSO(x, testfunc, _maxfuncs, lb, ub);
                case "PSO":
                    return distPSO(x, testfunc, _maxfuncs, lb, ub);
                default:
                    return new double[lb.Length];
            }
        }








        internal double HyperFunc_L_Rastrigin(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Rastrigin;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.12;
                ub[i] = 5.12;
            }
            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Levy(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Levy;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -10;
                ub[i] = 10;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Griewank(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Griewank;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -600;
                ub[i] = 600;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Ackley(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Ackley;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -32.768;
                ub[i] = 32.768;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Schwefel(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Schwefel;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -500.0;
                ub[i] = 500.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_P_Zhakarov(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.P_Zakharov;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.0;
                ub[i] = 10.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_V_Rosenbrock(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_Rosenbrock;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.0;
                ub[i] = 10.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_V_DixonPrice(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_DixonPrice;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -10.0;
                ub[i] = 10.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_Sphere(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Sphere;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.12;
                ub[i] = 5.12;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_SumSquares(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_SumSquares;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -10.0;
                ub[i] = 10.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_Trid(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Trid;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = Math.Pow(dvar, 2) * -1;
                ub[i] = Math.Pow(dvar, 2);
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_RotHypEll(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_RotHypEll;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -65.536;
                ub[i] = 65.536;
            }
            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_Perm0db(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Perm0db;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -dvar;
                ub[i] = dvar;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_O_PermDB(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.O_PermDB;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -dvar;
                ub[i] = dvar;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_O_StyblinskiTang(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.O_StyblinskiTang;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5;
                ub[i] = 5;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }








        internal double HyperFunc_L_Rastrigin_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Rastrigin;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 5.12 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Levy_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Levy;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 10 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
            //double m = switchMean(x, testfunc, _maxfuncs, lb, ub);
            //Console.WriteLine(m);
            //logs.Add(Convert.ToString(m));
            //return m;
        }

        internal double HyperFunc_L_Griewank_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Griewank;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 600 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Ackley_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Ackley;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 32.768 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_L_Schwefel_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Schwefel;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 500.0 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_P_Zhakarov_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.P_Zakharov;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 15.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_V_Rosenbrock_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_Rosenbrock;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 15.0;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_V_DixonPrice_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_DixonPrice;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 10.0 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_Sphere_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Sphere;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 5.12 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_SumSquares_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_SumSquares;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 10.0 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_B_RotHypEll_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_RotHypEll;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 65.536 * 2;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double HyperFunc_O_StyblinskiTang_Edge(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.O_StyblinskiTang;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -3;
                ub[i] = 8;
            }

            return switchMean(x, testfunc, _maxfuncs, lb, ub);
        }













        internal double[] HyperFunc_L_Rastrigin_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Rastrigin;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.12;
                ub[i] = 5.12;
            }
            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Levy_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Levy;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -10;
                ub[i] = 10;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Griewank_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Griewank;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -600;
                ub[i] = 600;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Ackley_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Ackley;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -32.768;
                ub[i] = 32.768;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Schwefel_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Schwefel;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -500.0;
                ub[i] = 500.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_P_Zhakarov_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.P_Zakharov;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.0;
                ub[i] = 10.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_V_Rosenbrock_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_Rosenbrock;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.0;
                ub[i] = 10.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_V_DixonPrice_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_DixonPrice;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -10.0;
                ub[i] = 10.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_Sphere_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Sphere;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5.12;
                ub[i] = 5.12;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_SumSquares_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_SumSquares;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -10.0;
                ub[i] = 10.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_Trid_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Trid;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = Math.Pow(dvar, 2) * -1;
                ub[i] = Math.Pow(dvar, 2);
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_RotHypEll_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_RotHypEll;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -65.536;
                ub[i] = 65.536;
            }
            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_Perm0db_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Perm0db;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -dvar;
                ub[i] = dvar;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_O_PermDB_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.O_PermDB;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -dvar;
                ub[i] = dvar;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_O_StyblinskiTang_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.O_StyblinskiTang;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -5;
                ub[i] = 5;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }








        internal double[] HyperFunc_L_Rastrigin_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Rastrigin;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 5.12 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Levy_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Levy;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 10 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
            //double m = switchMean(x, testfunc, _maxfuncs, lb, ub);
            //Console.WriteLine(m);
            //logs.Add(Convert.ToString(m));
            //return m;
        }

        internal double[] HyperFunc_L_Griewank_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Griewank;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 600 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Ackley_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Ackley;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 32.768 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_L_Schwefel_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.L_Schwefel;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 500.0 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_P_Zhakarov_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.P_Zakharov;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 15.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_V_Rosenbrock_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_Rosenbrock;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 15.0;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_V_DixonPrice_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.V_DixonPrice;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 10.0 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_Sphere_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_Sphere;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 5.12 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_SumSquares_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_SumSquares;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 10.0 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_B_RotHypEll_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.B_RotHypEll;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 65.536 * 2;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }

        internal double[] HyperFunc_O_StyblinskiTang_Edge_dist(double[] x)
        {
            int dvar = testfuncDim;

            Func<double[], double> testfunc = SO.O_StyblinskiTang;
            int _maxfuncs = (dvar + 1) * evalbudget;   //max func evals
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = -3;
                ub[i] = 8;
            }

            return switchDist(x, testfunc, _maxfuncs, lb, ub);
        }
    }
}
