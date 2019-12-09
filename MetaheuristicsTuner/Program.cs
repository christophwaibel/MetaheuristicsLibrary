using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MetaheuristicsLibrary.TestFunctions;
using MetaheuristicsLibrary;
using System.IO;


namespace MetaheuristicsTuner
{
    class Program
    {
        static void Main(string[] args)
        {
            //______________________________________________________________________________________
            ////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////

            //SGA PARAMETERS
            double[] hp = new double[7];
            hp[0] = 15.20115;   // popsize
            hp[1] = 90.24627;   // k
            hp[2] = 1;  // pcross
            hp[3] = 0.80054;    // pmut
            hp[4] = 1.80471;    // d
            hp[5] = 0.01;   // r
            hp[6] = 0.26437;	// elites


            ////ES parameters
            //double[] hp = new double[8];
            //// cluster a
            //hp[0] = 2;
            //hp[1] = 1;
            //hp[2] = 0.16085;
            //hp[3] = 0.07627;
            //hp[4] = 0.01;
            //hp[5] = 0.2288;
            //hp[6] = 0.01136;
            //hp[7] = 0.92639;


            ////PSO PARAmeters
            //double[] hp = new double[7];
            //// cluster  
            //hp[0] = 9.03316;
            //hp[1] = 0.001;
            //hp[2] = 18.24318;
            //hp[3] = 0;
            //hp[4] = 1;
            //hp[5] = 0;
            //hp[6] = 3.22193;



            ////FIPSO PARAMeters
            //double[] hp = new double[5];
            //// cluster a
            //hp[0] = 12.16184;
            //hp[1] = 0.55492;
            //hp[2] = 9.31087;
            //hp[3] = 0.59655;
            //hp[4] = 1;

            TestHyperParam("SGA", hp, 13, 100, 30);

            ////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////






            //______________________________________________________________________________________
            ////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////

            //HyperParameterClassifier.IdentifyTwoBestParameterSets("FIPS");    //"SGA", "ES", "PSO", "FIPS"

            ////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////







            //______________________________________________________________________________________
            ////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////

            //HyperParameterOptimization.TuneSolver("SGA", 20, 100);

            ////////////////////////////////////////////////////////////////////////////////////////
            ////////////////////////////////////////////////////////////////////////////////////////
        }




        
        /// <summary>
        /// Testing a hyper-parameter configuration for one test function
        /// </summary>
        /// <param name="args"></param>
        static void TestHyperParamOneFuncOnly()
        {
            string solver = "SGA";
            //double[] hp;
            int testfuncIndex = 19;
            int testfuncdim = 35;
            int evalbudgetmultipl = 35;
            int rerunsTestFunc = 35;

            //SGA PARAMETERS
            double[] hp = new double[7];
            hp[0] = 15.20115;   // popsize
            hp[1] = 90.24627;   // k
            hp[2] = 1;  // pcross
            hp[3] = 0.80054;    // pmut
            hp[4] = 1.80471;    // d
            hp[5] = 0.01;   // r
            hp[6] = 0.26437;	// elites


            ////ES parameters
            //double[] hp = new double[8];
            //// cluster a
            //hp[0] = 2;
            //hp[1] = 1;
            //hp[2] = 0.16085;
            //hp[3] = 0.07627;
            //hp[4] = 0.01;
            //hp[5] = 0.2288;
            //hp[6] = 0.01136;
            //hp[7] = 0.92639;


            ////PSO PARAmeters
            //double[] hp = new double[7];
            //// cluster  
            //hp[0] = 9.03316;
            //hp[1] = 0.001;
            //hp[2] = 18.24318;
            //hp[3] = 0;
            //hp[4] = 1;
            //hp[5] = 0;
            //hp[6] = 3.22193;



            ////FIPSO PARAMeters
            //double[] hp = new double[5];
            //// cluster a
            //hp[0] = 12.16184;
            //hp[1] = 0.55492;
            //hp[2] = 9.31087;
            //hp[3] = 0.59655;
            //hp[4] = 1;


            //TestHyperParamOneFuncOnly("FIPSO", hp, 19, 35, 30, 30);





            //int evalbudgetmultiplier = 100;
            //int testfuncdim = 20;                //tuned to this n
            //int rerunsTestFuncs = 30;            // 30
            HyperParameterFunctions hf = new HyperParameterFunctions(testfuncdim, rerunsTestFunc, solver, evalbudgetmultipl);   //problem dim; reruns of testfuncs; solver to be tuned
            Func<double[], double[]>[] manyhyperfuncs = new Func<double[], double[]>[20];
            manyhyperfuncs[0] = hf.HyperFunc_B_Perm0db_dist;
            manyhyperfuncs[1] = hf.HyperFunc_B_RotHypEll_dist;
            manyhyperfuncs[2] = hf.HyperFunc_B_Sphere_dist;
            manyhyperfuncs[3] = hf.HyperFunc_B_SumSquares_dist;
            manyhyperfuncs[4] = hf.HyperFunc_B_Trid_dist;
            manyhyperfuncs[5] = hf.HyperFunc_L_Ackley_dist;
            manyhyperfuncs[6] = hf.HyperFunc_L_Griewank_dist;
            manyhyperfuncs[7] = hf.HyperFunc_L_Levy_dist;
            manyhyperfuncs[8] = hf.HyperFunc_L_Rastrigin_dist;
            manyhyperfuncs[9] = hf.HyperFunc_L_Schwefel_dist;
            manyhyperfuncs[10] = hf.HyperFunc_O_PermDB_dist;
            manyhyperfuncs[11] = hf.HyperFunc_O_StyblinskiTang_dist;
            manyhyperfuncs[12] = hf.HyperFunc_P_Zhakarov_dist;
            manyhyperfuncs[13] = hf.HyperFunc_V_DixonPrice_dist;
            manyhyperfuncs[14] = hf.HyperFunc_V_Rosenbrock_dist;
            manyhyperfuncs[15] = hf.HyperFunc_L_Levy_Edge_dist;
            manyhyperfuncs[16] = hf.HyperFunc_L_Schwefel_Edge_dist;
            manyhyperfuncs[17] = hf.HyperFunc_O_StyblinskiTang_Edge_dist;
            manyhyperfuncs[18] = hf.HyperFunc_V_DixonPrice_Edge_dist;
            manyhyperfuncs[19] = hf.HyperFunc_V_Rosenbrock_Edge_dist;

            double [] dist = manyhyperfuncs[testfuncIndex](hp);
            foreach (double fx in dist)
            {
                Console.WriteLine(fx);
            }
            Console.WriteLine();
            Console.WriteLine("expected value: {0}", dist.Average());
            Console.ReadKey();


        }

        /// <summary>
        /// Testing a hyper-parameter configuration for all test functions
        /// </summary>
        /// <param name="args"></param>
        static void TestHyperParam(string solver, double [] hp, int testfuncdim, int evalbudgetmultipl, int rerunsTestFuncs)
        {
            // solver = "SGA";
            ////double[] hp;
            // testfuncdim = 13;
            //int evalbudgetmultipl = 100;
            //int rerunsTestFuncs = 30;


            ////SGA PARAMETERS
            //double[] hp = new double[7];
            //hp[0] = 15.20115;   // popsize
            //hp[1] = 90.24627;   // k
            //hp[2] = 1;  // pcross
            //hp[3] = 0.80054;    // pmut
            //hp[4] = 1.80471;    // d
            //hp[5] = 0.01;   // r
            //hp[6] = 0.26437;	// elites


            ////ES parameters
            //double[] hp = new double[8];
            //// cluster a
            //hp[0] = 2;
            //hp[1] = 1;
            //hp[2] = 0.16085;
            //hp[3] = 0.07627;
            //hp[4] = 0.01;
            //hp[5] = 0.2288;
            //hp[6] = 0.01136;
            //hp[7] = 0.92639;


            ////PSO PARAmeters
            //double[] hp = new double[7];
            //// cluster  
            //hp[0] = 9.03316;
            //hp[1] = 0.001;
            //hp[2] = 18.24318;
            //hp[3] = 0;
            //hp[4] = 1;
            //hp[5] = 0;
            //hp[6] = 3.22193;



            ////FIPSO PARAMeters
            //double[] hp = new double[5];
            //// cluster a
            //hp[0] = 12.16184;
            //hp[1] = 0.55492;
            //hp[2] = 9.31087;
            //hp[3] = 0.59655;
            //hp[4] = 1;



            //TestHyperParam("SGA", hp, 13, 100, 30);





            //int evalbudgetmultiplier = 100;
            //int testfuncdim = 20;                //tuned to this n
            //int rerunsTestFuncs = 30;            // 30
            HyperParameterFunctions hf = new HyperParameterFunctions(testfuncdim, rerunsTestFuncs, solver, evalbudgetmultipl);   //problem dim; reruns of testfuncs; solver to be tuned
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

            double[] means = new double[manyhyperfuncs.Length];
            for (int i = 0; i < manyhyperfuncs.Length; i++)
            {
                means[i] = manyhyperfuncs[i](hp);
                Console.WriteLine(means[i]);
            }
            Console.ReadKey();


        }


    }


 
}
