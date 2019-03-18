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

using System.Runtime.InteropServices;
using Excel = Microsoft.Office.Interop.Excel; 

namespace Tester
{
    class Program
    {

        /// <summary>
        /// Normalize cost values of BEO benchmarking
        /// </summary>
        /// <param name="args"></param>
        static void hhhMain(string[] args)
        {


            // load excel file and read tabs

            string Myfile = @"C:\_CHRIS\Test\Problem_15_180503.xlsx";
            //int n = 500; //p1 - 3
            //int n = 1400; // p 4- 7
            //int n = 1200; // p 8
            //int n = 1400; // p9-10
            //int n = 1900; // p 11
            //int n = 1800;   // p12
            //int n = 1100; // p13
            //int n = 2100; // p14
            int n = 1080; //p15
            int simgrad = 30; //simplex gradients , 30 for prob 15 only
            int runs = 20;
            string writenorm = @"C:\_CHRIS\Test\prob15.txt";

            Excel.Application xlApp ;
            Excel.Workbook xlWorkBook ;
            //Excel.Worksheet xlWorkSheet ;
            Excel.Range range ;



            xlApp = new Excel.Application();
            xlWorkBook = xlApp.Workbooks.Open(Myfile, 0, true, 5, "", "", true, Microsoft.Office.Interop.Excel.XlPlatform.xlWindows, "\t", false, false, 0, true, 1, 0);

            List<Excel.Worksheet> allsheets = new List<Excel.Worksheet>();
            allsheets.Add(xlWorkBook.Sheets["SA"]);
            allsheets.Add(xlWorkBook.Sheets["CMA_ES"]);
            allsheets.Add(xlWorkBook.Sheets["CRS2"]);
            allsheets.Add(xlWorkBook.Sheets["Direct"]);
            allsheets.Add(xlWorkBook.Sheets["ES_A"]);
            allsheets.Add(xlWorkBook.Sheets["ES_B"]);
            allsheets.Add(xlWorkBook.Sheets["ES"]);
            allsheets.Add(xlWorkBook.Sheets["FIPS_A"]);
            allsheets.Add(xlWorkBook.Sheets["FIPS_B"]);
            allsheets.Add(xlWorkBook.Sheets["FIPS"]);
            allsheets.Add(xlWorkBook.Sheets["GA"]);
            allsheets.Add(xlWorkBook.Sheets["PSO_A"]);
            allsheets.Add(xlWorkBook.Sheets["PSO_B"]);
            allsheets.Add(xlWorkBook.Sheets["PSO"]);
            allsheets.Add(xlWorkBook.Sheets["RBFOpt"]);
            allsheets.Add(xlWorkBook.Sheets["Sbplex"]);
            allsheets.Add(xlWorkBook.Sheets["SGA_WW"]);
            allsheets.Add(xlWorkBook.Sheets["SGA_A"]);
            allsheets.Add(xlWorkBook.Sheets["SGA_B"]);
            //allsheets.Add(xlWorkBook.Sheets["SilEye"]); 
            


            List<double[][]> allcost = new List<double[][]>();
            List<double[][]> costK = new List<double[][]>();
            allcost.Add(new double[runs][]); //SA
            allcost.Add(new double[runs][]); //CMA_ES
            allcost.Add(new double[runs][]); //CRS2
            allcost.Add(new double[runs][]); //Direct
            allcost.Add(new double[runs][]); //ES_A
            allcost.Add(new double[runs][]); //ES_B
            allcost.Add(new double[runs][]); //ES
            allcost.Add(new double[runs][]); //FIPS_A
            allcost.Add(new double[runs][]); //FIPS_B
            allcost.Add(new double[runs][]); //FIPS
            allcost.Add(new double[runs][]); //GA
            allcost.Add(new double[runs][]); //PSO_A
            allcost.Add(new double[runs][]); //PSO_B
            allcost.Add(new double[runs][]); //PSO
            allcost.Add(new double[runs][]); //RBFOpt
            allcost.Add(new double[runs][]); //Sbplex
            allcost.Add(new double[runs][]); //SGA_WW
            allcost.Add(new double[runs][]); //SGA_A
            allcost.Add(new double[runs][]); //SGA_B
            //allcost.Add(new double[runs][]); //SilEye
            for (int i = 0; i < allcost.Count; i++)
            {
                costK.Add(new double[runs][]);
                for (int c = 0; c < runs; c++)
                {
                    allcost[i][c] = new double[n];
                    costK[i][c] = new double[simgrad]; //100 simplex gradients only
                }
            }
            //get min and max of all
            double min = double.MaxValue;
            double max = double.MinValue;
            
            for (int i = 0; i < allsheets.Count; i++)
            {
                range = allsheets[i].UsedRange;
                for (int r = 0; r < n; r++) 
                {
                    for (int c = 0; c < runs; c++) 
                    {
                        allcost[i][c][r] = Convert.ToDouble((range.Cells[r+2, 7+(4*c)] as Excel.Range).Value2);      //+4
                        if (i == 3 && c > 0) allcost[i][c][r] = allcost[i][0][r];
                        if (r == n - 1)
                        {
                            if (allcost[i][c][r] < min) min = allcost[i][c][r];
                            if (allcost[i][c][r] > max) max = allcost[i][c][r];
                        }
                    }
                }
            }


            //reduce size to 100 simplex gradients
            // and normalize

            double k = n / simgrad;
            for (int i = 0; i < costK.Count; i++) //each solver
            {
                for (int c = 0; c < runs; c++) //each column (rerun)
                {
                    for (int s = 0; s < simgrad; s++) //each row (simplex gradient)
                    {
                        costK[i][c][s] = (allcost[i][c][Convert.ToInt32((s + 1) * k) - 1] - min) / (max - min);
                    }
                }
            }


            // write to text file

            string[] strsolver = new string[]{"SA","CMA-ES","CRS","Direct","ES_A","ES_B","ES","FIPS_A","FIPS_B","FIPS","GA","PSO_A","PSO_B","PSO","RBFOpt","Sbplex","SGA_WW","SGA_A","SGA_B","SilEye"};

            using (var writer = new StreamWriter(writenorm))
            {
                // Loop through ten numbers.
                for (int i = 0; i < costK.Count; i++)
                {
                    writer.WriteLine(strsolver[i]);
                    for (int r = 0; r < simgrad; r++) //rows, i.e. simplex gradients
                    {
                        string line = "";
                        for (int c = 0; c < runs; c++) //columns, i.e. reruns
                        {
                            line += Convert.ToString(costK[i][c][r]) + @";";
                            // Write format string to file.
                            //writer.Write("{0:0.0} ", i);
                        }
                        writer.WriteLine(line);
                    }
                }
            }

            xlWorkBook.Close(0);
            xlApp.Quit();
        }

        /// <summary>
        /// save only mins from Rhino Grasshopper runs.
        /// </summary>
        /// <param name="args"></param>
        static void grasshopperMain(string[] args)
        {


            //string func = "f16_LevyEdge";
            //string func = "f17_SchwefelEdge";
            //string func = "f18_StyblinskiEdge";
            //string func = "f19_DixonEdge";
            string func = "f20_RosenbrockEdge";
            int intN = 35;
            string n = "n" + intN;
            int k = 30; // simplex gradients. for n=35, make 30


            string[][] allPaths = new string[20][];
            allPaths[0] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\CRS2\", "*.txt", SearchOption.AllDirectories);
            allPaths[1] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\Direct\", "*.txt", SearchOption.AllDirectories);
            allPaths[2] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\ES_A\", "*.txt", SearchOption.AllDirectories);
            allPaths[3] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\ES_B\", "*.txt", SearchOption.AllDirectories);
            allPaths[4] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\ES_unt\", "*.txt", SearchOption.AllDirectories);
            allPaths[5] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\FIPS_A\", "*.txt", SearchOption.AllDirectories);
            allPaths[6] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\FIPS_B\", "*.txt", SearchOption.AllDirectories);
            allPaths[7] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\FIPS_unt\", "*.txt", SearchOption.AllDirectories);
            allPaths[8] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\GA\", "*.txt", SearchOption.AllDirectories);
            allPaths[9] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\PSO_A\", "*.txt", SearchOption.AllDirectories);
            allPaths[10] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\PSO_B\", "*.txt", SearchOption.AllDirectories);
            allPaths[11] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\PSO_unt\", "*.txt", SearchOption.AllDirectories);
            allPaths[12] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\SA\", "*.txt", SearchOption.AllDirectories);
            allPaths[13] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\SBplx\", "*.txt", SearchOption.AllDirectories);
            allPaths[14] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\Seye\", "*.txt", SearchOption.AllDirectories);
            allPaths[15] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\SGA_A\", "*.txt", SearchOption.AllDirectories);
            allPaths[16] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\SGA_B\", "*.txt", SearchOption.AllDirectories);
            allPaths[17] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\SGA_unt\", "*.txt", SearchOption.AllDirectories);
            allPaths[18] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\CMA_ES\", "*.txt", SearchOption.AllDirectories);
            allPaths[19] = Directory.GetFiles(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\RBFOpt\", "*.txt", SearchOption.AllDirectories);

            string[] solvernames = new string[20];
            solvernames[0] = "CRS2";
            solvernames[1] = "Direct";
            solvernames[2] = "ES_A";
            solvernames[3] = "ES_B";
            solvernames[4] = "ES_unt";
            solvernames[5] = "FIPS_A";
            solvernames[6] = "FIPS_B";
            solvernames[7] = "FIPS_unt";
            solvernames[8] = "GA";
            solvernames[9] = "PSO_A";
            solvernames[10] = "PSO_B";
            solvernames[11] = "PSO_unt";
            solvernames[12] = "SA";
            solvernames[13] = "SBplx";
            solvernames[14] = "Seye";
            solvernames[15] = "SGA_A";
            solvernames[16] = "SGA_B";
            solvernames[17] = "SGA_unt";
            solvernames[18] = "CMA_ES";
            solvernames[19] = "RBFOpt";

            List<double> f_CRS2 = new List<double>();
            List<double> f_Direct = new List<double>();
            List<double> f_ESn10A = new List<double>();
            List<double> f_ESn10B = new List<double>();
            List<double> f_ESunt = new List<double>();
            List<double> f_FIPSn10A = new List<double>();
            List<double> f_FIPSn10B = new List<double>();
            List<double> f_FIPSunt = new List<double>();
            List<double> f_GA = new List<double>();
            List<double> f_PSOn10A = new List<double>();
            List<double> f_PSOn10B = new List<double>();
            List<double> f_PSOunt = new List<double>();
            List<double> f_SA = new List<double>();
            List<double> f_SBplx = new List<double>();
            List<double> f_Seye = new List<double>();
            List<double> f_SGAn10A = new List<double>();
            List<double> f_SGAn10B = new List<double>();
            List<double> f_SGAunt = new List<double>();
            List<double> f_RBFopt = new List<double>();
            List<double> f_CMAES = new List<double>();

            List<List<string>> persolver = new List<List<string>>();


            int counter = 0;
            List<List<double>> all_fs = new List<List<double>>();
            foreach (string[] fileps in allPaths)
            {
                List<string> writelines = new List<string>();
                writelines.Add(solvernames[counter]);
                Console.WriteLine(solvernames[counter]);
                all_fs.Add(new List<double>());
                foreach (string filep in fileps)
                {
                    System.IO.StreamReader file = new System.IO.StreamReader(filep);
                    string line;
                    List<double> objs = new List<double>();

                    string[] text = new string[]{};
                    int linecount = 0;
                    while ((line = file.ReadLine()) != null && linecount < (intN + 1) * k)
                    {
                        text = line.Split(' ');
                        objs.Add(Convert.ToDouble(text[text.Length - 1]));
                        linecount++;
                    }
                    file.Close();

                    double min = objs.Min();
                    all_fs[counter].Add(min);
                    Console.WriteLine(min);
                    //string [] strleng = text[3].Split(',');
                    //Console.WriteLine(strleng.Length);
                    writelines.Add(min.ToString());
                }
                persolver.Add(writelines);
                counter++;
                Console.WriteLine();
                //writelines.Add("");
            }
            Console.WriteLine("DONE");
            Console.ReadKey();

            List<string> writeall = new List<string>();
            for (int i = 0; i < persolver[0].Count; i++ )
            {
                string line = "";
                for (int u = 0; u < persolver.Count; u++)
                {
                    line += persolver[u][i].ToString();
                    line += @";";
                }
                writeall.Add(line);
            }

            //write into text file
            System.IO.File.WriteAllLines(@"H:\PROJEKTE\16_OptimizerBenchmarking\4_CASESTUDY\Journal_Benchmarking\TEST_FUNCS\" + func + @"\" + n + @"\AllMins.txt", writeall.ToArray());
        }


        /// <summary>
        /// Running Sobol sequence (loaded from textfile) for mathematical test functions and saving output as textfile
        /// </summary>
        /// <param name="args"></param>
        static void sobolMain(string[] args)
        {
            int n = 20;
            string PathInput = @"H:\PROJEKTE\18_FitnessLandscapeAnalysis\04_CASESTUDY\TestFunctions\Sensitivity\N20_input.csv";

            List<Func<double[], double>> tf = new List<Func<double[], double>>();


            //n = 10 and n = 13
            tf.Add(MetaheuristicsTuner.Testfunctions.SO.B_Sphere);     // [-1, 1] for all x
            tf.Add(MetaheuristicsTuner.Testfunctions.SO.L_Ackley);     // [-32.768, 32.768] for all x
            tf.Add(MetaheuristicsTuner.Testfunctions.SO.L_Rastrigin);    // [-5.12, 5.12] for all x
            tf.Add(MetaheuristicsTuner.Testfunctions.SO.V_Rosenbrock); // [-2.048, 2.048] for all x
            double[] lb = new double[tf.Count];
            double[] ub = new double[tf.Count];
            lb[0] = -1;
            lb[1] = -32.768;
            lb[2] = -5.12;
            lb[3] = -2.048;
            ub[0] = 1;
            ub[1] = 32.768;
            ub[2] = 5.12;
            ub[3] = 2.048;


            //tf.Add(MetaheuristicsTuner.Testfunctions.SO.C_MichalewiczSchoenauer); // only for n=13. domain: [0,1] for i=0-8 and 12; [0,100] for i=9-11
            //double[] lb = new double[13];
            //double[] ub = new double[13];
            //for (int i = 0; i < 9; i++)
            //{
            //    lb[i] = 0;
            //    ub[i] = 1;
            //}
            //lb[12] = 0;
            //ub[12] = 1;
            //for (int i = 9; i < 12; i++)
            //{
            //    lb[i] = 0;
            //    ub[i] = 100;
            //}




            //load input sequence for sampling
            List<double[]> x = new List<double[]>();
            System.IO.StreamReader file = new System.IO.StreamReader(PathInput);
            string line;
            string[] text = new string[] { };
            int linecount = 0;
            while ((line = file.ReadLine()) != null)
            {
                x.Add(new double[n]);
                text = line.Split(';');
                int counter = 0;
                foreach(string t in text)
                {
                    x[linecount][counter] = Convert.ToDouble(t);
                    counter++;
                }
                linecount++;
            }
            file.Close();

            //compute outputs for x for every test function.
            List<double[]> ytf = new List<double[]>();
            int counterr = 0;
            foreach (Func<double[], double> f in tf)
            {
                ytf.Add(new double[x.Count]);
                for (int i = 0; i < x.Count; i++)
                {
                    double[] xin = new double[n];
                    for (int nn = 0; nn < n; nn++)
                    {
                        // NOT FOR MICHALEWICZ FUNCTION
                        xin[nn] = (x[i][nn] * (ub[counterr] - lb[counterr])) + lb[counterr];

                        //// ONLY MICHALEWICZ
                        //xin[nn] = (x[i][nn] * (ub[nn] - lb[nn])) + lb[nn];
                    }
                    
                    ytf[counterr][i] = f(xin);
                }
                counterr++;
            }


            //write outputs
            string PathOutput = @"H:\PROJEKTE\18_FitnessLandscapeAnalysis\04_CASESTUDY\TestFunctions\Sensitivity\";
            for (int fn = 0; fn < tf.Count; fn++)
            {
                List<string> writeall = new List<string>();
                for (int i = 0; i < ytf[fn].Length; i++)
                {
                    writeall.Add(ytf[fn][i].ToString());
                }
                
                
                //write into text file

                //ALL OTHER FUNCTIONS
                System.IO.File.WriteAllLines(PathOutput + "y" + fn + "_n" + n + ".txt", writeall.ToArray());

                ////MICHALEWICZ FUNCTION
                //System.IO.File.WriteAllLines(PathOutput + "y" + 4 + "_n" + n + ".txt", writeall.ToArray());

            }





            
        }



        #region Multi-Objective Tests
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
        static void MultiObjective_Main(string[] args)
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
        #endregion


        /// <summary>
        /// runtime analysis for test functions
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args)
        {
            int n = 20;

            double[] x = new double[n];
            for (int i = 0; i < n; i++) x[i] = 0;


            Func<double[], double>[] funcs = new Func<double[], double>[20];
            funcs[0] = SO.B_Perm0db;
            funcs[1] = SO.B_RotHypEll;
            funcs[2] = SO.B_Sphere;
            funcs[3] = SO.B_SumSquares;
            funcs[4] = SO.B_Trid;
            funcs[5] = SO.L_Ackley;
            funcs[6] = SO.L_Griewank;
            funcs[7] = SO.L_Levy;
            funcs[8] = SO.L_Rastrigin;
            funcs[9] = SO.L_Schwefel;
            funcs[10] = SO.O_PermDB;
            funcs[11] = SO.O_StyblinskiTang;
            funcs[12] = SO.P_Zakharov;
            funcs[13] = SO.V_DixonPrice;
            funcs[14] = SO.V_Rosenbrock;
            funcs[15] = SO.L_Levy;
            funcs[16] = SO.L_Schwefel;
            funcs[17] = SO.O_StyblinskiTang;
            funcs[18] = SO.V_DixonPrice;
            funcs[19] = SO.V_Rosenbrock;


            var watch = System.Diagnostics.Stopwatch.StartNew();
            // the code that you want to measure comes here

            //for (int u = 0; u < 5000; u++)
            //{
                for (int k = 0; k < 30; k++)
                {
                    for (int i = 0; i < (n + 1) * 100; i++)
                    {
                        foreach (Func<double[], double> fun in funcs)
                        {
                            double a = fun(x);
                        }
                    }
                }
            //}
            watch.Stop();
            Console.WriteLine(watch.ElapsedMilliseconds);
            Console.ReadKey();


            //around 85 ms for n=4
            //around 3670 ms for n=4

            //each times 5000 
            // times 20

            //n4 -> 7.08333 minutes x 20
            //n20 -> 305.83333 minutes x 20

            //for each solver PSO, FIPS, SGA, ES. ignoring overhead of each algorithm implementation (our ES implementaion is slower than the SGA)
        }


        #region Single-Objective Tests
        /// <summary>
        /// Testing a SO solver from MetaheuristicsLibrary.SovlersSO
        /// </summary>
        /// <param name="args"></param>
        static void hhMain(string[] args)
        {
           
            int seeds = 1;
            double[] optis = new double[seeds];
            int dvar = 2;
            int evalcount = (dvar) * 50;
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            bool[] xint = new bool[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = 0;
                ub[i] = 7;
                xint[i] = false;
            }
            lb[0] = -2.2;
            lb[1] = -0.75;
            ub[0] = 2.0;
            ub[1] = 2.5;


            Func<double[], double> testfunc = SO.V_Camel3Hump;


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
            //int simplegapop = Convert.ToInt32(6);   //x[0]
            //simplegapop += simplegapop % 2;
            //settings.Add("maxgen", Convert.ToInt32(Math.Floor(evalcount / Convert.ToDouble(simplegapop))));
            //settings.Add("popsize", simplegapop);
            //settings.Add("k", 35);               //x[1]  6
            //settings.Add("pcross", 1);          //x[2]  0.7
            //settings.Add("pmut", 1);          //x[3]  0.3
            //settings.Add("d", 0.2);             //x[4]  0.1
            //settings.Add("r", 2);             //x[5]  0.1
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



            //Dictionary<string, object> settingsES = new Dictionary<string, object>();
            //settingsES.Add("popsize", 3);          // ∈ {2,...,200}
            //settingsES.Add("lambda", 1);            // ∈ {1,...,200}
            //settingsES.Add("roh", 2);               // ∈ {1,...,popsize}  . in hyperoptimization, express as percentage of lambda
            //settingsES.Add("x0sampling", 0);        // ∈ {0,1}  0=uniform, 1=gaussian
            //settingsES.Add("stepsize0", 10);       // ∈ [0.01, 10]
            //settingsES.Add("stepsize", 1.8);        // ∈ [0.01, 10]
            //settingsES.Add("tauc", 0.01);              // ∈ [0.01, 50]
            //settingsES.Add("selmode", 1);
            ////settingsES.Add("pmut_int", 0.1);        // ∈ [0.01, 0.99] 
            //SimpleES[] es = new SimpleES[seeds];
            //for (int i = 0; i < seeds; i++)
            //{
            //    es[i] = new SimpleES(lb, ub, xint, evalcount, testfunc, i, settingsES);
            //    es[i].solve();
            //    optis[i] = es[i].get_fxoptimum();
            //}
            //Console.WriteLine("es average: {0}", optis.Average());
            //Console.ReadKey();





            //Dictionary<string, object> settingsPSO = new Dictionary<string, object>();
            //settingsPSO.Add("popsize", 10);        // popsize                           ∈ {4,..., 100}
            //settingsPSO.Add("chi", 0.39852);            // constriction coefficient         ∈ [0.001, 1]
            //settingsPSO.Add("phi", 9.24731);              // attraction to best particle     ∈ [0.01, 50]
            //settingsPSO.Add("v0max", 0.01699);          // max velocity at initialisation. fraction of domain. ∈ [0.01, 10]
            //settingsPSO.Add("x0samplingmode", 0);   // 0 = uniform, 1 = gaussian        ∈ [0.01, 10]
            //settingsPSO.Add("pxupdatemode", 0);     //0 = update after population. 1 = update after each evaluation ∈ [0.01, 10]
            //settingsPSO.Add("s0", 1);             //initial step size in case of gaussian x0 ∈ [0.01, 10]
            //settingsPSO.Add("psomode", 0);      //0 = fipso, 1 = inertia, 2 = constriction
            //settingsPSO.Add("phi1", 0);      //attraction own best
            //settingsPSO.Add("phi2", 0);      //attraction global best
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





            //Dictionary<string, object> settingsNM = new Dictionary<string, object>();
            //settingsNM.Add("alpha", 1);
            //settingsNM.Add("gamma", 1.5);
            //settingsNM.Add("rho", 0.25);
            //settingsNM.Add("sigma", 0.2);
            //settingsNM.Add("step0", 0.02);
            //NelderMead[] nm = new NelderMead[seeds];
            //for (int i = 0; i < seeds; i++)
            //{
            //    nm[i] = new NelderMead(lb, ub, xint, evalcount, testfunc, i, settingsNM);
            //    nm[i].solve();
            //    optis[i] = nm[i].get_fxoptimum();
            //}
            //Console.WriteLine("nm average: {0}", optis.Average());
            //Console.WriteLine("nm min: {0}", optis.Min());
            //Console.WriteLine("nm max: {0}", optis.Max());
            //Console.ReadKey();



            //double[] x0 = new double[dvar];
            //x0[0] = 1.5;
            //x0[1] = 2;

            //Dictionary<string, object> settingsRB = new Dictionary<string, object>();
            //settingsRB.Add("alpha", 3);
            //settingsRB.Add("beta", 0.5);
            //settingsRB.Add("stepsize", 0.125);
            //Rosenbrock[] rb = new Rosenbrock[seeds];
            //for (int i = 0; i < seeds; i++)
            //{
            //    rb[i] = new Rosenbrock(lb, ub, xint, evalcount, testfunc, i, settingsRB, x0);
            //    rb[i].solve();
            //    optis[i] = rb[i].get_fxoptimum();
            //}
            //Console.WriteLine("rb average: {0}", optis.Average());
            //Console.WriteLine("rb min: {0}", optis.Min());
            //Console.WriteLine("rb max: {0}", optis.Max());
            //Console.ReadKey();





            Direct dr = new Direct(lb,ub, xint, evalcount, testfunc, 0);
            dr.solve();
            Console.WriteLine("direct min: {0}", dr.get_fxoptimum());
            Console.ReadKey();

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
        #endregion

    }


}
