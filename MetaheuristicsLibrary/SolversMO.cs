using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MetaheuristicsLibrary.Misc;
namespace MetaheuristicsLibrary.SolversMO
{

    /// <summary>
    /// SPEA2
    ///  <para/>Strength Pareto Evolutionary Algorithm 2.
    ///  <para/>Paper: Zitzler E. et al. (2001). SPEA2: Improving the Strength Pareto Evolutionary Algorithm. In: Evolutionary Methods for Design Optimization and Control with Applications to Industrial Problems 2001: 95-100.
    ///  <para/>Source code: http://yarpiz.com/74/ypea122-spea2
    ///  <para/>Description:
    /// </summary>
    /// <remarks>
    /// ???
    /// </remarks>
    public class SPEA2
    {
        public bool flagXNormalized { get; protected set; }

        public int maxFncCalls { get; set; }
        public int FncCalls { get; set; }
        public double expectedBest { get; set; }
        public bool terminateNow { get; set; }

        public int nVar { get; protected set; }
        public double[] ub { get; protected set; }
        public double[] lb { get; protected set; }
        public bool[] integers { get; protected set; }

        public int rndSeed { get; protected set; }
        public RandomDistributions rnd { get; set; }

        public double[] x0 { get; protected set; }                          //x(0)
        public List<double[]> x0pop { get; protected set; }                       //x(0) population
        public double[] x { get; protected set; }                           //x(t)
        public double[] xBest { get; protected set; }                       //x(best)
        public double[] xArchive { get; protected set; }                    //x(archive). doesnt have to be the best
        public double[] xTest { get; protected set; }                       //x temporary

        public int PopSize { get; set; }                                    //population size
        protected int xTestCounter { get; set; }                            //counter for how many xTest are created within one iteration t
        public List<double[]> xPopulation0 { get; protected set; }          //xi(0)
        public List<double[]> xPopulation { get; protected set; }           //xi(t)
        public List<double[]> xPopulationBest { get; protected set; }       //xi(best)
        public List<double[]> xPopulationTest { get; protected set; }       //temporary x population 
        public List<double[]> xPopulationArchive { get; set; }              //archive. not trajectory based, as xPopulationBest

        public double objval0 { get; protected set; }                       //f(x0)
        public double objval { get; protected set; }                        //f(x)
        public double objvalBest { get; protected set; }                    //f(x(best))
        public double objvalArchive { get; protected set; }                 //f(x(archive))
        public double objvalTest { get; set; }                              //temporary f(x)
        public List<double> objvalPopulation0 { get; protected set; }       //f(xi0)
        public List<double> objvalPopulation { get; protected set; }        //f(xi))
        public List<double> objvalPopulationBest { get; set; }              //f(xi(best))
        public List<double> objvalPopulationTest { get; set; }              //temporary f(xi)
        public List<double> objvalPopulationArchive { get; protected set; } //f(xi(archive))

        public int mObj;                                      //number of objectives
        public double[] objvalTest_MO;                                //array of obj func values for each objective
        public List<double[]> objvalsTestPop_MO { get; set; }   //obj vals for each objective and solution candidate
        public List<double[]> objvalsBestPop_MO { get; set; }   //best obj vals for each obj. and best solution candidate

        public int?[] Rank { get; set; }                         //rank of the set. number indicating pareto fronts
        public int[] Strength { get; set; }                    //Strength of non domination. how many individuals does each individual i dominate
        public int[] RawFitness { get; set; }                    //Raw Fitness. sum of strengths of indivduals, which dominate individual i
        public double[] CrowdingDist { get; set; }              //crowding distance of the set
        public double[] KNNDensity { get; set; }                // Density, according to k-th nearest neighbour
        List<List<int>> Fronts;                                 // Pareto Fronts

        private Func<double[], double[]> evaluationFunction;       //obejctive function, returning multiple objectives

        //SPEA2 Main Loop (from Zitzler et al. 2001)
        //
        // Input:       N           population size
        //              Na          archive size
        //              T           maximum number of generations
        //
        // Output:      A           nondominated set
        //
        // Step 1:      Initialization:
        //              Generate an initial population Po and create the empty archive (external set) Pa0=∅. Set t=0.
        //
        // Step 2:      Fitness assignment:
        //              Calculate fitness values of individuals in Pt and Pat
        //
        // Step 3:      Environmental selection:
        //              copy all nondominated individuals in Pt and Pat to Pt+1. If size of Pat+1 exceeds Na then reduce Pt+1 by means of the truncation operator, 
        //              otherwise if size of Pt+1 is less than Na then fill Pt+1 with dominated individuals in Pt and Pat.
        //
        // Step 4:      Termination:
        //              If t>=T or another stopping criterion is satisfied then set A to the set of decision vectors represented by the nondominated individuals in Pat+1. Stop.
        //
        // Step 5:      Mating selection:
        //              Perform binary tournament selection with replacement on Pat1+ in order to fill the mating pool.
        //
        // Step 6:      Variation:
        //              Apply recombination and mutation operators to the mating pool and set Pt+1 to the resulting population. Increment generation counter (t=t+1) and go to Step 2.


        public int ga_nPop;                                     // population size
        protected int ga_nArchive;                              // archive size
        public int ga_gennow { get; protected set; }            // current gen
        public int ga_genmax;                                   // max generations

        public double ga_PCrossover;                            // probability for crossover
        protected int ga_nCrossover;                            // number of offspring generated by crossover
        protected int ga_nPairs;                                // number of pairs, which generate offspring
        public double ga_PMutation;                             // probability for mutation
        protected int ga_nMutation;                             // number of offspring generated by mutation

        public int[] ga_pairs { get; protected set; }           // pair of parents for recombination

        public int spea2_K { get; protected set; }              // KNN Parameter                                                K       = sqrt(nPop+nArchive)
        public int[] spea2_S { get; protected set; }            // strength value                                               S(i)    = |{j|j∈Pt+Parchive∧i≻j}|
        public int[] spea2_R { get; protected set; }            // raw fitness                                                  R(i)    = Sum(Sj) |j∈Pt+Parchive∧j≻i
        public double[] spea2_D { get; protected set; }         // density (K-Nearest Neighbours)                               D(i)    = 1 / (σk(i)+2)
        public double[,] spea2_sigma { get; protected set; }   // for each individual, distances to all other individuals      σ(i)    : vector
        public double[] spea2_sigmaK { get; protected set; }    // distances to kth-neighbour                                   σk(i)   = σ(i,K)

        public double M_InterRecomb_d;                          // Intermediate Recombination, extension of hypercube.          d∈[0,1]
        public double M_Mutation_r;                             // Mutation range. 0.1 means 10% of variable domain.            r∈[0,1]
        public int M_Mutation_k;                                // Mutation precision. The higher, the less the spread          k>=1, Integer


        //stuff for integer problem
        protected int ga_int_nOffspring;                                // number of offspring
        protected int ga_int_nParents;                                  // number of parents
        public int ga_int_toursize;                                 // tournament size
        public int ga_int_mu;                                       // distribution index for crossover
        public int ga_int_mum;                                      // distribution index for mutation
        public int[] ga_int_pairs { get; protected set; }          // pair of parents for recombination




        public SPEA2(int mObj, int nVar, double[] lb, double[] ub, Func<double[], double[]> evaluationFunction, int rndSeed, List<double[]> x0pop = null, bool[] integers = null)
        {
            if (x0pop != null)
            {
                this.x0pop = new List<double[]>(x0pop);
            }
            if (integers != null)
            {
                this.integers = new bool[nVar];
                integers.CopyTo(this.integers, 0);
            }


            this.rndSeed = rndSeed;
            this.rnd = new RandomDistributions(this.rndSeed);
            this.nVar = nVar;
            this.lb = new double[nVar];
            this.ub = new double[nVar];
            Array.Copy(lb, this.lb, nVar);
            Array.Copy(ub, this.ub, nVar);
            this.evaluationFunction = evaluationFunction;

            this.mObj = mObj;



            this.expectedBest = 0;
            this.maxFncCalls = 100;

            this.FncCalls = 0;

            this.xBest = new double[nVar];                  //need this for most single objective operators
            this.xTest = new double[nVar];                  // same
            this.xPopulationTest = new List<double[]>();
            this.xPopulationArchive = new List<double[]>();
            this.objvalsTestPop_MO = new List<double[]>();
            this.objvalsBestPop_MO = new List<double[]>();

            this.terminateNow = false;





            ga_nPop = 50;
            ga_gennow = 0;
            ga_genmax = 200;


            M_InterRecomb_d = 0.1;
            M_Mutation_r = 0.2;
            M_Mutation_k = 10;


            ga_int_toursize = 2;
            ga_int_mu = 20;
            ga_int_mum = 20;

            ga_PCrossover = 0.7;
            ga_PMutation = 1 - ga_PCrossover;
        }

        /// <summary>
        /// Initialize the solver.
        /// This includes evaluating the initial solutions
        /// </summary>
        public void initialize()
        {
            ga_int_nOffspring = ga_nPop / 2;
            ga_int_nParents = ga_int_nOffspring;


            ga_nArchive = ga_nPop;

            ga_nPairs = (int)((ga_PCrossover * ga_nPop) / 2);
            ga_nCrossover = ga_nPairs * 2;

            ga_nMutation = ga_nPop - ga_nCrossover;

            spea2_K = (int)(Math.Sqrt(ga_nPop + ga_nArchive));
            spea2_S = new int[ga_nPop + ga_nArchive];
            spea2_R = new int[ga_nPop + ga_nArchive];
            spea2_D = new double[ga_nPop + ga_nArchive];
            spea2_sigma = new double[ga_nPop + ga_nArchive, ga_nPop + ga_nArchive];
            spea2_sigmaK = new double[ga_nPop + ga_nArchive];


            I_uniformRandomPop(ga_nPop);
            E_MO_evaluateTestPopulation();


            //non domination ranks: Rank
            S_MO_NonDominatedFront();

            //calculate strength S(i): Strength
            S_MO_NonDomStrength();

            //raw fitness R(i): RawFitness
            S_MO_NonDomRawFitness();

            //Density D(i)
            S_MO_KNNDensity(spea2_K);

            //spea2 fitness F(i)
            objvalPopulationTest = new List<double>();
            for (int i = 0; i < xPopulationTest.Count; i++)
            {
                objvalPopulationTest.Add(RawFitness[i] + KNNDensity[i]);
            }



            //first population is also current best
            xPopulationArchive = new List<double[]>(xPopulationTest);
            objvalsBestPop_MO = new List<double[]>(objvalsTestPop_MO);
            objvalPopulationBest = new List<double>(objvalPopulationTest);
        }

        /// <summary>
        /// Run the solver until termination criterion is fulfilled (max generations).
        /// </summary>
        public void Solve()
        {
            while (this.terminateNow == false)
            {
                move();            //should maybe be called "explore"? and then another one "localsearch", or "exploit"...... or just .search()... explore and exploit are parts of .search
                evaluate();          //this could include function approximations. RBF, Kriging, ANN, ...
                collect();              //.best should be part of select?...
                select();            // this would include sorting?...
                terminate();          // could contain RESTART operator
            }
        }

        /// <summary>
        /// Run the solver for a specific number of iterations.
        /// </summary>
        /// <param name="iterations">number or iterations.</param>
        public void Solve(int iterations)
        {
            for (int i = 0; i < iterations; i++)
            {
                if (this.terminateNow == false)
                {
                    move();
                    evaluate();
                    collect();
                    select();
                    terminate();
                }
                else
                {
                    break;
                }
            }
        }

        /// <summary>
        /// Get the name of the instantiated solver.
        /// </summary>
        /// <returns>String</returns>
        public string GetSolverName() { return this.GetType().Name; }




        /// <summary>
        /// Generate new candidate solutions using recombination and mutation.
        /// </summary>
        private void move()
        {
            M_InterRecomb(S_BinTournSelPairs(ga_nPairs), M_InterRecomb_d);
            M_MutationPohlheim(S_BinTournSel(ga_nMutation), M_Mutation_r, M_Mutation_k, false);

            ga_gennow++;
        }

        /// <summary>
        /// Evaluate the candidate solutions with the evaluation function.
        /// </summary>
        private void evaluate()
        {
            //F(X)... ->objvalsTestPop_MO
            E_MO_evaluateTestPopulation();
        }

        /// <summary>
        /// Collect solutions into an archive. Also, the multi-objective quality indicators are calculated here.
        /// </summary>
        private void collect()
        {
            //fill new xBestPopulation (Archive).
            //      sort according to PF(i)... objvalsTestPopulation
            //      then go through one by one, but only take those with pareto dominance S(i) on first rank (non-dominated). 
            //              if there are more than nArchive non-dominated, do archive truncation. (p.8)
            //              if its not full (while <nArchive), go to next rank, next rank, and so on.



            //****************************************************************************
            //unify test and best pop:
            xPopulationTest.AddRange(xPopulationArchive);
            objvalsTestPop_MO.AddRange(objvalsBestPop_MO);

            //pareto dominance S(i)
            S_MO_NonDominatedFront();
            S_MO_NonDomStrength();

            //raw fitness R(i)
            S_MO_NonDomRawFitness();

            //sigma
            //sigmaK
            //density D(i)    -> propose this sum thing
            S_MO_KNNDensity(spea2_K);


            //pareto Fitness PF(i)... ->objvalsTestPopulation  (from single objective)
            objvalPopulationTest = new List<double>();
            for (int i = 0; i < xPopulationTest.Count; i++)
            {
                objvalPopulationTest.Add(RawFitness[i] + KNNDensity[i]);
            }

        }

        /// <summary>
        /// Selecting solutions for the next iteration.
        /// </summary>
        private void select()
        {
            //environmental selection
            //sort according to F(i) and truncate
            B_MO_ArchiveTruncation(ga_nArchive, false);
        }

        /// <summary>
        /// Checking termination criterion (e.g. maximum generations)
        /// </summary>
        private void terminate()
        {
            T_MaxAlgoSpecificIter(ga_gennow, ga_genmax);
        }









        /// <summary>
        /// Creates a random initial population with uniform random distribution.
        /// </summary>
        /// <param name="popsize"></param>
        private void I_uniformRandomPop(int popsize)
        {
            this.xPopulationTest = new List<double[]>();

            int iStart = 0;
            if (this.x0 != null)
            {
                this.xPopulationTest.Add(this.x0);
                iStart++;
            }
            else if (this.x0pop != null)
            {
                for (int i = 0; i < x0pop.Count && i < popsize; i++)
                {
                    this.xPopulationTest.Add(this.x0pop[i]);
                    iStart++;
                }
            }

            for (int i = iStart; i < popsize; i++)
            {
                I_uniformRandom();
                double[] copy = new double[nVar];
                Array.Copy(xTest, copy, nVar);
                xPopulationTest.Add(copy);
            }
        }

        /// <summary>
        /// Creates a random solution with uniform random distribution.
        /// </summary>
        private void I_uniformRandom()
        {
            if (this.flagXNormalized)
            {
                for (int i = 0; i < this.nVar; i++)
                {
                    this.xTest[i] = this.rnd.NextDouble();
                }
            }
            else
            {
                for (int i = 0; i < this.nVar; i++)
                {
                    this.xTest[i] = this.rnd.NextDouble() * (this.ub[i] - this.lb[i]) + this.lb[i];
                    if (this.integers != null)
                    {
                        if (this.integers[i]) this.xTest[i] = Math.Round(this.xTest[i], 0);
                    }
                }
            }
        }


        /// <summary>
        /// Evaluate f(x) for all solution candidates in solver.xTestPopulation.
        /// </summary>
        private void E_MO_evaluateTestPopulation()
        {
            this.objvalsTestPop_MO = new List<double[]>();

            foreach (double[] xTest in this.xPopulationTest)
            {
                this.objvalTest_MO = this.evaluationFunction(xTest);
                this.objvalsTestPop_MO.Add(this.objvalTest_MO);
                this.FncCalls++;
            }
        }


        /// <summary>
        /// Intermediate Recombiation.
        /// <para/>Creates from each parents pair 1 offspring. New variables constructed on a (extended) hypercube between parents variables.
        /// <para/>working on xBestPopulation.
        /// </summary>
        /// <param name="parents">pair of parents. indicating index in the Population array.</param>
        /// <param name="d">if d>0, then hypercube is extended. if d=0 then hypercube not extended.</param>
        private void M_InterRecomb(int[][] parents, double d)
        {
            xPopulationTest = new List<double[]>();

            //pair of parents have to be chosen before, e.g. BinaryTournamentSelection
            int nParents = parents.Length;
            for (int i = 0; i < nParents; i++)
            {
                double[] child_1 = new double[nVar];
                double[] child_2 = new double[nVar];
                double[] xPar_1 = xPopulationArchive[parents[i][0]].ToArray();
                double[] xPar_2 = xPopulationArchive[parents[i][1]].ToArray();
                for (int n = 0; n < nVar; n++)
                {

                    //if (this.integers != null)
                    //{
                    if (this.integers[n])
                    {
                        double[] children = M__BinCrossover(xPar_1[n], xPar_2[n], this.ub[n], this.lb[n]);
                        child_1[n] = children[0];
                        child_2[n] = children[1];
                    }
                    else
                    {
                        double alpha = rnd.NextDouble() * ((1 + d) + d) - d;
                        child_1[n] = xPar_1[n] * alpha + xPar_2[n] * (1 - alpha);
                        child_2[n] = xPar_2[n] * alpha + xPar_1[n] * (1 - alpha);
                    }
                    //}
                    if (child_1[n] < lb[n]) child_1[n] = lb[n];
                    if (child_2[n] < lb[n]) child_2[n] = lb[n];
                    if (child_1[n] > ub[n]) child_1[n] = ub[n];
                    if (child_2[n] > ub[n]) child_2[n] = ub[n];
                }
                xPopulationTest.Add(child_1);
                xPopulationTest.Add(child_2);
            }
        }

        /// <summary>
        /// Single Point Crossover of Bitstrings.
        /// </summary>
        /// <param name="xpar1"></param>
        /// <param name="xpar2"></param>
        /// <param name="ub"></param>
        /// <param name="lb"></param>
        /// <returns></returns>
        private double[] M__BinCrossover(double xpar1, double xpar2, double ub, double lb)
        {
            bool[] intx = new bool[] { true };
            double[] _lb = new double[] { lb };
            double[] _ub = new double[] { ub };
            int reqBitLen = Misc.Misc.BinReqLength(intx, _lb, _ub);
            if (reqBitLen < 2) reqBitLen = 2;

            string binx1 = Misc.Misc.Dec2Bin(Convert.ToInt32(xpar1 - lb), reqBitLen);
            string binx2 = Misc.Misc.Dec2Bin(Convert.ToInt32(xpar2 - lb), reqBitLen);
            char[] binxarr1 = binx1.ToCharArray();
            char[] binxarr2 = binx2.ToCharArray();
            int pointer;
            if (reqBitLen > 2)
            {
                pointer = Convert.ToInt32(rnd.NextDouble() * ((double)reqBitLen - 2)) + 1;
            }
            else
            {
                pointer = 1;
            }


            char[] newbinxarr1 = new char[reqBitLen];
            char[] newbinxarr2 = new char[reqBitLen];
            for (int i = 0; i < reqBitLen; i++)
            {
                if (i < pointer)
                {
                    newbinxarr1[i] = binxarr1[i];
                    newbinxarr2[i] = binxarr2[i];
                }
                else
                {
                    newbinxarr1[i] = binxarr2[i];
                    newbinxarr2[i] = binxarr1[i];
                }
            }

            int xdual1 = Misc.Misc.Bin2Dec(new string(newbinxarr1));
            int xdual2 = Misc.Misc.Bin2Dec(new string(newbinxarr2));

            double xdbl1 = Convert.ToDouble(xdual1) + lb;
            double xdbl2 = Convert.ToDouble(xdual2) + lb;
            double rndmutate = rnd.NextDouble();
            if (rndmutate < 0.2)
            {
                xdbl1 = Convert.ToInt32(rnd.NextDouble() * (ub - lb) + lb);
                xdbl2 = Convert.ToInt32(rnd.NextDouble() * (ub - lb) + lb);
            }

            return new double[] { xdbl1, xdbl2 };
        }

        /// <summary>
        /// Mutation of real numbers.
        /// <para/>working on xBestPopulation.
        /// <para/>Pohlheim (1999). Evolutionäre Algorithmen.
        /// </summary>
        /// <param name="x2mutate">indices of xBestPopulation to be mutated</param>
        /// <param name="r">mutation range. r∈[0,1]. 0.1 means 10% of the variable domain.</param>
        /// <param name="k">mutation precision. Higher number means less spread.</param>
        /// <param name="eraseTest">erase xTestPopulataion and start new list now?</param>
        private void M_MutationPohlheim(int[] x2mutate, double r, int k, bool eraseTest)
        {
            if (eraseTest) xPopulationTest = new List<double[]>();

            int nMutate = x2mutate.Length;
            for (int i = 0; i < nMutate; i++)
            {
                double[] xMutated = new double[nVar];
                for (int n = 0; n < nVar; n++)
                {
                    if (this.integers[n])
                    {
                        xMutated[n] = M__BinMutateUnif(xPopulationArchive[x2mutate[i]][n], r + 0.1, this.ub[n], this.lb[n]);
                    }
                    else
                    {
                        xMutated[n] = xPopulationArchive[x2mutate[i]][n] + (rnd.NextDouble() * 2 - 1) * (r * (Math.Abs(ub[n] - lb[n]))) * (Math.Pow(2, (rnd.NextDouble() * -1) * k));
                    }
                    //if (this.integers != null)
                    //{

                    //}
                    if (xMutated[n] < lb[n]) xMutated[n] = lb[n];
                    if (xMutated[n] > ub[n]) xMutated[n] = ub[n];
                }
                xPopulationTest.Add(xMutated);
            }

        }

        private double M__BinMutateUnif(double x2mutate, double prob, double ub, double lb)
        {
            int xMutated;
            bool[] intx = new bool[] { true };
            double[] _lb = new double[] { lb };
            double[] _ub = new double[] { ub };
            int reqBitLen = Misc.Misc.BinReqLength(intx, _lb, _ub);
            if (reqBitLen < 2) reqBitLen = 2;

            string binx = Misc.Misc.Dec2Bin(Convert.ToInt32(x2mutate - lb), reqBitLen);
            char[] binxArray = binx.ToCharArray();
            for (int i = 0; i < binxArray.Length; i++)
            {
                double rn = rnd.NextDouble();
                if (rn <= prob)
                {
                    if (binxArray[i].Equals('1'))
                    {
                        binxArray[i] = '0';
                    }
                    else
                    {
                        binxArray[i] = '1';
                    }
                }
            }
            xMutated = Misc.Misc.Bin2Dec(new string(binxArray));
            return (Convert.ToDouble(xMutated) + lb);
        }

        /// <summary>
        /// Checking termination criterion. For example maximum Generations in a Genetic Algorithm.
        /// </summary>
        private void T_MaxAlgoSpecificIter(int currIter, int maxIter)
        {
            if (currIter >= maxIter)
            {
                this.terminateNow = true;
            }
            else
            {
                this.terminateNow = false;
            }
        }

        /// <summary>
        /// Binary Tournament Selection.
        /// <para/>Creating pairs. 
        /// </summary>
        /// <param name="pairs">amount of pairs to be drawn</param>
        /// <returns>pairs of parents, drawn from xBestPopulation</returns>
        private int[][] S_BinTournSelPairs(int pairs)
        {
            int[][] couples = new int[pairs][];

            int[] parents = new int[pairs];
            for (int i = 0; i < pairs; i++)
            {
                couples[i] = new int[2];
                int drawn1 = 0; //making sure, the couple won't consist of the same individuals
                int drawn2 = 0;
                for (int u = 0; u < 2; u++)
                {
                    int i1 = rnd.Next(xPopulationArchive.Count);
                    if (u == 0)
                    {
                        drawn1 = i1;
                    }
                    else
                    {
                        while (drawn1 == i1 || drawn2 == i1) i1 = rnd.Next(xPopulationArchive.Count);
                    }

                    int i2 = rnd.Next(xPopulationArchive.Count);
                    while (i2 == i1) i2 = rnd.Next(xPopulationArchive.Count);      //making sure, you don't compete against yourself
                    if (u == 0)
                    {
                        drawn2 = i2;
                    }
                    else
                    {
                        while (drawn2 == i2 || drawn1 == i2) i2 = rnd.Next(xPopulationArchive.Count);
                    }

                    if (objvalPopulationBest[i1] < objvalPopulationBest[i2])
                    {
                        couples[i][u] = i1;
                    }
                    else
                    {
                        couples[i][u] = i2;
                    }
                }
            }


            return couples;
        }

        /// <summary>
        /// Binary Tournamet Selection.
        /// </summary>
        /// <param name="nSelect">amount of candidates to be selected</param>
        /// <returns>indices of xBestPopulation</returns>
        private int[] S_BinTournSel(int nSelect)
        {
            int[] chosen = new int[nSelect];
            for (int i = 0; i < nSelect; i++)
            {

                int i1 = rnd.Next(xPopulationArchive.Count);
                int i2 = rnd.Next(xPopulationArchive.Count);
                while (i2 == i1) i2 = rnd.Next(xPopulationArchive.Count);      //making sure, you don't compete against yourself


                if (objvalPopulationBest[i1] < objvalPopulationBest[i2])
                {
                    chosen[i] = i1;
                }
                else
                {
                    chosen[i] = i2;
                }
            }

            return chosen;
        }







        #region MOoperators
        /// <summary>
        /// returns this.Rank and sorts xTestPopulation, objvalsTestPop_MO
        /// source: https://www.mathworks.com/matlabcentral/fileexchange/10429-nsga-ii--a-multi-objective-optimization-algorithm
        /// </summary>
        private void S_MO_NonDominatedFront()
        {
            /*  Non-Dominated sort. 
            The initialized population is sorted based on non-domination. The fast
            sort algorithm [1] is described as below for each

            • for each individual p in main population P do the following
                – Initialize Sp = []. This set would contain all the individuals that is
                    being dominated by p.
                – Initialize np = 0. This would be the number of individuals that domi-
                    nate p.
                – for each individual q in P
                    * if p dominated q then
                        · add q to the set Sp i.e. Sp = Sp ? {q}
                    * else if q dominates p then
                        · increment the domination counter for p i.e. np = np + 1
                – if np = 0 i.e. no individuals dominate p then p belongs to the first
                    front; Set rank of individual p to one i.e prank = 1. Update the first
                    front set by adding p to front one i.e F1 = F1 ? {p}
            • This is carried out for all the individuals in main population P.
            • Initialize the front counter to one. i = 1
            • following is carried out while the ith front is nonempty i.e. Fi != []
                – Q = []. The set for storing the individuals for (i + 1)th front.
                – for each individual p in front Fi
                    * for each individual q in Sp (Sp is the set of individuals
                      dominated by p)
                        · nq = nq?1, decrement the domination count for individual q.
                        · if nq = 0 then none of the individuals in the subsequent
                          fronts would dominate q. Hence set qrank = i + 1. Update
                          the set Q with individual q i.e. Q = Q ? q.
                – Increment the front counter by one.
                – Now the set Q is the next front and hence Fi = Q.
    
            This algorithm is better than the original NSGA ([2]) since it utilize
            the informatoion about the set that an individual dominate (Sp) and
            number of individuals that dominate the individual (np).*/

            this.Rank = new int?[xPopulationTest.Count];

            int[] dominated_amount = new int[xPopulationTest.Count];     //number of individuals that dominate this individual
            List<int>[] dominate_who = new List<int>[xPopulationTest.Count];     //indices of individuals, this individual dominates

            int index_front = 0;
            Fronts = new List<List<int>>();

            int[] idx = new int[xPopulationTest.Count];

            Fronts.Add(new List<int>());
            for (int i = 0; i < xPopulationTest.Count; i++)
            {
                idx[i] = i;         //storing indices for sorting later on

                dominated_amount[i] = 0;
                dominate_who[i] = new List<int>();

                for (int j = 0; j < xPopulationTest.Count; j++)
                {
                    int dom_less = 0;
                    int dom_equal = 0;
                    int dom_more = 0;
                    for (int m = 0; m < mObj; m++)
                    {
                        if (objvalsTestPop_MO[i][m] < objvalsTestPop_MO[j][m])
                        {
                            dom_less++;
                        }
                        else if (objvalsTestPop_MO[i][m] == objvalsTestPop_MO[j][m])
                        {
                            dom_equal++;
                        }
                        else
                        {
                            dom_more++;
                        }
                    }
                    if (dom_less == 0 && dom_equal != mObj)
                    {
                        dominated_amount[i]++;
                    }
                    else if (dom_more == 0 && dom_equal != mObj)
                    {
                        dominate_who[i].Add(j);
                    }
                }
                if (dominated_amount[i] == 0)
                {
                    Rank[i] = 0;
                    Fronts[index_front].Add(i);     //individuals on the first pareto front (rank 0)
                }
                else
                {
                    Rank[i] = null;
                }
            }



            //find the subsequent fronts
            while (Fronts[index_front].Count > 0)
            {
                List<int> Q = new List<int>();
                for (int i = 0; i < Fronts[index_front].Count; i++)
                {
                    if (dominate_who[Fronts[index_front][i]].Count > 0)      //if this individual dominates other
                    {
                        for (int j = 0; j < dominate_who[Fronts[index_front][i]].Count; j++)
                        {
                            dominated_amount[dominate_who[Fronts[index_front][i]][j]]--;
                            if (dominated_amount[dominate_who[Fronts[index_front][i]][j]] == 0)
                            {
                                Rank[dominate_who[Fronts[index_front][i]][j]] = index_front + 1;
                                Q.Add(dominate_who[Fronts[index_front][i]][j]);
                            }
                        }
                    }
                }
                index_front++;
                Fronts.Add(Q);
            }



            Array.Sort(Rank, idx);      //sorting indices according to parate fronts

            //sort :
            //      xTestPopulation
            //      xBestPopulation
            //      objvalsTestPopulation
            //double[][] copyXBest = xBestPopulation.ToArray();
            double[][] copyXTest = xPopulationTest.ToArray();
            double[][] copyObjVals = objvalsTestPop_MO.ToArray();

            xPopulationTest = new List<double[]>();
            //xBestPopulation = new List<double[]>();
            objvalsTestPop_MO = new List<double[]>();

            for (int i = 0; i < idx.Length; i++)
            {
                xPopulationTest.Add(copyXTest[idx[i]]);
                //xBestPopulation.Add(copyXTest[idx[i]]);
                objvalsTestPop_MO.Add(copyObjVals[idx[i]]);
            }

        }

        /// <summary>
        /// S(i)    = |{j|j∈Pt+Parchive∧i≻j}|
        /// <para/>returns Strength for each individual in xTestPopulation 
        /// <para/>i.e. number of individuals it dominates
        /// <para/>requires base.Rank[] to be calculated already
        /// </summary>
        private void S_MO_NonDomStrength()
        {
            int fronts = (int)Rank.Last();
            Strength = new int[xPopulationTest.Count];
            int[] counts = new int[fronts + 1];
            int u = 0;
            for (int i = 0; i < fronts + 1; i++)
            {
                do
                {
                    if ((int)Rank[u] == i) counts[i]++;
                    else break;
                    u++;
                } while (u < xPopulationTest.Count);
            }
            for (int i = 0; i < Strength.Length; i++)
            {
                int sum = 0;
                for (u = (int)Rank[i] + 1; u < counts.Length; u++)
                {
                    sum += counts[u];
                }
                Strength[i] = sum;
            }
        }

        /// <summary>
        /// R(i)    = Sum(Sj) |j∈Pt+Parchive∧j≻i
        /// <para/>returns RawFitness for each individual in xTestPopulation.
        /// <para/>i.e. sum of S(j)
        /// <para/>requires base.Rank and base.Strength
        /// </summary>
        private void S_MO_NonDomRawFitness()
        {
            int fronts = (int)Rank.Last();
            RawFitness = new int[xPopulationTest.Count];
            int[] counts = new int[fronts + 1];

            int u = 0;
            for (int i = 0; i < fronts + 1; i++)
            {
                do
                {
                    if ((int)Rank[u] == i) counts[i]++;
                    else break;
                    u++;
                } while (u < xPopulationTest.Count);
            }
            for (int i = 0; i < RawFitness.Length; i++)
            {
                if ((int)Rank[i] > 0)
                {
                    int prod = 0;
                    for (u = (int)Rank[i]; u > 0; u--)
                    {
                        prod += counts[u - 1] * Strength[counts[u - 1] - 1];
                    }
                    RawFitness[i] = prod;
                }
                else
                {
                    RawFitness[i] = 0;
                }
            }
        }

        /// <summary>
        /// D(i)    = 1 / (σk(i)+2)
        /// <para/>Density, k-th nearest neighbour method (Silvermann 1986)
        /// </summary>
        /// <param name="k"></param>
        private void S_MO_KNNDensity(int k)
        {
            KNNDensity = new double[xPopulationTest.Count];

            //euclidean distance to all neighbours, in objective space
            //      therefore, normalize objective space with current min and max... really? what if i have penalty constraints, 
            //          where a violated solution gets a super bad value? well its not on the nondom front anyway. 
            //          so also has a high rawfitness. and density values are super small 0-1


            //normalization
            double[] max = new double[mObj];
            double[] min = new double[mObj];
            for (int m = 0; m < mObj; m++)
            {
                min[m] = double.MaxValue;
                max[m] = double.MinValue;
            }
            for (int m = 0; m < mObj; m++)
            {
                for (int i = 0; i < objvalsTestPop_MO.Count; i++)
                {
                    if (objvalsTestPop_MO[i][m] < min[m]) min[m] = objvalsTestPop_MO[i][m];
                    else if (objvalsTestPop_MO[i][m] > max[m]) max[m] = objvalsTestPop_MO[i][m];
                }
            }

            double[,] normcost = new double[xPopulationTest.Count, mObj];
            for (int m = 0; m < mObj; m++)
            {
                for (int i = 0; i < xPopulationTest.Count; i++)
                {
                    normcost[i, m] = (objvalsTestPop_MO[i][m] - min[m]) / (max[m] - min[m]);
                    //normcost[i, m] = objvalsTestPop_MO[i][m];
                }
            }

            // sigmas
            double[][] sigma = new double[xPopulationTest.Count][];
            int[][] indices = new int[xPopulationTest.Count][];
            double[] sigmak = new double[xPopulationTest.Count];        //only taking k-th distance
            double[] sigmakSUM = new double[xPopulationTest.Count];     //summing all k distances

            for (int i = 0; i < xPopulationTest.Count; i++)
            {
                sigma[i] = new double[xPopulationTest.Count];
                indices[i] = new int[xPopulationTest.Count];
                for (int j = 0; j < xPopulationTest.Count; j++)
                {
                    indices[i][j] = j;
                    if (i != j)
                    {
                        for (int m = 0; m < mObj; m++)
                        {
                            sigma[i][j] += Math.Pow(normcost[i, m] - normcost[j, m], 2);
                        }
                        sigma[i][j] = Math.Sqrt(sigma[i][j]);
                    }
                    else sigma[i][j] = 0;
                }
                Array.Sort(sigma[i], indices[i]);

                sigmak[i] = sigma[i][k];
                for (int j = 1; j < k + 1; j++)
                {
                    sigmakSUM[i] += sigma[i][j];
                }

                //KNNDensity[i] = 1.0 / (sigmak[i] + 2);
                KNNDensity[i] = 1.0 / (sigmakSUM[i] + 1);     //chris version
            }

        }

        /// <summary>
        /// Archive truncation to maintain best solutions in the set.
        /// </summary>
        /// <param name="nArchive">Archive size to maintain.</param>
        /// <param name="mergeTestBest">Merging test population with current best population?</param>
        private void B_MO_ArchiveTruncation(int nArchive, bool mergeTestBest)
        {

            //first sort according to F(i) ...-> objvalTestPopulation
            //B_ArchiveOfBest(xTestPopulation.Count, false);

            //merge with TestPopulaton
            if (mergeTestBest)
            {
                xPopulationTest.AddRange(xPopulationArchive);
                objvalPopulationTest.AddRange(objvalPopulationBest);
                objvalsTestPop_MO.AddRange(objvalsBestPop_MO);
            }



            //Sort x and f_MO according to objvals
            double[] fAll = objvalPopulationTest.ToArray();
            double[][] items = new double[xPopulationTest.Count][];
            for (int i = 0; i < xPopulationTest.Count; i++)
            {
                items[i] = new double[xPopulationTest[i].Length + objvalsTestPop_MO[i].Length];
                xPopulationTest[i].CopyTo(items[i], 0);
                objvalsTestPop_MO[i].CopyTo(items[i], xPopulationTest[i].Length);
            }
            Array.Sort(fAll, items);


            //new archive with best n solutions
            xPopulationArchive = new List<double[]>();
            objvalPopulationBest = new List<double>();
            objvalsBestPop_MO = new List<double[]>();
            for (int n = 0; n < xPopulationTest.Count; n++)
            {
                double[] xparse = new double[xPopulationTest[n].Length];
                Array.Copy(items[n], 0, xparse, 0, xparse.Length);
                xPopulationArchive.Add(xparse);

                double[] fmoparse = new double[mObj];
                Array.Copy(items[n], xparse.Length, fmoparse, 0, mObj);
                objvalsBestPop_MO.Add(fmoparse);

                objvalPopulationBest.Add(fAll[n]);
            }



            //  Rank is implicitly in the right order. even after sorting. F(i) is correlated to Rank
            //now truncate
            int count = 0;
            for (int i = 0; i < xPopulationTest.Count; i++)
            {
                if ((int)Rank[i] == 0) count++;
            }


            //make another MO operator, which truncates.
            //      first: if |rank(0)|>nArchive
            //                  delete all non-rank(0)
            //                  truncate the rank(0)s with this density stuff
            //             else
            //                  keep only nArchive in the archive
            if (count <= nArchive)
            {
                xPopulationArchive.RemoveRange(nArchive, xPopulationArchive.Count - nArchive);
                objvalPopulationBest.RemoveRange(nArchive, objvalPopulationBest.Count - nArchive);
                objvalsBestPop_MO.RemoveRange(nArchive, objvalsBestPop_MO.Count - nArchive);
            }
            else
            {
                //deleting all non-0 ranks
                xPopulationArchive.RemoveRange(count, xPopulationArchive.Count - count);
                objvalPopulationBest.RemoveRange(count, objvalPopulationBest.Count - count);
                objvalsBestPop_MO.RemoveRange(count, objvalsBestPop_MO.Count - count);

                //truncate 0-ranks with spea2 archive truncation
                //while archive too full
                //      remove rank-0 item with lowest sigmak
                //      recalculate sigmaks for all remaining items
                //      repeat

                while (xPopulationArchive.Count > nArchive)
                {
                    //normalization
                    double[] max = new double[mObj];
                    double[] min = new double[mObj];
                    for (int m = 0; m < mObj; m++)
                    {
                        min[m] = double.MaxValue;
                        max[m] = double.MinValue;
                    }
                    for (int m = 0; m < mObj; m++)
                    {
                        for (int i = 0; i < objvalsBestPop_MO.Count; i++)
                        {
                            if (objvalsBestPop_MO[i][m] < min[m]) min[m] = objvalsBestPop_MO[i][m];
                            if (objvalsBestPop_MO[i][m] > max[m]) max[m] = objvalsBestPop_MO[i][m];
                        }
                    }

                    double[,] normcost = new double[xPopulationArchive.Count, mObj];
                    for (int m = 0; m < mObj; m++)
                    {
                        for (int i = 0; i < xPopulationArchive.Count; i++)
                        {
                            normcost[i, m] = (objvalsBestPop_MO[i][m] - min[m]) / (max[m] - min[m]);
                            //normcost[i, m] = objvalsBestPop_MO[i][m];
                        }
                    }


                    double[][] sigma = new double[xPopulationArchive.Count][];
                    double[] smallest = new double[xPopulationArchive.Count];
                    //int[][] indices = new int[xBestPopulation.Count][];
                    int[] indicess = new int[xPopulationArchive.Count];
                    double[] sigmaSUM = new double[xPopulationArchive.Count];     //summing all distances per individual

                    for (int i = 0; i < xPopulationArchive.Count; i++)
                    {
                        indicess[i] = i;
                        //indices[i] = new int[xBestPopulation.Count];
                        sigma[i] = new double[xPopulationArchive.Count];

                        for (int j = 0; j < xPopulationArchive.Count; j++)
                        {
                            //indices[i][j] = j;
                            if (i != j)
                            {
                                for (int m = 0; m < mObj; m++)
                                {
                                    sigma[i][j] += Math.Pow(normcost[i, m] - normcost[j, m], 2);
                                }
                                sigma[i][j] = Math.Sqrt(sigma[i][j]);
                            }
                            else sigma[i][j] = double.MaxValue;
                        }
                        //Array.Sort(sigma[i], indices[i]);
                        Array.Sort(sigma[i]);
                        smallest[i] = sigma[i][0];
                        for (int j = 0; j < xPopulationArchive.Count; j++)
                        {
                            sigmaSUM[i] += sigma[i][j];
                        }
                    }

                    //Array.Sort(sigmaSUM, indicess);
                    Array.Sort(smallest, indicess);
                    //sort according to sigmaSUM and delete that one with the smallest. but only one.

                    //deleting all non-0 ranks
                    xPopulationArchive.RemoveAt(indicess[0]);
                    objvalPopulationBest.RemoveAt(indicess[0]);
                    objvalsBestPop_MO.RemoveAt(indicess[0]);

                }

            }
        }
        #endregion
    }

















}
