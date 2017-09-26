using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MetaheuristicsLibrary.Misc;

namespace MetaheuristicsLibrary.SolversSO
{
    public abstract class SO_Solver
    {
        /// <summary>
        /// Decision variable count
        /// </summary>
        public int n { get; private set; }
        /// <summary>
        /// Lower bound for each variable.
        /// </summary>
        public double[] lb { get; private set; }
        /// <summary>
        /// Upper bound for each variable.
        /// </summary>
        public double[] ub { get; private set; }
        /// <summary>
        /// Indicate, which variable is integer
        /// </summary>
        public bool[] xint { get; private set; }
        /// <summary>
        /// Maximum function evaluations.
        /// </summary>
        public int itermax { get; private set; }
        /// <summary>
        /// Evaluation function.
        /// </summary>
        public Func<double[], double> evalfnc { get; private set; }
        /// <summary>
        /// Variable vector of final solution.
        /// </summary>
        public double[] xopt { get; protected set; }
        /// <summary>
        /// Cost of final solution.
        /// </summary>
        public double fxopt { get; protected set; }

        protected RandomDistributions rnd;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="lb"></param>
        /// <param name="ub"></param>
        /// <param name="xinteger">indicate, which variable is integer</param>
        /// <param name="itermax"></param>
        /// <param name="evalfnc"></param>
        /// <param name="seed"></param>
        public SO_Solver(double[] lb, double[] ub, bool[] xint, int itermax, Func<double[], double> evalfnc, int seed)
        {
            this.lb = lb;
            this.ub = ub;
            this.xint = xint;
            this.itermax = itermax;
            this.evalfnc = evalfnc;
            this.n = lb.Length;
            this.rnd = new RandomDistributions(seed);
        }



        /// <summary>
        /// Solve. Override this.
        /// </summary>
        public abstract void solve();



        /// <summary>
        /// Get the variable vector of the final solution. Call this after solving complete.
        /// </summary>
        /// <returns>Variable vector.</returns>
        public double[] get_Xoptimum()
        {
            return this.xopt;
        }

        /// <summary>
        /// Get the cost value of the final solution. Call this after solving is complete.
        /// </summary>
        /// <returns>Cost value.</returns>
        public double get_fxoptimum()
        {
            return this.fxopt;
        }




    }

    public class Hillclimber : SO_Solver
    {
        //source: Stochastic Hill-Climbing, in: Clever Algorithms: Nature-Inspired Programming Recipes (Jason Brownlee)
        //
        //Input: Itermax, ProblemSize 
        //Output: Current 
        //Current  <- RandomSolution(ProblemSize)
        //For (iteri ∈ Itermax )
        //    Candidate  <- RandomNeighbor(Current)
        //    If (Cost(Candidate) >= Cost(Current))
        //        Current  <- Candidate
        //    End
        //End
        //Return (Current)

        /// <summary>
        /// Stepsize.
        /// </summary>
        public double stepsize { get; private set; }


        /// <summary>
        /// Initialize a stochastic hill climber optimization algorithm. Assuming minimization problems.
        /// </summary>
        /// <param name="lb">Lower bound for each variable.</param>
        /// <param name="ub">Upper bound for each variable.</param>
        /// <param name="stepsize">Stepsize.</param>
        /// <param name="itermax">Maximum iterations.</param>
        /// <param name="evalfnc">Evaluation function.</param>
        /// <param name="seed">Seed for random number generator.</param>
        public Hillclimber(double[] lb, double[] ub, bool[] xint, int itermax, Func<double[], double> evalfnc, int seed, double stepsize) :
            base(lb, ub, xint, itermax, evalfnc, seed)
        {
            this.stepsize = stepsize;
        }

        /// <summary>
        /// Minimizes an evaluation function using stochastic hill climbing.
        /// </summary>
        public override void solve()
        {
            int n = lb.Length;
            double[] x = new double[n];
            double[] stdev = new double[n];

            for (int i = 0; i < n; i++)
            {
                x[i] = rnd.NextDouble() * (ub[i] - lb[i]) + lb[i];
                stdev[i] = stepsize * (ub[i] - lb[i]);
            }
            double fx = evalfnc(x);

            for (int t = 0; t < itermax; t++)
            {
                double[] xtest = new double[n];
                for (int i = 0; i < n; i++)
                {
                    xtest[i] = rnd.NextGaussian(x[i], stdev[i]);
                    if (xtest[i] > ub[i]) xtest[i] = ub[i];
                    else if (xtest[i] < lb[i]) xtest[i] = lb[i];
                }
                double fxtest = evalfnc(xtest);

                if (Double.IsNaN(fxtest)) return;

                if (fxtest < fx)
                {
                    xtest.CopyTo(x, 0);
                    fx = fxtest;

                    xopt = new double[n];
                    x.CopyTo(xopt, 0);
                    fxopt = fx;
                }
            }



        }

    }



    public class SimpleGA : SO_Solver
    {
        /// <summary>
        /// Population size
        /// </summary>
        public int popsize { get; private set; }
        /// <summary>
        /// Current generation
        /// </summary>
        public int gen { get; private set; }
        /// <summary>
        /// Maximum generations
        /// </summary>
        public int maxgen { get; private set; }
        /// <summary>
        /// current count of function evaluation calls
        /// </summary>
        public int currentiter { get; private set; }
        /// <summary>
        /// length of bitstring
        /// </summary>
        private int[] lchrom;

        /// <summary>
        /// crossover probability
        /// </summary>
        public double pcross { get; private set; }
        /// <summary>
        /// mutation probability
        /// </summary>
        public double pmutation { get; private set; }

        /// <summary>
        /// extension of hypercube, for intermediate recombination
        /// </summary>
        public double d { get; private set; }

        /// <summary>
        /// 
        /// </summary>
        public double r { get; private set; }

        /// <summary>
        /// 
        /// </summary>
        public double k { get; private set; }

        /// <summary>
        /// sum of population fitness. required for roulette wheel
        /// </summary>
        private double sumfitness;


        /// <summary>
        /// statistics
        /// </summary>
        private double stat_avg, stat_max, stat_min;


        /// <summary>
        /// cost values of previous population
        /// </summary>
        private double[] fx_pop_old;
        /// <summary>
        /// decision variables of previous population
        /// </summary>
        private double[][] x_pop_old;
        /// <summary>
        /// cost values of new population
        /// </summary>
        public double[] fx_pop { get; private set; }
        /// <summary>
        /// decision variables of new population
        /// </summary>
        public double[][] x_pop { get; private set; }

        /// <summary>
        /// cost values of initial solutions
        /// </summary>
        private double[] fx_0;
        /// <summary>
        /// initial solutions
        /// </summary>
        private double[][] x_0;



        /// <summary>
        /// Simple GA, according to Goldberg 1989, chapter 3.
        /// </summary>
        /// <param name="lb"></param>
        /// <param name="ub"></param>
        /// <param name="itermax"></param>
        /// <param name="evalfnc"></param>
        /// <param name="seed"></param>
        /// <param name="settings">Dictionary. Should contain: population size ("popsize", int), crossover probability ("pcross", double), mutation probability ("pmut", double).</param>
        /// <param name="x0">Decision variables of initial population</param>
        /// <param name="fx0">Cost values of initial population</param>
        public SimpleGA(double[] lb, double[] ub, bool[] xint, int itermax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[][] x0 = null, double[] fx0 = null)
            : base(lb, ub, xint, itermax, evalfnc, seed)
        {
            this.x_0 = x0 ?? new double[0][];
            this.fx_0 = fx0 ?? new double[0];


            gen = 0;

            this.lchrom = new int[base.n];
            for (int i = 0; i < base.n; i++)
            {
                if (base.xint[i])
                {
                    bool[] xint_ = new bool[1];
                    xint_[0] = base.xint[i];
                    double[] lb_ = new double[1];
                    double[] ub_ = new double[1];
                    lb_[0] = base.lb[i];
                    ub_[0] = base.ub[i];

                    this.lchrom[i] = Misc.Misc.BinReqLength(xint_, lb_, ub_);
                }
                else
                {
                    this.lchrom[i] = 0;
                }
            }





            if (settings.ContainsKey("maxgen"))
            {
                maxgen = Convert.ToInt32(settings["maxgen"]);
            }
            else
            {
                maxgen = 100;
            }

            if (settings.ContainsKey("popsize"))
            {
                popsize = Convert.ToInt32(settings["popsize"]);
                if (popsize > 100) popsize = 100;
            }
            else
            {
                popsize = 20;
            }

            if (settings.ContainsKey("pcross"))
            {
                pcross = Convert.ToDouble(settings["pcross"]);
            }
            else
            {
                pcross = 0.7;
            }

            if (settings.ContainsKey("pmut"))
                pmutation = Convert.ToDouble(settings["pmut"]);
            else
                pmutation = 0.3;


            if (settings.ContainsKey("d"))
                d = Convert.ToDouble(settings["d"]);
            else
                d = 0.1;

            if (settings.ContainsKey("r"))
                r = Convert.ToDouble(settings["r"]);
            else
                r = 0.1;

            if (settings.ContainsKey("k"))
                k = Convert.ToDouble(settings["k"]);
            else
                k = 16;
        }

        /// <summary>
        /// Minimizes an evaluation function using stochastic hill climbing.
        /// </summary>
        public override void solve()
        {
            //initialize
            initialize();

            this.x_pop = new double[this.popsize][];
            this.fx_pop = new double[this.popsize];
            this.x_0.CopyTo(this.x_pop, 0);
            this.fx_0.CopyTo(this.fx_pop, 0);
            //main loop
            do
            {
                double[] fitness_pop = CostToFitness(this.fx_pop);
                this.sumfitness = 0;
                for (int p = 0; p < this.popsize; p++)
                {
                    this.sumfitness += fitness_pop[p];
                }


                int nOffspring = 0;
                double[][] x_new = new double[this.popsize][];
                double[] fx_new = new double[this.popsize];
                do
                {
                    double[] child1;
                    double[] child2;
                    double[] parent1 = new double[base.n];
                    double[] parent2 = new double[base.n];
                    parent1 = this.x_pop[select(this.sumfitness, fitness_pop)];
                    parent2 = this.x_pop[select(this.sumfitness, fitness_pop)];
                    crossovermutate(out child1, out child2, parent1, parent2, d, r, k);



                    x_new[nOffspring] = child1;
                    x_new[nOffspring + 1] = child2;
                    fx_new[nOffspring] = evalfnc(child1);
                    this.currentiter++;
                    if (CheckIfNaN(fx_new[nOffspring]))
                    {
                        return;
                    }
                    fx_new[nOffspring + 1] = evalfnc(child2);
                    this.currentiter++;
                    if (CheckIfNaN(fx_new[nOffspring + 1]))
                    {
                        return;
                    }


                    nOffspring += 2;
                } while (nOffspring < this.popsize);

                x_new.CopyTo(this.x_pop, 0);
                fx_new.CopyTo(this.fx_pop, 0);



                //get the best
                for (int p = 0; p < this.popsize; p++)
                {
                    if (this.fx_pop[p] < base.fxopt)
                    {
                        base.fxopt = this.fx_pop[p];
                        this.x_pop[p].CopyTo(base.xopt, 0);
                    }
                }


                gen++;
            } while (gen < maxgen && currentiter < itermax);










        }

        private bool CheckIfNaN(double fxtest)
        {
            bool stop = false;
            if (Double.IsNaN(fxtest))
            {
                //get the best
                for (int p = 0; p < this.popsize; p++)
                {
                    if (this.fx_pop[p] < base.fxopt)
                    {
                        base.fxopt = this.fx_pop[p];
                        this.x_pop[p].CopyTo(base.xopt, 0);
                    }
                }

                stop = true;
            }
            return stop;
        }


        /// <summary>
        /// initialize algorithm, i.e. initial population. Uniform random sampling over domain.
        /// </summary>
        private void initialize()
        {
            base.fxopt = double.MaxValue;
            base.xopt = new double[base.n];
            this.sumfitness = 0;
            if (this.x_0.Length > 0)
            {
                for (int p = 0; p < this.popsize; p++)
                {
                    this.sumfitness += fx_0[p];
                }
                return;
            }
            this.x_0 = new double[this.popsize][];
            this.fx_0 = new double[this.popsize];

            for (int p = 0; p < this.popsize; p++)
            {
                this.x_0[p] = new double[base.n];
                for (int i = 0; i < base.n; i++)
                {
                    this.x_0[p][i] = base.rnd.NextDouble() * (base.ub[i] - base.lb[i]) + base.lb[i];
                    if (base.xint[i])
                    {
                        this.x_0[p][i] = Math.Round(this.x_0[p][i], 0);
                    }
                }
                this.fx_0[p] = base.evalfnc(this.x_0[p]);
                this.currentiter++;

            }

            double[] fitness_pop = CostToFitness(this.fx_0);
            for (int p = 0; p < this.popsize; p++)
            {
                this.sumfitness += fitness_pop[p];
            }

            //get the best
            for (int p = 0; p < this.popsize; p++)
            {
                if (this.fx_0[p] < base.fxopt)
                {
                    base.fxopt = this.fx_0[p];
                    this.x_0[p].CopyTo(base.xopt, 0);
                }
            }
        }

        private double[] CostToFitness(double[] fx_pop)
        {
            int[] items = new int[this.popsize];
            double[] fitness = new double[this.popsize];
            double[] keys = new double[this.popsize];
            fx_pop.CopyTo(keys,0);
            //double fxsum = 0;
            //for (int p = 0; p < this.popsize; p++)
            //{
            //    fxsum += fx_pop[p];
            //}
            //for (int p = 0; p < this.popsize; p++)
            //{
            //    fitness[p] = 1 - (fx_pop[p] / fxsum);
            //    //this.sumfitness += fitness[p];
            //}

            for (int p = 0; p < this.popsize; p++)
            {
                items[p] = p;
            }
            Array.Sort(keys, items);
            for (int p = 0; p < this.popsize; p++)
            {
                fitness[items[p]] = Convert.ToDouble(this.popsize - p); 
            }

            return fitness;
        }


        /// <summary>
        /// Roulette wheel selection
        /// </summary>
        /// <param name="popsize"></param>
        /// <param name="sumfitness"></param>
        /// <param name="fitness_pop"></param>
        /// <returns></returns>
        private int select(double sumfitness, double[] fitness_pop)
        {
            int j = -1;
            double partsum = 0.0;
            double rand = base.rnd.NextDouble() * sumfitness;

            do
            {
                j++;
                partsum += fitness_pop[j];
            } while (partsum < rand && j != this.popsize - 1);

            return j;
        }



        private void crossovermutate(out double[] child1, out double[] child2, double[] parent1, double[] parent2, double d, double r, double k)
        {
            int jcross;  //crossover point

            child1 = new double[base.n];
            child2 = new double[base.n];

            for (int i = 0; i < base.n; i++)
            {
                if (base.rnd.NextDouble() <= this.pcross)
                {
                    if (base.xint[i])   //integer variable
                    {
                        if (this.lchrom[i] > 1)
                        {
                            jcross = base.rnd.Next(0, lchrom[i]);
                            string strP1 = Misc.Misc.Dec2Bin(Convert.ToInt32(parent1[i] - base.lb[i]), this.lchrom[i]);
                            string strP2 = Misc.Misc.Dec2Bin(Convert.ToInt32(parent2[i] - base.lb[i]), this.lchrom[i]);
                            char[] charP1 = strP1.ToCharArray();
                            char[] charP2 = strP2.ToCharArray();

                            char[] charC1 = new char[charP1.Length];
                            char[] charC2 = new char[charP1.Length];

                            for (int j = 0; j < jcross; j++)
                            {
                                charC1[j] = charP1[j];
                                charC2[j] = charP2[j];
                            }
                            for (int j = jcross; j < this.lchrom[i]; j++)
                            {
                                charC1[j] = charP2[j];
                                charC2[j] = charP1[j];
                            }
                            for (int j = 0; j < this.lchrom[i]; j++)
                            {
                                if (base.rnd.NextDouble() <= this.pmutation)
                                {
                                    if (charC1[j].Equals('0')) charC1[j] = '1';
                                    else charC1[j] = '0';
                                }
                                if (base.rnd.NextDouble() <= this.pmutation)
                                {
                                    if (charC2[j].Equals('0')) charC2[j] = '1';
                                    else charC2[j] = '0';
                                }
                            }

                            string strC1 = new string(charC1);
                            string strC2 = new string(charC2);
                            child1[i] = Misc.Misc.Bin2Dec(strC1) + base.lb[i];
                            if (child1[i] > base.ub[i]) child1[i] = base.ub[i]; //< lb not happening
                            child2[i] = Misc.Misc.Bin2Dec(strC2) + base.lb[i];
                            if (child2[i] > base.ub[i]) child2[i] = base.ub[i];
                        }
                        else    //swap 0 and 1.
                        {
                            child1[i] = parent2[i];
                            child2[i] = parent1[i];
                        }

                    }
                    else                //real variable... intermediate recombination
                    {
                        double alpha = this.rnd.NextDouble() * ((1 + d) + d) - d;
                        child1[i] = parent1[i] * alpha + parent2[i] * (1 - alpha);
                        child2[i] = parent2[i] * alpha + parent1[i] * (1 - alpha);


                        if (base.rnd.NextDouble() <= this.pmutation)
                        {
                             child1[i] = child1[i] + (base.rnd.NextDouble() * 2 - 1) * (r * (Math.Abs(base.ub[i] - base.lb[i]))) * (Math.Pow(2, (base.rnd.NextDouble() * -1) * k));
                             child2[i] = child2[i] + (base.rnd.NextDouble() * 2 - 1) * (r * (Math.Abs(base.ub[i] - base.lb[i]))) * (Math.Pow(2, (base.rnd.NextDouble() * -1) * k));
                        }
                    }
                }
                else
                {
                    child1[i] = parent1[i];
                    child2[i] = parent2[i];
                }

            }

            checkBounds(ref child1, base.lb, base.ub);
            checkBounds(ref child2, base.lb, base.ub);

        }

        private void checkBounds(ref double[] x, double[] _lb, double[] _ub)
        {
            for (int i = 0; i < base.n; i++)
            {
                if (x[i] < _lb[i]) x[i] = _lb[i];
                else if (x[i] > _ub[i]) x[i] = _ub[i];
            }
        }





    }
}
