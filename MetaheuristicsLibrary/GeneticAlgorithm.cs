using System;
using System.Collections.Generic;


namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Simple Genetic Algorithm
    /// Goldberg (1989). Genetic Algorithms in search, optimization and machine learning.
    /// </summary>
    public class GeneticAlgorithm : SingleObjective
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
        /// cost values of new population
        /// </summary>
        public double[] fx_pop { get; private set; }
        /// <summary>
        /// decision variables of new population
        /// </summary>
        public double[][] x_pop { get; private set; }

        /// <summary>
        /// initial solutions
        /// </summary>
        private double[][] x0;
        private double[] fx0;

        /// <summary>
        /// how many elites to be maintained per generation
        /// </summary>
        private int elite;

        /// <summary>
        /// Simple GA, according to Goldberg 1989, chapter 3.
        /// </summary>
        /// <param name="lb"></param>
        /// <param name="ub"></param>
        /// <param name="evalmax"></param>
        /// <param name="evalfnc"></param>
        /// <param name="seed"></param>
        /// <param name="settings">Dictionary. Should contain: population size ("popsize", int), crossover probability ("pcross", double), mutation probability ("pmut", double).</param>
        /// <param name="x0">Decision variables of initial population</param>
        /// <param name="fx0">Cost values of initial population</param>
        public GeneticAlgorithm(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[][] x0 = null)
            : base(lb, ub, xint, evalmax, evalfnc, seed)
        {
            this.x0 = x0 ?? new double[0][];

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
            }
            else
            {
                popsize = 20;
            }

            if (settings.ContainsKey("elite"))
            {
                elite = Convert.ToInt32(settings["elite"]);
            }
            else
            {
                elite = 1;
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
        /// Minimizes an evaluation function using a simple Genetic Algorithm.
        /// </summary>
        public override void solve()
        {
            //initialize
            initialize();

            this.x_pop = new double[this.popsize][];
            this.fx_pop = new double[this.popsize];
            this.x0.CopyTo(this.x_pop, 0);
            this.fx0.CopyTo(this.fx_pop, 0);
            //main loop
            do
            {
                Array.Sort(this.fx_pop, this.x_pop);
                double[] fitness_pop = CostToFitness(this.fx_pop);
                this.sumfitness = 0;
                for (int p = 0; p < this.popsize; p++)
                {
                    this.sumfitness += fitness_pop[p];
                }


                int nOffspring = 0;
                double[][] x_new = new double[this.popsize][];
                double[] fx_new = new double[this.popsize];
                for (int e = 0; e < this.elite; e++)
                {
                    x_new[nOffspring] = new double[base.n];
                    this.x_pop[e].CopyTo(x_new[nOffspring], 0);
                    fx_new[nOffspring] = this.fx_pop[e];
                    nOffspring++;
                }
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
                    fx_new[nOffspring] = evalfnc(child1);
                    //Console.WriteLine("fx,x1,x2,{0},{1},{2}", fx_new[nOffspring], child1[0], child1[1]);
                    base.evalcount++;
                    if (CheckIfNaN(fx_new[nOffspring]))
                    {
                        return;
                    }
                    if (fx_new[nOffspring] < base.fxopt)
                    {
                        base.fxopt = fx_new[nOffspring];
                        x_new[nOffspring].CopyTo(base.xopt, 0);
                    }


                    if (nOffspring + 1 < this.popsize)
                    {
                        x_new[nOffspring + 1] = child2;
                        fx_new[nOffspring + 1] = evalfnc(child2);
                        //Console.WriteLine("fx,x1,x2,{0},{1},{2}", fx_new[nOffspring + 1], child2[0], child2[1]);
                        base.evalcount++;
                        if (CheckIfNaN(fx_new[nOffspring + 1]))
                        {
                            return;
                        }
                        if (fx_new[nOffspring + 1] < base.fxopt)
                        {
                            base.fxopt = fx_new[nOffspring + 1];
                            x_new[nOffspring + 1].CopyTo(base.xopt, 0);
                        }
                        nOffspring++;
                    }
                    nOffspring++;
                } while (nOffspring < this.popsize);

                x_new.CopyTo(this.x_pop, 0);
                fx_new.CopyTo(this.fx_pop, 0);


                gen++;
            } while (gen < maxgen && base.evalcount < evalmax);
        }


        protected override void storeCurrentBest()
        {
            for (int p = 0; p < this.popsize; p++)
            {
                if (this.fx_pop[p] < base.fxopt)
                {
                    base.fxopt = this.fx_pop[p];
                    this.x_pop[p].CopyTo(base.xopt, 0);
                }
            }
        }

        protected override bool CheckIfNaN(double fxtest)
        {
            bool stop = false;
            if (Double.IsNaN(fxtest))
            {
                //get the best
                storeCurrentBest();

                stop = true;
            }
            return stop;
        }


        /// <summary>
        /// initialize algorithm, i.e. initial population. Uniform random sampling over domain.
        /// </summary>
        private void initialize()
        {
            this.sumfitness = 0;

            int existing_p = this.x0.Length;
            this.fx0 = new double[this.popsize];
            if (this.x0.Length > 0)
            {
                double[][] _xpop = new double[this.popsize][];
                for (int p = 0; p < this.x0.Length; p++)
                {
                    _xpop[p] = new double[base.n];
                    this.x0.CopyTo(_xpop[p], 0);
                    this.fx0[p] = base.evalfnc(this.x0[p]);
                    //Console.WriteLine("fx,x1,x2,{0},{1},{2}", this.fx0[p], this.x0[p][0], this.x0[p][1]);
                    base.evalcount++;
                }
                this.x0 = new double[this.popsize][];   //i'm doing this, because x0.Length at initialisation could be < popsize
                _xpop.CopyTo(this.x0, 0);
            }
            else
            {
                this.x0 = new double[this.popsize][];
            }

            for (int p = existing_p; p < this.popsize; p++)
            {
                this.x0[p] = new double[base.n];
                for (int i = 0; i < base.n; i++)
                {
                    this.x0[p][i] = base.rnd.NextDouble() * (base.ub[i] - base.lb[i]) + base.lb[i];
                    if (base.xint[i])
                    {
                        this.x0[p][i] = Math.Round(this.x0[p][i], 0);
                    }
                }
                this.fx0[p] = base.evalfnc(this.x0[p]);
                //Console.WriteLine("fx,x1,x2,{0},{1},{2}", this.fx0[p], this.x0[p][0], this.x0[p][1]);
                base.evalcount++;
            }


            double[] fitness_pop = CostToFitness(this.fx0);
            for (int p = 0; p < this.popsize; p++)
            {
                this.sumfitness += fitness_pop[p];
            }

            //get the best
            for (int p = 0; p < this.popsize; p++)
            {
                if (this.fx0[p] < base.fxopt)
                {
                    base.fxopt = this.fx0[p];
                    this.x0[p].CopyTo(base.xopt, 0);
                }
            }
        }

        private double[] CostToFitness(double[] fx_pop)
        {
            int[] items = new int[this.popsize];
            double[] fitness = new double[this.popsize];
            double[] keys = new double[this.popsize];
            fx_pop.CopyTo(keys, 0);

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

            base.checkBounds(ref child1, base.lb, base.ub);
            base.checkBounds(ref child2, base.lb, base.ub);

        }


    }
}
