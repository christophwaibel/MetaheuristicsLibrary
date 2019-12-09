using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Simple Evolution Strategy
    /// Beyer and Schwefel (2002). Evolution strategies - A comprehensive introduction.
    /// </summary>
    public class EvolutionStrategy : SingleObjective
    {

        //(mu/roh ,+ lambda) - ES

        //Begin
        //g:=0;
        //initialize(Pop0:={(x0_i,s0_i,fx0_i), i=1,...,mu}
        //Repeat
        // For l:=1 to lambda Do Begin
        //     Fl := marriage(Pg, roh);
        //     sl := s_recombination(Fl);
        //     yl := y_recombination(Fl);
        //     stildel := s_mutation(sl);
        //     ytildel := y_mutation(yl, stildel);
        //     Ftildel := F(ytildel);
        //  End;
        //  P0g := {(ytildel, stildel, Ftildel), l=1,...,lambda};
        //  Case selection_type Of
        //     (mu , lambda) : Pg+1 := selection(P0g, mu);
        //     (mu + lambda) : Pg+1 := selection(P0g, Pg, mu);
        //  End;
        //  g := g+1;
        //Until termination_condition
        //End





        /// <summary>
        /// Population size
        /// </summary>
        public int popsize { get; private set; }


        /// <summary>
        /// cost values of new population
        /// </summary>
        public double[] fx_pop { get; private set; }
        /// <summary>
        /// decision variables of new population
        /// </summary>
        public double[][] x_pop { get; private set; }


        private double[][] x0;
        private double[] fx0;

        /// <summary>
        /// strategy parameters. step-size.
        /// </summary>
        private double[][] s;
        /// <summary>
        /// step-size s0 only used at initial sampling
        /// </summary>
        private double[] s0;

        /// <summary>
        /// learning rate
        /// </summary>
        private double tau0, tau, tauc;

        /// <summary>
        /// mu = parents, lambda = offspring, roh = mixing number
        /// </summary>
        private int lambda, roh;

        /// <summary>
        /// initiale population, sampling mode. 0=uniform sampling; 1=gaussian sampling around a base point.
        /// </summary>
        private int x0samplingmode;

        /// <summary>
        /// mutation probability per discrete variable x.
        /// </summary>
        private double pmut_int;

        /// <summary>
        /// selection mode for recombination. 0 = random, 1 = roulette wheel,... to do 2 = stochastic universal sampling
        /// </summary>
        private int selmode;

        public EvolutionStrategy(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[][] x0 = null)
            : base(lb, ub, xint, evalmax, evalfnc, seed)
        {
            this.x0 = x0 ?? new double[0][];


            //popsize. same as mu
            if (settings.ContainsKey("popsize"))
                this.popsize = Convert.ToInt16(settings["popsize"]);
            else
                this.popsize = 20;


            //offspring
            if (settings.ContainsKey("lambda"))
                this.lambda = Convert.ToInt16(settings["lambda"]);
            else
                this.lambda = this.popsize;


            //mixing nr., i.e. how many parents involved in creating one offspring. roh=1 is only mutation.
            if (settings.ContainsKey("roh"))
                this.roh = Convert.ToInt16(settings["roh"]);
            else
                this.roh = 2;
            if (this.roh > this.popsize) this.roh = this.popsize;


            //selection mode for recombination. 0 = random, 1 = roulette, 2 = SUS 
            if (settings.ContainsKey("selmode"))
                this.selmode = Convert.ToInt16(settings["selmode"]);
            else
                this.selmode = 1;


            //stepsize s
            this.s = new double[this.popsize][];
            for (int p = 0; p < this.popsize; p++)
            {
                this.s[p] = new double[base.n];
                if (settings.ContainsKey("stepsize"))
                    for (int i = 0; i < base.n; i++)
                        this.s[p][i] = Convert.ToDouble(settings["stepsize"]);
                else
                    for (int i = 0; i < base.n; i++)
                        this.s[p][i] = 0.5;
            }


            //initial stepsize s0
            this.s0 = new double[base.n];
            if (settings.ContainsKey("stepsize0"))
                for (int i = 0; i < base.n; i++)
                    this.s0[i] = Convert.ToDouble(settings["stepsize0"]);
            else
                for (int i = 0; i < base.n; i++)
                    this.s0[i] = 0.5;


            //learning rate tau
            if (settings.ContainsKey("tauc"))
                this.tauc = Convert.ToDouble(settings["tauc"]);
            else
                this.tauc = 1;

            this.tau = this.tauc / (Math.Sqrt(2 * Math.Sqrt(base.n)));
            this.tau0 = this.tauc / (Math.Sqrt(2 * base.n));


            //mutation probability, only for integer
            if (settings.ContainsKey("pmut_int"))
                this.pmut_int = Convert.ToDouble(settings["pmut_int"]);
            else
                this.pmut_int = 0.5;









            if (settings.ContainsKey("x0sampling"))
            {
                this.x0samplingmode = Convert.ToInt16(settings["x0sampling"]);
                if (this.x0samplingmode > 1) this.x0samplingmode = 1;               //gaussian sampling around a point, which is uniformly sampled (unless an x0 is given)
            }
            else
            {
                this.x0samplingmode = 0;    //uniform sampling within the search domain
            }

        }


        private void initialize()
        {
            int existing_p = this.x0.Length;
            this.fx0 = new double[this.popsize];
            bool x0exists = false;
            if (this.x0.Length > 0)
            {
                x0exists = true;
                double[][] _xpop = new double[this.popsize][];
                for (int p = 0; p < this.x0.Length; p++)
                {
                    _xpop[p] = new double[base.n];
                    this.x0[p].CopyTo(_xpop[p], 0);
                    this.fx0[p] = base.evalfnc(this.x0[p]);
                    this.evalcount++;
                }
                this.x0 = new double[this.popsize][];   //i'm doing this, because x0.Length at initialisation could be < popsize
                _xpop.CopyTo(this.x0, 0);
            }
            else
            {
                this.x0 = new double[this.popsize][];
            }

            double[] xbasepoint = new double[base.n];
            if (x0exists)
            {
                this.x0[0].CopyTo(xbasepoint, 0);       //always take the first entry, no matter if more than 1 has been inputted
            }
            else
            {
                for (int i = 0; i < base.n; i++)
                {
                    xbasepoint[i] = base.rnd.NextDouble() * (base.ub[i] - base.lb[i]) + base.lb[i]; //if no x0 exists, uniform sampling
                    if (base.xint[i])
                    {
                        xbasepoint[i] = Math.Round(xbasepoint[i], 0);
                    }
                }
            }
            for (int p = existing_p; p < this.popsize; p++)
            {
                this.x0[p] = new double[base.n];
                if (this.x0samplingmode == 0)
                {
                    // uniform sampling withing the search domain
                    for (int i = 0; i < base.n; i++)
                    {
                        this.x0[p][i] = base.rnd.NextDouble() * (base.ub[i] - base.lb[i]) + base.lb[i];
                        if (base.xint[i])
                        {
                            this.x0[p][i] = Math.Round(this.x0[p][i], 0);
                        }
                    }
                }
                else
                {
                    //gaussian sampling around a point
                    for (int i = 0; i < base.n; i++)
                    {
                        this.x0[p][i] = (base.rnd.NextGaussian(0, this.s0[i]) * (base.ub[i] - base.lb[i])) + xbasepoint[i];
                        //this.x0[p][i] = (base.rnd.NextGaussian(0, this.s0[i]) * (base.ub[i] - base.lb[i])) + xbasepoint[i];
                        if (base.xint[i])
                        {
                            this.x0[p][i] = Math.Round(this.x0[p][i], 0);
                        }
                    }
                }
                base.checkBounds(ref this.x0[p], base.lb, base.ub);
                this.fx0[p] = base.evalfnc(this.x0[p]);
                this.evalcount++;
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


            this.x_pop = new double[this.popsize][];
            this.x0.CopyTo(this.x_pop, 0);
            this.fx_pop = new double[this.popsize];
            this.fx0.CopyTo(this.fx_pop, 0);

            Array.Sort(this.fx_pop, this.x_pop);
        }


        public override void solve()
        {
            this.initialize();

            //Begin
            //g:=0;     dont need generations, im using eval count
            while (base.evalcount < base.evalmax)
            {
                int[] int_family;
                double[][] s_new = new double[this.lambda][];
                double[][] x_new = new double[this.lambda][];
                double[] fx_new = new double[this.lambda];
                for (int l = 0; l < this.lambda; l++)
                {
                    this.marriage(out int_family, this.popsize, this.roh, this.selmode);
                    this.s_recombination(out s_new[l], this.s, int_family);
                    this.x_recombination(out x_new[l], this.x_pop, int_family);   //interm. recomb. for R, or coordinate-wise recomb for N
                    this.s_mutation(ref s_new[l], this.tau);
                    this.x_mutation(ref x_new[l], s_new[l], this.tau);
                    fx_new[l] = base.evalfnc(x_new[l]);
                    base.evalcount++;
                    if (this.CheckIfNaN(fx_new[l]))
                    {
                        if (fx_new[l] < base.fxopt)
                        {
                            base.fxopt = fx_new[l];
                            x_new[l].CopyTo(base.xopt, 0);
                        }
                        return;
                    }
                }
                double[][] x_sel;
                double[] fx_sel;
                double[][] s_sel;
                this.selection(out x_sel, out fx_sel, out s_sel, this.x_pop, this.fx_pop, this.s, x_new, fx_new, s_new);
                x_sel.CopyTo(this.x_pop, 0);
                fx_sel.CopyTo(this.fx_pop, 0);
                s_sel.CopyTo(this.s, 0);
                this.storeCurrentBest();
                //Console.WriteLine("eval: {0} with fx: {1}", base.evalcount, base.get_fxoptimum());
            }
        }


        /// <summary>
        /// pure random selection. No fitness-proportionate selection.
        /// </summary>
        /// <param name="_xfamily"></param>
        /// <param name="_sfamily"></param>
        /// <param name="_xpop"></param>
        /// <param name="_spop"></param>
        /// <param name="_roh"></param>
        private void marriage(out int[] _intfamily, int _popsize, int _roh, int selmode)
        {
            _intfamily = new int[this.roh];
            bool[] selected = new bool[_popsize];

            switch (selmode)
            {
                case 1: //roulette
                    List<int> lengths = new List<int>();
                    for (int p = 0; p < _popsize; p++) lengths.Add(_popsize - p); //first entry gets biggest value
                    for (int i = 0; i < this.roh; i++)
                    {
                        int ranksum = lengths.Sum();
                        int rand = base.rnd.Next(ranksum);
                        int sel = -1;
                        int partsum = 0;
                        do
                        {
                            sel++;
                            partsum += lengths[sel];
                        } while (partsum < rand && sel != _popsize - 1);

                        _intfamily[i] = popsize - lengths[sel];
                        lengths.RemoveAt(sel);
                    }
                    break;
                default:    //random, no fitness proportionate
                    List<int> indices = new List<int>();
                    for (int p = 0; p < _popsize; p++) indices.Add(p);
                    for (int i = 0; i < this.roh; i++)
                    {
                        int rndind = rnd.Next(indices.Count());
                        int sel = indices[rndind];
                        indices.RemoveAt(rndind);
                        _intfamily[i] = sel;
                    }
                    break;
            }

        }

        private void s_recombination(out double[] _s_new, double[][] _s, int[] _intfamily)
        {
            int _roh = _intfamily.Length;
            _s_new = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                for (int u = 0; u < _roh; u++) //length of roh, mixing number
                {
                    _s_new[i] += _s[_intfamily[u]][i];
                }
                _s_new[i] /= _roh;
            }
        }

        /// <summary>
        /// intermediate recombination for Real-valued. coordinate-wise random selection for discrete parameters.
        /// </summary>
        /// <param name="_x_new"></param>
        /// <param name="_intfamily"></param>
        private void x_recombination(out double[] _x_new, double[][] _xpop, int[] _intfamily)
        {
            int _roh = _intfamily.Length;
            _x_new = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                if (base.xint[i])   //integer
                {
                    _x_new[i] = _xpop[_intfamily[rnd.Next(roh)]][i];
                }
                else                //real
                {
                    for (int u = 0; u < _roh; u++) //length of roh, mixing number
                    {
                        _x_new[i] += _xpop[_intfamily[u]][i];
                    }
                    _x_new[i] /= _roh;
                }
            }
        }


        private void s_mutation(ref double[] _snew, double _tau)
        {
            double tau0exp = Math.Exp(this.tau0 * rnd.NextGaussian(0, 1));
            for (int i = 0; i < base.n; i++)
            {
                _snew[i] = tau0exp * _snew[i] * Math.Exp(this.tau * rnd.NextGaussian(0, 1));
            }
        }

        private void x_mutation(ref double[] _xnew, double[] _snew, double _tau)
        {
            for (int i = 0; i < base.n; i++)
            {
                if (base.xint[i])
                {
                    if (rnd.NextDouble() < this.pmut_int)
                    {
                        _xnew[i] = Convert.ToDouble(rnd.Next(Convert.ToInt16(base.lb[i]), Convert.ToInt16(base.ub[i]) + 1));
                    }
                }
                else
                {
                    _xnew[i] = _xnew[i] + (_snew[i] * rnd.NextGaussian(0, 1));
                }
            }
            base.checkBounds(ref _xnew, base.lb, base.ub);
        }

        private void selection(out double[][] _x_sel, out double[] _fx_sel, out double[][] _s_sel,
            double[][] _x_old, double[] _fx_old, double[][] _s_old, double[][] _x_new, double[] _fx_new, double[][] _s_new)
        {
            int popsize = _x_old.Length;
            int lambda = _x_new.Length;

            double[][] _x_merged = new double[popsize + lambda][];
            double[] _fx_merged = new double[popsize + lambda];
            double[][] _s_merged = new double[popsize + lambda][];
            _x_old.CopyTo(_x_merged, 0);
            _x_new.CopyTo(_x_merged, popsize);
            _fx_old.CopyTo(_fx_merged, 0);
            _fx_new.CopyTo(_fx_merged, popsize);
            _s_old.CopyTo(_s_merged, 0);
            _s_new.CopyTo(_s_merged, popsize);

            _x_old = new double[popsize][];
            _fx_old = new double[popsize];
            _s_old = new double[popsize][];
            double[][][] _xands = new double[popsize + lambda][][];
            for (int i = 0; i < _xands.Length; i++)
                _xands[i] = new double[2][];

            for (int i = 0; i < popsize + lambda; i++)
            {
                _xands[i][0] = new double[_x_merged[0].Length];
                _xands[i][1] = new double[_s_merged[0].Length];
                _x_merged[i].CopyTo(_xands[i][0], 0);
                _s_merged[i].CopyTo(_xands[i][1], 0);
            }

            //rank according to fitness
            Array.Sort(_fx_merged, _xands);

            _x_sel = new double[popsize][];
            _fx_sel = new double[popsize];
            _s_sel = new double[popsize][];
            for (int i = 0; i < popsize; i++)
            {
                _x_sel[i] = new double[_x_merged[0].Length];
                _xands[i][0].CopyTo(_x_sel[i], 0);

                _fx_sel[i] = _fx_merged[i];

                _s_sel[i] = new double[_s_merged[0].Length];
                _xands[i][1].CopyTo(_s_sel[i], 0);
            }
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
                this.storeCurrentBest();

                stop = true;
            }
            return stop;
        }

    }
}
