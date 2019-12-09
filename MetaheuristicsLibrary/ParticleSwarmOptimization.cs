using System;
using System.Collections.Generic;


namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Particle Swarm Optimizer (PSO)
    /// Fully informed PSO, which uses multiple neighbouring particle vectors for updating velocity.
    /// The canonical PSO only considers 2 vectors for updating: a particle's own best and the global / local best particle.
    /// Source:     Mendes R., Kennedy J., Neves J. (2004). Fully Informed Particle Swarm: Simpler, Maybe Better.
    /// 2nd Source: Poli R., Kennedy J., Blackwell T. (2007). Particle swarm optimization - An overview.
    /// </summary>
    public class ParticleSwarmOptimization : SingleObjective
    {
        //Poli wt al. (2007), eqt. 5
        //
        //v_i <- χ (v_i + 1/K (Σ_n∈K (U(0,ϕ) ⨂ (x_n,i,best - x_i) ) ) )
        //x_i <- x_i + v_i
        //
        // with: 
        // x_n,i,best   = best foud so far of neighbour n of particle i
        // K            = number of neighbours
        // χ            = constriction coefficient, typically 0.7298
        // ϕ            = 4.1?




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
        /// population size
        /// </summary>
        private int popsize;
        /// <summary>
        /// constriction coefficient. typically 0.7298
        /// </summary>
        private double chi;
        /// <summary>
        /// upper bound for random number for attraction of particle i to best particle j. 4.1?
        /// </summary>
        private double phi;
        /// <summary>
        /// constriction coefficient for random initial velocity. fraction of the search domain.
        /// </summary>
        private double v0max;
        /// <summary>
        /// velocities
        /// </summary>
        private double[][] v;

        /// <summary>
        /// best solution vector found so far, per particle
        /// </summary>
        private double[][] px_best;
        /// <summary>
        /// best cost value found so far, per particle
        /// </summary>
        private double[] pfx_best;
        /// <summary>
        /// Per particle, 4 indices of neighbours. excluding own best. cp. Mendes et al. (2004).
        /// </summary>
        private int[][] indK;

        /// <summary>
        /// 0 = uniform, 1= gaussian
        /// </summary>
        private int x0samplingmode;

        /// <summary>
        /// updating pxbest mode. mode 0 = update after looping trough entire population. mode 1 = update after each function evaluation.
        /// </summary>
        private int pxupdatemode;


        /// <summary>
        /// initial stepsize, in case of gaussian initial sampling
        /// </summary>
        private double[] s0;

        /// <summary>
        /// phi1 = own best, phi2 = global best
        /// </summary>
        private double phi1, phi2;


        /// <summary>
        /// false means inertia weight PSO, true is fipso
        /// </summary>
        private int psomode;

        /// <summary>
        /// probability for uniform random mutation of integer variable
        /// </summary>
        private double intmut;

        private double mutstdev;



        /// <summary>
        /// using von Neumann topology.
        /// </summary>
        /// <param name="lb"></param>
        /// <param name="ub"></param>
        /// <param name="xint"></param>
        /// <param name="evalmax"></param>
        /// <param name="evalfnc"></param>
        /// <param name="seed"></param>
        /// <param name="settings"></param>
        /// <param name="x0"></param>
        public ParticleSwarmOptimization(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[][] x0 = null)
            : base(lb, ub, xint, evalmax, evalfnc, seed)
        {
            this.x0 = x0 ?? new double[0][];

            if (settings.ContainsKey("psomode"))
            {
                this.psomode = Convert.ToInt16(settings["psomode"]);
            }
            else
            {
                this.psomode = 0;       //fipso... 1=pso inertia, 2 = pso constriction coeff
            }

            if (settings.ContainsKey("phi1"))
            {
                this.phi1 = Convert.ToDouble(settings["phi1"]);
            }
            else
            {
                this.phi1 = 2.05;
            }

            if (settings.ContainsKey("phi2"))
            {
                this.phi2 = Convert.ToDouble(settings["phi2"]);
            }
            else
            {
                this.phi2 = 2.05;
            }


            if (settings.ContainsKey("popsize"))
            {
                this.popsize = Convert.ToInt32(settings["popsize"]);
            }
            else
            {
                this.popsize = 24;      // population size.
            }

            if (settings.ContainsKey("chi"))
            {
                this.chi = Convert.ToDouble(settings["chi"]);
            }
            else
            {
                this.chi = 0.1;      // constriction coefficient
            }


            if (settings.ContainsKey("phi"))
            {
                this.phi = Convert.ToDouble(settings["phi"]);
            }
            else
            {
                this.phi = 4;         // attraction to best particle 
            }

            if (settings.ContainsKey("v0max"))
            {
                this.v0max = Convert.ToDouble(settings["v0max"]);
            }
            else
            {
                this.v0max = 0.2;       // max velocity at initialisation. fraction of domain.
            }

            if (settings.ContainsKey("x0samplingmode"))
            {
                this.x0samplingmode = Convert.ToInt16(settings["x0samplingmode"]);
            }
            else
            {
                this.x0samplingmode = 0;// 0 = uniform, 1 = gaussian
            }

            if (settings.ContainsKey("pxupdatemode"))
            {
                this.pxupdatemode = Convert.ToInt16(settings["pxupdatemode"]);
            }
            else
            {
                this.pxupdatemode = 0;  //0 = update after population. 1 = update after each evaluation
            }

            this.s0 = new double[base.n];

            if (settings.ContainsKey("s0"))
            {
                for (int i = 0; i < base.n; i++) s0[i] = Convert.ToDouble(settings["s0"]);
            }
            else
            {
                for (int i = 0; i < base.n; i++) s0[i] = 1.0;       //initial step size in case of gaussian sampling
            }



            if (settings.ContainsKey("intmut"))
            {
                this.intmut = Convert.ToDouble(settings["intmut"]);
            }
            else
            {
                this.intmut = 0.5;
            }

            if (settings.ContainsKey("mutstdev"))
            {
                this.mutstdev = Convert.ToDouble(settings["mutstdev"]);
            }
            else
            {
                this.mutstdev = 0.3;
            }
        }

        public override void solve()
        {
            if (psomode == 0)
            {
                this.initializeNeighbourhood();
                this.initialSamples();

                int pUpdate = this.popsize;
                if (this.pxupdatemode == 1) pUpdate = 1;

                int pNow = 0;

                while (base.evalcount < base.evalmax)
                {
                    for (int p = 0; p < pUpdate; p++)
                    {
                        double[] v_new;
                        double[] x_new;
                        double fx_new;
                        this.updateParticleFIPSO(out v_new, out x_new, this.v[pNow], this.x_pop[pNow], this.px_best, this.indK[pNow]);
                        fx_new = base.evalfnc(x_new);
                        this.v[pNow] = new double[base.n];
                        this.x_pop[pNow] = new double[base.n];
                        v_new.CopyTo(this.v[pNow], 0);
                        x_new.CopyTo(this.x_pop[pNow], 0);
                        this.fx_pop[pNow] = fx_new;

                        this.updateParticlesBest(ref this.px_best[pNow], ref this.pfx_best[pNow], x_new, fx_new);

                        //Console.WriteLine("fx,x1,x2,{0},{1},{2}", fx_new, x_new[0], x_new[1]);
                        base.evalcount++;
                        pNow++;
                        if (pNow == this.popsize) pNow = 0;

                        this.storeCurrentBest();
                        if (this.CheckIfNaN(fx_new))
                        {
                            return;
                        }
                    }
                }
            }
            else
            {
                this.initialSamples();
                int pUpdate = this.popsize;
                if (this.pxupdatemode == 1) pUpdate = 1;

                int pNow = 0;

                while (base.evalcount < base.evalmax)
                {
                    for (int p = 0; p < pUpdate; p++)
                    {
                        double[] v_new;
                        double[] x_new;
                        double fx_new;
                        this.updateParticle(out v_new, out x_new, this.v[pNow], this.x_pop[pNow], this.xopt, this.px_best[pNow], this.chi, this.phi1, this.phi2, this.psomode);
                        fx_new = base.evalfnc(x_new);
                        this.v[pNow] = new double[base.n];
                        v_new.CopyTo(this.v[pNow], 0);
                        this.x_pop[pNow] = new double[base.n];
                        x_new.CopyTo(this.x_pop[pNow], 0);      //magically, its also copied to pxbest???!!!
                        this.fx_pop[pNow] = fx_new;

                        this.updateParticlesBest(ref this.px_best[pNow], ref this.pfx_best[pNow], x_new, fx_new);

                        base.evalcount++;
                        pNow++;
                        if (pNow == this.popsize) pNow = 0;

                        this.storeCurrentBest();
                        if (this.CheckIfNaN(fx_new))
                        {
                            return;
                        }
                    }
                }
            }
        }


        /// <summary>
        /// von Neumann topology, excluding own particle
        /// </summary>
        private void initializeNeighbourhood()
        {
            double sqrtPop = Math.Sqrt(this.popsize);
            int cols = (int)Math.Floor(sqrtPop);
            int rows = (int)Math.Round((Convert.ToDouble(this.popsize) / Convert.ToDouble(cols)), 0);
            this.popsize = cols * rows;         //making sure popsize can be divided into a cols*rows matrix

            this.indK = new int[this.popsize][];
            for (int i = 0; i < this.popsize; i++)
                this.indK[i] = new int[4];

            int count = 0;
            for (int i = 0; i < cols; i++)
            {
                for (int j = 0; j < rows; j++)
                {
                    this.indK[count][0] = count - cols;
                    this.indK[count][1] = count - 1;
                    this.indK[count][2] = count + 1;
                    this.indK[count][3] = count + cols;

                    if ((count + 1) % cols == 0) this.indK[count][2] = count - cols + 1;
                    else if ((count + 1) % cols == 1) this.indK[count][1] = count + cols - 1;

                    if (count - cols < 0) this.indK[count][0] = count + this.popsize - cols;
                    if (count + cols > this.popsize - 1) this.indK[count][3] = count - this.popsize + cols;

                    count++;
                }
            }
        }

        /// <summary>
        /// sampling initial population x0 and fx0, initial velocity v0, and store as current best particles.
        /// </summary>
        private void initialSamples()
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
                    for (int i = 0; i < base.n; i++)
                    {
                        if (base.xint[i])
                        {
                            this.x0[p][i] = Math.Round(this.x_pop[0][i], 0);
                        }
                    }
                    _xpop[p] = new double[base.n];
                    this.x0[p].CopyTo(_xpop[p], 0);
                    this.fx0[p] = base.evalfnc(this.x0[p]);
                    //Console.WriteLine("fx,x1,x2,{0},{1},{2}", this.fx0[p], this.x0[p][0], this.x0[p][1]);
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
                    // uniform sampling within the search domain
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
                this.checkBounds(ref this.x0[p], base.lb, base.ub);
                this.fx0[p] = base.evalfnc(this.x0[p]);
                //Console.WriteLine("fx,x1,x2,{0},{1},{2}", this.fx0[p], this.x0[p][0], this.x0[p][1]);
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
            this.fx_pop = new double[this.popsize];
            this.px_best = new double[this.popsize][];
            this.pfx_best = new double[this.popsize];
            this.x0.CopyTo(this.x_pop, 0);
            this.x0.CopyTo(this.px_best, 0);
            this.fx0.CopyTo(this.fx_pop, 0);
            this.fx0.CopyTo(this.pfx_best, 0);



            this.v = new double[this.popsize][];
            for (int p = 0; p < this.popsize; p++)
            {
                this.v[p] = new double[base.n];
                for (int i = 0; i < base.n; i++)
                {
                    //!!!!!!!!!!!!!!!! BUG
                    // v0 should be able to go into both directions, negative and positive!
                    // here, its only going in one direction basically, and then shifted to lb. if thats -, it can happen
                    // that we have - velocity.
                    //but it should be:
                    this.v[p][i] = base.rnd.NextGaussian(0, base.ub[i] - base.lb[i]) * this.v0max;

                    //this.v[p][i] = base.rnd.NextDouble() * (base.ub[i] - base.lb[i]) + base.lb[i];
                    //this.v[p][i] *= this.v0max;
                }
            }
        }

        private void updateParticleFIPSO(out double[] v_new, out double[] x_new, double[] v_old, double[] x_old, double[][] px_best, int[] neighbours)
        {
            v_new = new double[base.n];
            x_new = new double[base.n];
            //v_i <- χ (v_i + 1/K (Σ_n∈K (U(0,ϕ) ⨂ (x_n,i,best - x_i) ) ) )
            //x_i <- x_i + v_i
            int K = neighbours.Length;
            for (int i = 0; i < base.n; i++)
            {
                double sumXi = 0;
                for (int k = 0; k < K; k++)
                {
                    sumXi += (rnd.NextDouble() * this.phi) * (px_best[neighbours[k]][i] - x_old[i]);
                }
                v_new[i] = v_old[i] + (sumXi / K);
                v_new[i] *= this.chi;
                if (Math.Round(v_new[i], 5) == 0)
                {
                    double rndgauss = base.rnd.NextGaussian(0, this.chi * (base.ub[i] - base.lb[i]));
                    v_new[i] = rndgauss;
                }

                x_new[i] = x_old[i] + v_new[i];
                if (base.xint[i])
                {
                    x_new[i] = Math.Round(x_new[i], 0);
                    if (base.rnd.NextDouble() < this.intmut)
                    {
                        double rndgauss = base.rnd.NextGaussian(0, this.mutstdev * (base.ub[i] - base.lb[i]));
                        //double rndgauss = base.rnd.NextGaussian(0, v_new[i] * (base.ub[i] - base.lb[i]));
                        x_new[i] = Convert.ToInt32(x_new[i] + rndgauss);
                    }
                }
            }

            base.checkBounds(ref x_new, base.lb, base.ub);
        }

        private void updateParticle(out double[] v_new, out double[] x_new, double[] v_old, double[] x_old,
            double[] gx_best, double[] ownx_best,
            double inertia, double phi1, double phi2,
            int psomode)
        {
            v_new = new double[base.n];
            x_new = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                if (psomode == 1)
                {
                    v_new[i] = inertia * v_old[i] +
                        (rnd.NextDouble() * phi1) * (ownx_best[i] - x_old[i]) +
                        (rnd.NextDouble() * phi2) * (gx_best[i] - x_old[i]);
                }
                else
                {
                    v_new[i] = inertia * (v_old[i] +
                        (rnd.NextDouble() * phi1) * (ownx_best[i] - x_old[i]) +
                        (rnd.NextDouble() * phi2) * (gx_best[i] - x_old[i]));
                }
                if (Math.Round(v_new[i], 5) == 0)
                {
                    double rndgauss = base.rnd.NextGaussian(0, inertia * (base.ub[i] - base.lb[i]));
                    v_new[i] = rndgauss;
                }
                x_new[i] = x_old[i] + v_new[i];
                if (base.xint[i])
                {
                    x_new[i] = Math.Round(x_new[i], 0);
                    if (base.rnd.NextDouble() < this.intmut)
                    {
                        double rndgauss = base.rnd.NextGaussian(0, this.mutstdev * (base.ub[i] - base.lb[i]));
                        //double rndgauss = base.rnd.NextGaussian(0, v_new[i] * (base.ub[i] - base.lb[i]));
                        x_new[i] = Convert.ToInt32(x_new[i] + rndgauss);
                    }
                }
            }

            base.checkBounds(ref x_new, base.lb, base.ub);
        }


        private void updateParticlesBest(ref double[] px_best, ref double pfx_best, double[] x_now, double fx_now)
        {
            if (fx_now < pfx_best)
            {
                x_now.CopyTo(px_best, 0);
                pfx_best = fx_now;
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
                stop = true;
            }
            return stop;
        }


    }

}
