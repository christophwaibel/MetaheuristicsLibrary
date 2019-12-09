using System;
using System.Collections.Generic;
using System.Linq;


namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Nelder-Mead Simplex.
    ///  <para/>Nelder J.A. and Mead R. (1965). A Simplex method for function minimization.
    ///  <para/>Lagarias J.C., Reeds J.A., Wright M.H. and Wright P.E. (1998). Convergence Properties of the Nelder-Mead Simplex Method in Low Dimensions.
    ///  <para/>Based on source code from: MATLAB 2014 fmnisearch.m
    ///  <para/>Description: Local hill-climber. Derivative-free direct search method. Constructs a n+1 simplex, whose vertices are reflected, expanded and contracted. 
    /// </summary>
    public class NelderMead : SingleObjective
    {
        /// <summary>
        /// Population size. In Nelder Mead it is n+1
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


        private double[] x0;
        private double fx0;

        /// <summary>
        /// Strategy parameters for reflection, expansion, contraction, and shrinkage coefficients.
        /// </summary>
        private double alpha, gamma, rho, sigma;
        /// <summary>
        /// Stepsize for initial sampling, as standard deviation of a normal distribution around a point
        /// </summary>
        private double step0;

        public NelderMead(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[] x0 = null)
            : base(lb, ub, xint, evalmax, evalfnc, seed)
        {
            this.x0 = x0;//?? new double[base.n];

            //popsize always n+1
            this.popsize = base.n + 1;

            //reflection
            if (settings.ContainsKey("alpha"))
                this.alpha = Convert.ToDouble(settings["alpha"]);
            else
                this.alpha = 1;

            //expansion
            if (settings.ContainsKey("gamma"))
                this.gamma = Convert.ToDouble(settings["gamma"]);
            else
                this.gamma = 1.5;

            //contraction
            if (settings.ContainsKey("rho"))
                this.rho = Convert.ToDouble(settings["rho"]);
            else
                this.rho = 0.25;

            //shrinkage
            if (settings.ContainsKey("sigma"))
                this.sigma = Convert.ToDouble(settings["sigma"]);
            else
                this.sigma = 0.2;

            //stepsize for initial sampling
            if (settings.ContainsKey("step0"))
                this.step0 = Convert.ToDouble(settings["step0"]);
            else
                this.step0 = 0.02;
        }

        private void initialize()
        {
            this.x_pop = new double[this.popsize][];
            this.fx_pop = new double[this.popsize];

            double[] center = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                center[i] = (base.lb[i] + base.ub[i]) / 2;
            }

            if (this.x0 == null || this.x0.Length == 0)
            {
                // sample new x0
                //gaussian sampling around a point
                this.x0 = new double[base.n];
                for (int i = 0; i < base.n; i++)
                {
                    this.x0[i] = (base.rnd.NextGaussian(0, this.step0) * (base.ub[i] - base.lb[i])) + center[i];
                    if (base.xint[i])
                    {
                        this.x0[i] = Math.Round(this.x0[i], 0);
                    }
                }
            }
            base.checkBounds(ref this.x0, base.lb, base.ub);
            this.x_pop[0] = new double[base.n];
            this.x0.CopyTo(this.x_pop[0], 0);
            this.fx_pop[0] = base.evalfnc(this.x_pop[0]);
            base.evalcount++;

            for (int p = 1; p < this.popsize; p++)
            {
                this.x_pop[p] = new double[base.n];
                for (int i = 0; i < base.n; i++)
                {
                    this.x_pop[p][i] = (base.rnd.NextGaussian(0, this.step0) * (base.ub[i] - base.lb[i])) + this.x0[i];
                    if (base.xint[i])
                    {
                        this.x_pop[p][i] = Math.Round(this.x_pop[p][i], 0);
                    }
                }
                base.checkBounds(ref this.x_pop[p], base.lb, base.ub);
                this.fx_pop[p] = base.evalfnc(this.x_pop[p]);
                base.evalcount++;
            }

            Array.Sort(this.fx_pop, x_pop);
            this.storeCurrentBest();
        }


        public override void solve()
        {
            this.initialize();


            while (base.evalcount < base.evalmax)
            {
                double[][] Simplex_cop = this.x_pop.ToArray();
                double[][] Simplex = new double[Simplex_cop.Length][];
                for (int i = 0; i < Simplex.Length; i++)
                {
                    Simplex[i] = new double[Simplex_cop[i].Length];
                    Simplex_cop[i].CopyTo(Simplex[i], 0);
                }

                double[] SimplexObjval = this.fx_pop.ToArray();


                // create xc - centroid of all vertices, except n+1
                double[] xc = this.NM_centroid(Simplex);

                // reflect, expand, contract, shrink
                double[] xr = this.NM_Reflect(Simplex[base.n], xc, this.alpha);

                // evaluate xr
                for (int i = 0; i < base.n; i++)
                {
                    if (xr[i] < 0) xr[i] = 0;
                    else if (xr[i] > 1) xr[i] = 1;
                }

                double xrObjval = base.evalfnc(xr);
                base.evalcount++;

                // check, what to do and when to end this loop
                // if its better then the 2nd worst but not better then the best,
                if (xrObjval < SimplexObjval[base.n - 1] && xrObjval > SimplexObjval[0])
                {
                    //replace the worst and start from beginning ->end sub
                    this.NM_replaceVertex(xr, xrObjval, ref Simplex, ref SimplexObjval);
                }
                // if its the new best
                else if (xrObjval < SimplexObjval[0])
                {
                    // expand the reflected point further
                    double[] xs = this.NM_expand(xr, xc, gamma);
                    this.checkBounds(ref xs, base.lb, base.ub);
                    double xsObjval = base.evalfnc(xs);
                    base.evalcount++;

                    //if expanded point is better, then reflected
                    if (xsObjval < xrObjval)
                    {
                        // take it as new one to replace the current worst
                        this.NM_replaceVertex(xs, xsObjval, ref Simplex, ref SimplexObjval);
                    }
                    else
                    {
                        // keep the reflected point 
                        this.NM_replaceVertex(xr, xrObjval, ref Simplex, ref SimplexObjval);
                    }
                }
                // if reflected is worse then 2nd worst
                else if (xrObjval > SimplexObjval[base.n - 1])
                {
                    double[] xcon = new double[base.n];
                    double xconObjval = 0;
                    bool shrink = false;
                    // if reflected is better then current worst
                    if (xrObjval < SimplexObjval[base.n])
                    {
                        //contract outside
                        this.NM_contract(ref xcon, ref xconObjval, ref shrink, xr, xrObjval, xc, rho);
                    }
                    // if reflected is worse then current worst
                    else if (xrObjval >= SimplexObjval[base.n])
                    {
                        // contract inside
                        this.NM_contract(ref xcon, ref xconObjval, ref shrink, Simplex[base.n], SimplexObjval[base.n], xc, rho);
                    }
                    if (shrink == true)
                    {
                        this.NM_shrink(ref Simplex, ref SimplexObjval, sigma);
                    }
                    else
                    {
                        this.NM_replaceVertex(xcon, xconObjval, ref Simplex, ref SimplexObjval);
                    }
                }
                //e.g. if its worse then the worst or has the same objval as the worst
                else
                {
                    Console.WriteLine("Degeneracy! at evalfunccalls {0}", base.evalcount);
                    return;
                }
                Array.Sort(SimplexObjval, Simplex);
                this.x_pop = new double[this.popsize][];//Simplex);
                this.fx_pop = new double[this.popsize]; //(SimplexObjval);
                Simplex.CopyTo(this.x_pop, 0);
                SimplexObjval.CopyTo(this.fx_pop, 0);



                this.storeCurrentBest();
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

        private double[] NM_centroid(double[][] Simplex)
        {
            double[][] SimplexCen = new double[Simplex.Length - 1][];
            for (int i = 0; i < Simplex.Length - 1; i++)
            {
                SimplexCen[i] = new double[base.n];
                SimplexCen[i] = Simplex[i];
            }
            return Misc.Vector.Centroid(SimplexCen);
        }

        private double[] NM_Reflect(double[] xLast, double[] xc, double alpha)
        {
            //r = 2xc -x(n+1)
            //xc = sum(x(i))/n .... the centroid 
            double[] xr = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                xr[i] = xc[i] + alpha * (xc[i] - xLast[i]);
            }

            return xr;
        }

        private void NM_replaceVertex(double[] xReplace, double objValReplace, ref double[][] Simplex, ref double[] SimplexObjval)
        {
            Simplex[base.n] = xReplace;
            SimplexObjval[base.n] = objValReplace;
        }

        private double[] NM_expand(double[] xReflected, double[] xc, double gamma)
        {
            //xs = xc + 2(xc -x(n+1))
            //xc = sum(x(i))/n .... the centroid 
            double[] xs = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                xs[i] = xc[i] + gamma * (xReflected[i] - xc[i]);
            }

            return xs;
        }

        private void NM_contract(ref double[] xcon, ref double xconObjval, ref bool shrink,
           double[] xCompare, double xCompareObjval, double[] xc, double rho)
        {
            for (int i = 0; i < base.n; i++)
            {
                xcon[i] = xc[i] + rho * (xCompare[i] - xc[i]);
            }
            this.checkBounds(ref xcon, base.lb, base.ub);
            xconObjval = base.evalfnc(xcon);
            base.evalcount++;

            // outside: if the contracted point is better then the reflected, take it, otherwise shrink
            // inside: if the contracted point is better then the worst, take it, otherwise shrink
            if (xconObjval < xCompareObjval)
            {
                shrink = false;
            }
            else
            {
                shrink = true;
            }
        }

        private void NM_shrink(ref double[][] Simplex, ref double[] SimplexObjval, double sigma)
        {
            // keep the best vertex
            // shrink the rest
            for (int i = 1; i < Simplex.Length; i++)
            {
                double[] xShrinked = new double[base.n];
                for (int j = 0; j < base.n; j++)
                {
                    xShrinked[j] = Simplex[0][j] + sigma * (Simplex[i][j] - Simplex[0][j]);
                }
                this.checkBounds(ref xShrinked, base.lb, base.ub);
                double fxShrinked = base.evalfnc(xShrinked);
                base.evalcount++;

                Simplex[i] = new double[base.n];
                xShrinked.CopyTo(Simplex[i], 0);
                SimplexObjval[i] = fxShrinked;
            }
        }
    }
}
