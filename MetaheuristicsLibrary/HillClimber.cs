using System;


namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Stochastic Hill-Climbing,
    /// <para/>Pseudocode from: Clever Algorithms: Nature-Inspired Programming Recipes (Jason Brownlee)
    /// </summary>
    public class Hillclimber : SingleObjective
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


        private double[] xtest;
        private double[] x;
        private double fxtest;
        private double fx;
        private double[] x0;


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
        /// <param name="evalmax">Maximum iterations.</param>
        /// <param name="evalfnc">Evaluation function.</param>
        /// <param name="seed">Seed for random number generator.</param>
        public Hillclimber(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed, double stepsize, double[] x0 = null) :
            base(lb, ub, xint, evalmax, evalfnc, seed)
        {
            this.stepsize = stepsize;


            this.x0 = x0 ?? new double[0];
        }

        /// <summary>
        /// Minimizes an evaluation function using stochastic hill climbing.
        /// </summary>
        public override void solve()
        {
            int n = lb.Length;
            this.x = new double[n];

            double[] stdev = new double[n];

            if (this.x0.Length == base.n)
            {
                this.x0.CopyTo(this.x, 0);
            }
            else
            {
                for (int i = 0; i < n; i++)
                {
                    this.x[i] = rnd.NextDouble() * (ub[i] - lb[i]) + lb[i];
                    stdev[i] = stepsize * (ub[i] - lb[i]);
                }
            }
            this.fx = evalfnc(this.x);

            for (base.evalcount = 0; base.evalcount < evalmax; base.evalcount++)
            {
                this.xtest = new double[n];
                for (int i = 0; i < n; i++)
                {
                    this.xtest[i] = rnd.NextGaussian(this.x[i], stdev[i]);
                    if (this.xtest[i] > ub[i]) this.xtest[i] = ub[i];
                    else if (this.xtest[i] < lb[i]) this.xtest[i] = lb[i];
                }
                this.fxtest = evalfnc(this.xtest);

                if (CheckIfNaN(this.fxtest)) return;


                storeCurrentBest();
            }



        }


        protected override void storeCurrentBest()
        {
            if (this.fxtest < this.fx)
            {
                this.xtest.CopyTo(this.x, 0);
                this.fx = this.fxtest;

                base.xopt = new double[n];
                this.x.CopyTo(base.xopt, 0);
                base.fxopt = this.fx;
            }
        }

        protected override bool CheckIfNaN(double fxtest)
        {
            if (Double.IsNaN(fxtest))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
    }
}
