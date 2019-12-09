using MetaheuristicsLibrary.Misc;
using System;
using System.Collections.Generic;


namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Rosenbrock Search Algorithm.
    /// <para/>Rosenbrock H.H. (1960). An Automatic Method for finding the Greatest or Least Value of a Function.
    /// <para/>Palmer J.R. (1969). An improved procedure for orthogonalising the search vectors in Rosenbrock's and Swann's direct search optimisation methods.
    /// <para/>Based on source code from: Kniaz.net 
    /// <para/>Description: Local optimizer. Conducts orthogonal steps in a coordinate system that is rotated in order to approximate the gradient.
    /// </summary>
    public class Rosenbrock : SingleObjective
    {

        private double[] x0;
        private double fx0;
        private double[] xtest;
        private double fxtest;


        //Rosenbrock variables
        private bool[] flag_success;    //for every directon (dimension), indicating if there was a successfull move in the iteration
        private bool[] flag_fail;       //same, but for fail moves
        private int flags_sum;          // summing up TRUEs. so I dont always need to loop through all flags to check. I know: sum=2n, then all true
        private int flag_currentDirection; //indicating, in which direction to move in this move-iteration
        private double[] xOrigin;       //starting point of iteration (i.e. until all flags fullfilled). then xOrigin changes xBest - xOrigin
        private double[][] vMoves;      //n orthogonal vectors


        /// <summary>
        /// Strategy parameters for extension, shrinking, and initial steplength
        /// </summary>
        private double alpha, beta, stepsize;

        public Rosenbrock(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[] x0 = null)
            : base(lb, ub, xint, evalmax, evalfnc, seed)
        {
            this.x0 = x0;//?? new double[base.n];


            //extension
            if (settings.ContainsKey("alpha"))
                this.alpha = Convert.ToDouble(settings["alpha"]);
            else
                this.alpha = 2;

            //shrinking
            if (settings.ContainsKey("beta"))
                this.beta = Convert.ToDouble(settings["beta"]);
            else
                this.beta = 0.5;

            //initial step size
            if (settings.ContainsKey("stepsize"))
                this.stepsize = Convert.ToDouble(settings["stepsize"]);
            else
                this.stepsize = 0.125;


        }


        private void initialize()
        {

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
                    this.x0[i] = (base.rnd.NextGaussian(0, this.stepsize) * (base.ub[i] - base.lb[i])) + center[i];
                    if (base.xint[i])
                    {
                        this.x0[i] = Math.Round(this.x0[i], 0);
                    }
                }
            }
            base.checkBounds(ref this.x0, base.lb, base.ub);
            this.xtest = new double[base.n];
            this.x0.CopyTo(this.xtest, 0);
            this.fx0 = base.evalfnc(this.xtest);
            this.fxtest = this.fx0;
            base.evalcount++;

            this.xtest = this.NormalizeX(this.xtest);


            Console.WriteLine("x0, x1, fx, {0},{1},{2}", this.ReverseNormalization(this.xtest)[0], this.ReverseNormalization(this.xtest)[1], this.fxtest);
            this.storeCurrentBest();






            //Initialize Rosenbrock strategy parameters and variables
            List<bool[]> rb_flag_fail = new List<bool[]>();
            List<bool[]> rb_flag_success = new List<bool[]>();
            List<int> rb_flags_sum = new List<int>();
            List<int> rb_flag_currentDirection = new List<int>();
            List<double[][]> rb_vMoves = new List<double[][]>();
            List<double[]> rb_xOrigin = new List<double[]>();

            I_Rosenbrock(ref rb_flag_fail, ref rb_flag_success, ref rb_flags_sum, ref rb_flag_currentDirection,
               ref rb_vMoves, ref rb_xOrigin, this.stepsize);

            this.flag_fail = rb_flag_fail[0];
            this.flag_success = rb_flag_success[0];
            this.flags_sum = rb_flags_sum[0];
            this.flag_currentDirection = rb_flag_currentDirection[0];
            this.vMoves = rb_vMoves[0];

            this.xOrigin = new double[base.n];
            rb_xOrigin[0].CopyTo(this.xOrigin, 0);
        }


        public override void solve()
        {
            this.initialize();

            while (base.evalcount < base.evalmax)
            {
                this.M_Rosenbrock(ref this.flag_fail, ref this.flag_success, ref this.flags_sum, ref this.flag_currentDirection, ref this.xOrigin, ref this.vMoves, this.alpha, this.beta, this.stepsize);
                if (this.CheckIfNaN(this.fxtest)) return;
                this.storeCurrentBest();
            }
        }

        protected override void storeCurrentBest()
        {
            if (this.fxtest < base.fxopt)
            {
                base.xopt = new double[n];
                this.ReverseNormalization(this.xtest).CopyTo(base.xopt, 0);
                base.fxopt = this.fxtest;
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

        /// <summary>
        /// Rosenbrock moves.
        /// </summary>
        /// <param name="flag_fail"></param>
        /// <param name="flag_success"></param>
        /// <param name="flags_sum"></param>
        /// <param name="flag_currentDirection"></param>
        /// <param name="xOrigin"></param>
        /// <param name="vMoves"></param>
        /// <param name="alpha"></param>
        /// <param name="beta"></param>
        private void M_Rosenbrock(ref bool[] flag_fail, ref bool[] flag_success,
            ref int flags_sum, ref int flag_currentDirection,
            ref double[] xOrigin, ref double[][] vMoves,
            double alpha, double beta, double steplength)
        {
            // **********************************************************************************
            // 0.
            // **********************************************************************************
            //check if 4 conditions fullfilled to make new orthogonaliztion
            if (flags_sum == 2 * base.n)
            {
                double[] basevector = new double[base.n];
                for (int i = 0; i < base.n; i++)
                {
                    basevector[i] = xOrigin[i] - this.NormalizeX(base.xopt)[i];
                }
                Array.Copy(this.NormalizeX(base.xopt), xOrigin, base.n);   //new origin point
                double[][] vMoves_copy = Transformation.GramSchmidt(basevector, base.n, vMoves);     //new move vectors
                for (int i = 0; i < base.n; i++)
                {
                    for (int u = 0; u < base.n; u++)
                    {
                        vMoves_copy[i][u] *= steplength;
                    }
                }
                vMoves_copy.CopyTo(vMoves, 0);

                flag_fail = new bool[base.n];
                flag_success = new bool[base.n];
                flags_sum = 0;
                flag_currentDirection = 0;        //dunno, should this be also reseted?
            }
            //if yes, make a vector between xbase(previous starting point) and xbest(current best value)
            //use this as base vector to make gram-schmidt. basically, is the eigenvector




            // **********************************************************************************
            // 1.
            // **********************************************************************************
            // n-searches sequentially in direction of base vectors
            //      so 1st search only change xTest = x0 + [Alpha, 0, 0, ..., 0n]
            //      then 2nd search           xTest' = xTest + [0, Alpha, 0, ..., 0n]
            //      and that for every n
            if (flag_currentDirection == base.n) flag_currentDirection = 0;

            for (int i = 0; i < base.n; i++)
            {
                this.xtest[i] = this.NormalizeX(base.xopt)[i] + vMoves[flag_currentDirection][i];
                if (this.xtest[i] < 0) this.xtest[i] = 0;
                else if (this.xtest[i] > 1) this.xtest[i] = 1;
            }
            this.fxtest = base.evalfnc(this.ReverseNormalization(this.xtest));
            base.evalcount++;
            Console.WriteLine("x0, x1, fx, {0},{1},{2}", this.ReverseNormalization(this.xtest)[0], this.ReverseNormalization(this.xtest)[1], this.fxtest);
            if (this.CheckIfNaN(this.fxtest)) return;


            if (this.fxtest < base.fxopt)
            {
                //its better, so take it and increase stepsize with alpha
                for (int i = 0; i < base.n; i++)
                    vMoves[flag_currentDirection][i] *= alpha;

                //if flag is not 1 already, make it 1 and add to flags_sum
                if (flag_success[flag_currentDirection] != true)
                {
                    flag_success[flag_currentDirection] = true;
                    flags_sum++;
                }

                //make it new best point
                Array.Copy(this.ReverseNormalization(this.xtest), base.xopt, base.n);
                base.fxopt = this.fxtest;
            }
            else
            {
                //its worse, so reverse vMove for this direction and contract steplength with beta
                for (int i = 0; i < base.n; i++)
                    vMoves[flag_currentDirection][i] *= -1 * beta;

                //if fail flag not 1, make it 1 an add to flags sum
                if (flag_fail[flag_currentDirection] != true)
                {
                    flag_fail[flag_currentDirection] = true;
                    flags_sum++;
                }
            }



            flag_currentDirection++;
        }

        /// <summary>
        /// initializes and returns stuff for Rosenbrock moves
        /// </summary>
        /// <param name="flag_fail"></param>
        /// <param name="flag_success"></param>
        /// <param name="flags_sum"></param>
        /// <param name="flag_currentDirection"></param>
        /// <param name="vMoves"></param>
        /// <param name="xOrigin"></param>
        /// <param name="steplength"></param>
        private void I_Rosenbrock(ref List<bool[]> rb_flag_fail, ref List<bool[]> rb_flag_success,
            ref List<int> rb_flags_sum, ref List<int> rb_flag_currentDirection,
            ref List<double[][]> rb_vMoves, ref List<double[]> rb_xOrigin, double steplength)
        {
            rb_flag_fail = new List<bool[]>();
            rb_flag_success = new List<bool[]>();
            rb_flags_sum = new List<int>();
            rb_flag_currentDirection = new List<int>();
            rb_vMoves = new List<double[][]>();
            rb_xOrigin = new List<double[]>();

            //conditions for new orthogonalization
            bool[] _flag_fail = new bool[base.n];
            bool[] _flag_success = new bool[base.n];
            int _flags_sum = 0;
            int _flag_currentDirection = 0;
            double[][] _vMoves = new double[base.n][];
            double[] _xOrigin = new double[base.n];


            //origin point for current iteration( base)
            rb_xOrigin.Add(this.xtest);

            //initial orthogonal move vectors
            //rb_vMoves[i] = new double[this.base.n][];
            for (int u = 0; u < base.n; u++)
            {
                _vMoves[u] = new double[base.n];
                for (int k = 0; k < base.n; k++)
                {
                    if (u == k) _vMoves[u][k] = steplength;
                    else _vMoves[u][k] = 0;
                }
            }




            rb_flag_fail.Add(_flag_fail);
            rb_flag_success.Add(_flag_success);
            rb_flags_sum.Add(_flags_sum);
            rb_flag_currentDirection.Add(_flag_currentDirection);
            rb_vMoves.Add(_vMoves);

        }


        private double[] NormalizeX(double[] _x)
        {
            double[] _xNorm = new double[_x.Length];
            //Array.Copy(x0, xTest, base.n);
            for (int i = 0; i < base.n; i++)
            {
                _xNorm[i] = (_x[i] - base.lb[i]) / Math.Abs(base.ub[i] - base.lb[i]);
            }
            return _xNorm;
        }

        private double[] ReverseNormalization(double[] _xNorm)
        {
            double[] _x = new double[_xNorm.Length];
            for (int i = 0; i < _xNorm.Length; i++)
            {
                _x[i] = Math.Abs(base.ub[i] - base.lb[i]) * _xNorm[i] + base.lb[i];
            }
            return _x;
        }
    }
}
