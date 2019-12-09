using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MetaheuristicsLibrary.Misc;

namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// Interface for Single objective solver.
    /// </summary>
    public abstract class SingleObjective
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
        public int evalmax { get; private set; }
        /// <summary>
        /// current count of function evaluation calls
        /// </summary>
        public int evalcount { get; protected set; }
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
        /// <param name="evalmax"></param>
        /// <param name="evalfnc"></param>
        /// <param name="seed"></param>
        public SingleObjective(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed)
        {
            this.lb = lb;
            this.ub = ub;
            this.xint = xint;
            this.evalcount = 0;
            this.evalmax = evalmax;
            this.evalfnc = evalfnc;
            this.n = lb.Length;
            this.rnd = new RandomDistributions(seed);

            this.fxopt = double.MaxValue;
            this.xopt = new double[this.n];
        }



        /// <summary>
        /// Solve. Override this.
        /// </summary>
        public abstract void solve();

        /// <summary>
        /// Storing current best solution.
        /// Could be, that the algorithm doesn't have memory.
        /// </summary>
        protected abstract void storeCurrentBest();

        /// <summary>
        /// stopping criterion for FrOG
        /// </summary>
        /// <returns></returns>
        protected abstract bool CheckIfNaN(double fxtest);


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

        /// <summary>
        /// checking bounds
        /// </summary>
        /// <param name="x"></param>
        /// <param name="_lb"></param>
        /// <param name="_ub"></param>
        protected void checkBounds(ref double[] x, double[] _lb, double[] _ub)
        {
            for (int i = 0; i < this.n; i++)
            {
                if (x[i] < _lb[i]) x[i] = _lb[i];
                else if (x[i] > _ub[i]) x[i] = _ub[i];
            }
        }



    }




    /*
    /// <summary>
    /// 
    /// Glasmachers et al (2010). Exponential Natural Evolution Strategies. GECCO
    /// </summary>
    public class xNES : SO_Solver
    {
        //function [xopt,fopt] = xnes(f,d,x,timeout)

        //% Written by Sun Yi (yi@idsia.ch).

        //% parameters
        //L = 4+3*floor(log(d));
        //etax = 1; etaA = 0.5*min(1.0/d,0.25);
        //shape = max(0.0, log(L/2+1.0)-log(1:L)); shape = shape / sum(shape);

        //% initialize
        //xopt = x; fopt = f(x);
        //A = zeros(d);
        //weights = zeros(1,L);
        //fit = zeros(1,L);
        //tm = cputime;

        //while cputime - tm < timeout
        //    expA = expm(A);

        //    % step 1: sampling & importance mixing
        //    Z = randn(d,L); X = repmat(x,1,L)+expA*Z;
        //    for i = 1 : L, fit(i) = f(X(:,i)); end

        //    % step 2: fitness reshaping
        //    [~, idx] = sort(fit); weights(idx) = shape;
        //    if fit(idx(1)) < fopt
        //        xopt = X(:,idx(1)); fopt = fit(idx(1));
        //    end

        //    % step 3: compute the gradient for C and x
        //    G = (repmat(weights,d,1).*Z)*Z' - sum(weights)*eye(d);
        //    dx = etax * expA * (Z*weights');
        //    dA = etaA * G;

        //    % step 4: compute the update  
        //    x = x + dx; A = A + dA;

        //    if trace(A)/d < -10*log(10), break; end
        //end

        public xNES(double[] lb, double[] ub, bool[] xint, int itermax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[][] x0 = null)
            : base(lb, ub, xint, itermax, evalfnc, seed)
        {

        }


        public override void solve()
        {
            throw new NotImplementedException();
        }

        protected override void storeCurrentBest()
        {
            throw new NotImplementedException();
        }

        protected override bool CheckIfNaN(double fxtest)
        {
            throw new NotImplementedException();
        }
    }


    public class sNES : SO_Solver
    {

        //        function[mean_vec, var_vec, mean_vec_fit] = snes(fitness, num_iter, dim, pop_size, learn_rates)
        //%[mean_vec, var_vec, mean_vec_fit] = 
        //%           snes(dim, mean_vec, var_vec, pop_size, learn_rates, @fitness)
        //%INPUTS:
        //% @fitness: fitness function (e.g., @rosenfit for rosenfit.m)
        //%           which takes an individual as input and outputs the fitness
        //% num_iter: number of iterations
        //% dim: dimension
        //% **pop_size: number of units to evaluate 
        //% **learn_rates: vector of two, for mean and variance
        //% defaults if pop_size and learn_rates not supplied
        //%OUTPUTS
        //% mean_vec: final mean
        //% var_vec: final variance
        //% mean_vec_fit: the fitness of the mean, measured at each iteration
        //%
        //%code by Matt Luciw (matt.luciw at gmail)
        //%contact Tom Schaul with all your SNES questions

        //%initial mean and variance
        //mean_vec = rand(dim,1);
        //var_vec = ones(dim,1);  %ones!  <-- this is important

        //if (nargin < 4)
        //    %abra cadabra
        //    pop_size = 4 + floor(3 * log(dim));

        //    learn_rates = [1 (3 + log(dim))/(5 * sqrt(dim))];
        //end

        //mean_vec_fit = zeros(1,num_iter);

        //%outer loop: number of population evaluations
        //for i = 1:num_iter

        //    if (mod(i,500)==0) 
        //        fprintf(1, '\nGeneration %d...', i)
        //    end

        //    %draw from standard normal distribution
        //    curr_samples = randn(pop_size,dim);

        //    %add the input mean and variance
        //    curr_members = (curr_samples .* repmat(var_vec',pop_size,1)) + ...
        //        repmat(mean_vec',pop_size,1);

        //    %store samples
        //    S = curr_samples';

        //    %inner loop: number of population members
        //    for j = 1 : pop_size

        //        %fitness evaluated here for this sample (and stored)
        //        fit(j) = fitness(curr_members(j,:));

        //    end

        //    %sort by fitness so most fit guys are last
        //    [dummy order] = sort(fit);

        //    %ordered set of samples
        //    S = S(:,order);

        //    %utilities which must sum to one
        //    %first half of the population has zero utility
        //    threshold = floor(pop_size / 2);
        //    step_size = 1 / threshold; 
        //    U = zeros(1,pop_size);
        //    U(end-threshold+1:end) = step_size:step_size:1;
        //    U = U ./ sum(U);

        //    %compute gradients
        //    %one for mean
        //    mean_grad = U*S';

        //    %variance gradient
        //    S_sq_minus = S.^2 - 1;
        //    var_grad = U * S_sq_minus';

        //    %update parameters
        //    mean_vec = mean_vec + learn_rates(1) * var_vec .* mean_grad';

        //    var_vec = var_vec .* exp(learn_rates(2) / 2 * var_grad)';

        //    %evaluate fitness of mean (for plotting)
        //    mean_vec_fit(i) = fitness(mean_vec);

        //    %uncomment for spiffy updating plot
        //    %plot(mean_vec_fit)
        //    %drawnow
        //end


        public sNES(double[] lb, double[] ub, bool[] xint, int itermax, Func<double[], double> evalfnc, int seed, Dictionary<string, object> settings, double[][] x0 = null)
            : base(lb, ub, xint, itermax, evalfnc, seed)
        {

        }



        public override void solve()
        {
            throw new NotImplementedException();
        }

        protected override void storeCurrentBest()
        {
            throw new NotImplementedException();
        }

        protected override bool CheckIfNaN(double fxtest)
        {
            throw new NotImplementedException();
        }
    }



    
    /// <summary>
    /// options: deterministic or stochastic (metropolis) acceptance rules
    ///  Moscato & Fontanari 1990
    /// implement re-annealing, with x* of previous run as new x0
    /// </summary>
    public class SimpleSA : SO_Solver
    {

        public override void solve()
        {
            throw new NotImplementedException();
        }
    }




    /// <summary>
    /// Differential Evolution.
    /// Storn and Price (1997). Differential Evolution - A Simple and Efficient Heuristic for Global Optimization over Contiuous Spaces.
    /// </summary>
    public class SimpleDE : SO_Solver
    {

        public override void solve()
        {
            throw new NotImplementedException();
        }
    }









    /// <summary>
    /// with probabilistic restart...
    /// Luersen & Le Riche (2004). Globalized Nelder-Mead method for engineering optimization.
    /// </summary>
    public class NelderMead : SO_Solver
    {

        public override void solve()
        {
            throw new NotImplementedException();
        }
    }
    */
}
