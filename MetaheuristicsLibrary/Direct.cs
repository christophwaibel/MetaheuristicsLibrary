using System;
using System.Collections.Generic;
using System.Linq;


namespace MetaheuristicsLibrary.SingleObjective
{
    /// <summary>
    /// DIRECT.
    /// <para/>DIviding RECTangles algorithm.
    /// <para/>Jones, Perttunen, Stuckman (1993). Lipschitzian Optimization Without the Lipschitz Constant.
    /// <para/>Based on Matlab code: http://www4.ncsu.edu/~ctk/Finkel_Direct/
    /// </summary>
    public class Direct : SingleObjective
    {

        public bool[] xInteger { get; private set; }    //indicating, which variable is an Integer. 0:double, 1:Integer


        // parameters DIRECT
        double par_ep = 1e-4;
        bool par_testflag = false;                      // terminate if within a relative tolerence of f_opt?
        double par_globalmin = 0;                       // minimum value of function. if known.



        double[] thirds = new double[100];              // thirds of the hypercube. just an array storing the side lengths of thirds of thirds of thirds...
        List<int>[] lengths;                            // # of slices per dimension-direction. one slice means division into 3, or 2 cuts. so 1 means the hypercube is cut into 3 parts in that dimension
        List<double> szes = new List<double>();         // size of each hyperrectangle. its the normvector ||x||, x∈{x1,x2,x3,...}. each vector is measured perpendicular from center of hyperrectangle to its edge
        List<double>[] c;                               // center coordinates of all hyperrectangles. array of list, because each element per variable dimension. list then extended if new boxes are added
        List<double> fc = new List<double>();           // function values at each box
        double minval;                                  // minimum function value so far
        double[] xatmin;                                // x of the best (minval) so far
        double perror;                                  // error... difference to expected minimum




        List<Tuple<int, int, double>> history = new List<Tuple<int, int, double>>();       // history. item1:  iteration, item2: function evaluations, item3: min value at that iteration
        Dictionary<int, double> S;                      // Pareto optimal solutions per iteration. int=> index of box, double=>size of box (sorted)


        double fxAbort;     //if NaN, then abort


        public Direct(double[] lb, double[] ub, bool[] xint, int evalmax, Func<double[], double> evalfnc, int seed)
            : base(lb, ub, xint, evalmax, evalfnc, seed)
        {

        }


        private void initialize()
        {
            this.DIRini();
            if (this.CheckIfNaN(this.fxAbort)) return;
            this.S = new Dictionary<int, double>();  // purge List of potential boxes. 1 is index (make (int)...), 2 is szes size of box
            this.find_po();          // find pareto optimal boxes... dictionary S
        }


        public override void solve()
        {
            this.initialize();
            if (this.CheckIfNaN(this.fxAbort)) return;

            while (base.evalcount < base.evalmax)
            {
                // %-- Loop through the potentially optimal hrectangles -----------%
                // %-- and divide -------------------------------------------------%
                foreach (KeyValuePair<int, double> potentialBox in this.S)
                {
                    this.DIRdivide(potentialBox.Key);
                    if (this.CheckIfNaN(this.fxAbort)) return;
                }

                this.DIRbest();
                S = new Dictionary<int, double>();  // purge List of potential boxes. 1 is index (make (int)...), 2 is szes size of box
                find_po();          // find pareto optimal boxes... dictionary S


            }
        }

        protected override void storeCurrentBest()
        {
            // being taken care of in DIRbest()
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



        void DIRbest()
        {
            //            %-- update minval, xatmin --------------------------------------%
            double min = fc[0];
            int minIndex = 0;

            for (int i = 1; i < fc.Count; i++)
            {
                if (fc[i] < min)
                {
                    min = fc[i];
                    minIndex = i;
                }
            }
            minval = fc[minIndex];

            for (int u = 0; u < base.n; u++)
            {
                xatmin[u] = c[u][minIndex];
            }


            base.xopt = xatmin;
            base.fxopt = fc[minIndex];
        }



        /// <summary>
        /// Divides rectangle i that is passed in
        /// </summary>
        void DIRdivide(int index)
        {
            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 1. Determine which sides are the largest
            int[] li = new int[base.n];
            for (int i = 0; i < li.Length; i++)
            {
                li[i] = lengths[i][index];
            }
            int biggy = li.Min();
            List<int> ls = new List<int>();
            for (int i = 0; i < li.Length; i++)
            {
                if (li[i] == biggy) ls.Add(i);
            }
            int lssize = ls.Count;
            //////////////////////////////////////////////////////////////////////////////////////



            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 2. Evaluate function in directions of biggest size
            //%--    to determine which direction to make divisions
            double[] oldc = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                oldc[i] = c[i][index];
            }

            double delta = thirds[biggy];


            double[][] newc_left = new double[lssize][];
            double[][] newc_right = new double[lssize][];
            for (int i = 0; i < lssize; i++)
            {
                newc_left[i] = new double[base.n];
                newc_right[i] = new double[base.n];
            }

            for (int i = 0; i < base.n; i++)
            {
                for (int u = 0; u < lssize; u++)
                {
                    newc_left[u][i] = oldc[i];
                    newc_right[u][i] = oldc[i];
                }
            }



            double[] f_left = new double[lssize + 1];
            double[] f_right = new double[lssize + 1];


            int oldcounter = base.evalcount;
            for (int i = 0; i < lssize; i++)
            {
                int lsi = ls[i];
                newc_left[i][lsi] = newc_left[i][lsi] - delta;
                newc_right[i][lsi] = newc_right[i][lsi] + delta;
                f_left[i] = base.evalfnc(ReverseNormalization(newc_left[i]));
                Console.WriteLine("fx,x1,x2,{0},{1},{2}", f_left[i], ReverseNormalization(newc_left[i])[0], ReverseNormalization(newc_left[i])[1]);
                if (this.CheckIfNaN(f_left[i]))
                {
                    this.fxAbort = double.NaN;
                    return;
                }
                f_right[i] = base.evalfnc(ReverseNormalization(newc_right[i]));
                Console.WriteLine("fx,x1,x2,{0},{1},{2}", f_right[i], ReverseNormalization(newc_right[i])[0], ReverseNormalization(newc_right[i])[1]);
                if (this.CheckIfNaN(f_right[i]))
                {
                    this.fxAbort = double.NaN;
                    return;
                }
                base.evalcount += 2;
            }
            //////////////////////////////////////////////////////////////////////////////////////


            DIRsort(f_left, f_right, ls.ToArray(), oldcounter, index, newc_left, newc_right);
        }


        void DIRsort(double[] f_left, double[] f_right, int[] ls, int oldcounter, int index, double[][] newc_left, double[][] newc_right)
        {
            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 3. Sort w for division order
            double[,] w = new double[ls.Length, 2];
            for (int i = 0; i < ls.Length; i++)
            {
                w[i, 1] = ls[i];
                double min_f;
                if (f_left[i] < f_right[i])
                {
                    min_f = f_left[i];
                }
                else
                {
                    min_f = f_right[i];
                }
                w[i, 0] = min_f;
            }

            double[,] V = new double[ls.Length, 2];
            int[,] order = new int[ls.Length, 2];

            int[] order1 = new int[ls.Length];
            int[] order2 = new int[ls.Length];
            double[] dim1_f = new double[ls.Length];
            int[] dim2_ls = new int[ls.Length];
            for (int i = 0; i < ls.Length; i++)
            {
                dim1_f[i] = w[i, 0];
                dim2_ls[i] = (int)w[i, 1];
                order1[i] = i;
                order2[i] = i;
            }
            Array.Sort(dim1_f, order1);
            Array.Sort(dim2_ls, order2);

            for (int i = 0; i < ls.Length; i++)
            {
                V[i, 0] = dim1_f[i];
                V[i, 1] = dim2_ls[i];
                order[i, 0] = order1[i];
                order[i, 1] = order2[i];
            }
            //////////////////////////////////////////////////////////////////////////////////////




            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 4. Make divisions in order specified by order
            for (int i = 0; i < order.GetLength(0); i++)
            {
                int newleftindex = oldcounter + 2 * i;
                int newrightindex = oldcounter + 2 * i + 1;

                //   %-- 4.1 create new rectangles identical to the old one
                int[] oldrect = new int[base.n];
                for (int u = 0; u < base.n; u++)
                {
                    int copy = lengths[u][index];
                    oldrect[u] = copy;
                }

                for (int u = 0; u < base.n; u++)
                {
                    lengths[u].Add(0);
                    lengths[u].Add(0);
                    int copy = oldrect[u];
                    lengths[u][newleftindex] = copy;
                    lengths[u][newrightindex] = copy;
                }

                //   %-- old, and new rectangles have been sliced in order(i) direction
                lengths[ls[order[i, 1]]][newleftindex] = lengths[ls[order[i, 1]]][index] + 1;
                lengths[ls[order[i, 1]]][newrightindex] = lengths[ls[order[i, 1]]][index] + 1;
                lengths[ls[order[i, 1]]][index] = lengths[ls[order[i, 1]]][index] + 1;

                //   %-- add new columns to c
                for (int u = 0; u < base.n; u++)
                {
                    c[u].Add(0);
                    c[u].Add(0);
                    double copyleft = newc_left[order[i, 0]][u];
                    double copyright = newc_right[order[i, 0]][u];
                    c[u][newleftindex] = copyleft;
                    c[u][newrightindex] = copyright;
                }

                //   %-- add new values to fc
                fc.Add(f_left[order[i, 0]]);
                fc.Add(f_right[order[i, 0]]);

                //   %-- add new values to con
                //   %-- add new flag values to feas_flags
                //   %-- store sizes of each rectangle
                double[] xleft = new double[base.n];
                double[] xright = new double[base.n];
                for (int u = 0; u < base.n; u++)
                {
                    xleft[u] = Math.Pow((double)1 / 3, lengths[u][newleftindex]);
                    xright[u] = Math.Pow((double)1 / 3, lengths[u][newrightindex]);
                }
                double szesleft = 0.5 * Misc.Vector.Norm(xleft);
                double szesright = 0.5 * Misc.Vector.Norm(xright);
                szes.Add(szesleft);
                szes.Add(szesright);
            }


            double[] xindex = new double[base.n];
            for (int u = 0; u < base.n; u++)
            {
                xindex[u] = Math.Pow((double)1 / 3, lengths[u][index]);
            }
            double szesindex = 0.5 * Misc.Vector.Norm(xindex);
            szes[index] = szesindex;
            //////////////////////////////////////////////////////////////////////////////////////
        }



        /// <summary>
        /// Initialization of Direct to eliminate storing floating numbers???...
        /// </summary>
        public void DIRini()
        {
            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //1. thirds
            //%-- start by calculating the thirds array
            //%-- here we precalculate (1/3)^i which we will use frequently
            thirds[0] = (double)1 / 3;
            for (int i = 1; i < 100; i++)
            {
                thirds[i] = (1.0 / 3.0) * thirds[i - 1];
            }
            //////////////////////////////////////////////////////////////////////////////////////


            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //2. lengths
            //%-- length array will store # of slices in each dimension for        
            //%-- each rectangle. dimension will be rows; each rectangle
            //%-- will be a column.
            //             ...one slice is 2 cuts (making a box into 3 boxes)
            //%-- first rectangle is the whole unit hyperrectangle

            //4. c, center coordinates
            //%-- first element of c is the center of the unit hyperrectangle
            lengths = new List<int>[base.n];
            c = new List<double>[base.n];
            for (int i = 0; i < base.n; i++)
            {
                lengths[i] = new List<int>();
                lengths[i].Add(0);

                c[i] = new List<double>();
                c[i].Add((double)1 / 2);
            }
            //////////////////////////////////////////////////////////////////////////////////////


            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //3. Size of hyperrectangles
            //%-- store size of hyperrectangle in vector szes
            szes.Add(1);
            //////////////////////////////////////////////////////////////////////////////////////








            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //6. Eval function
            //%-- first element of f is going to be the function evaluated
            //%-- at the center of the unit hyper-rectangle.
            double[] xIn = new double[base.n];
            for (int i = 0; i < base.n; i++)
            {
                xIn[i] = c[i][0];
            }
            double fx = base.evalfnc(ReverseNormalization(xIn));
            if (this.CheckIfNaN(fx))
            {
                this.fxAbort = double.NaN;
                return;
            }
            fc.Add(fx);
            Console.WriteLine("fx,x1,x2,{0},{1},{2}", fx, ReverseNormalization(xIn)[0], ReverseNormalization(xIn)[1]);
            base.evalcount++;
            //////////////////////////////////////////////////////////////////////////////////////




            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //7. minval and xatmin
            //  that is: best so far, obj fnuc val and its x
            //            %-- initialize minval and xatmin to be center of hyper-rectangle
            xatmin = new double[base.n];
            Array.Copy(xIn, xatmin, base.n);
            minval = fc[0];

            if (par_testflag == true)
            {
                if (par_globalmin != 0)
                {
                    perror = 100 * (minval - par_globalmin) / Math.Abs(par_globalmin);
                }
                else
                {
                    perror = 100 * minval;
                }
            }
            else
            {
                perror = 2;
            }
            //////////////////////////////////////////////////////////////////////////////////////



            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //8. history
            //%-- initialize history
            history.Add(Tuple.Create(0, base.evalcount, minval));
            //////////////////////////////////////////////////////////////////////////////////////

        }


        /// <summary>
        /// Return list of PO hyperrectangles
        /// </summary>
        void find_po()
        {
            //function rects = find_po(fc,lengths,minval,ep,szes)


            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //            %-- 1. Find all rects on hub
            int[] diff_szes = new int[lengths[0].Count];

            for (int u = 0; u < diff_szes.Length; u++)
            {
                diff_szes[u] = 0;
                for (int i = 0; i < base.n; i++)
                {
                    diff_szes[u] = diff_szes[u] + lengths[i][u];
                }
            }

            int tmp_max = diff_szes.Max();

            int j = 0;

            int[] sum_lengths = new int[diff_szes.Length];
            Array.Copy(diff_szes, sum_lengths, diff_szes.Length);

            List<int> hull = new List<int>();

            for (int i = 0; i <= tmp_max; i++)
            {
                List<int> tmp_idx = new List<int>();
                for (int u = 0; u < sum_lengths.Length; u++)
                {
                    if (sum_lengths[u] == i) tmp_idx.Add(u);
                }

                if (tmp_idx.Count > 0)
                {
                    double tmp_n = 0;
                    int hullidx = 0;
                    for (int u = 0; u < tmp_idx.Count; u++)
                    {
                        if (u == 0)
                        {
                            tmp_n = fc[tmp_idx[u]];
                            hullidx = tmp_idx[u];
                        }
                        else if (fc[tmp_idx[u]] < tmp_n)
                        {
                            tmp_n = fc[tmp_idx[u]];
                            hullidx = tmp_idx[u];
                        }
                    }

                    hull.Add(hullidx);
                    j++;

                    //        %-- 1.5 Check for ties
                    List<int> ties = new List<int>();
                    for (int u = 0; u < tmp_idx.Count; u++)
                    {
                        if (Math.Abs(fc[tmp_idx[u]] - tmp_n) <= 1e-13) ties.Add(u);
                    }

                    if (ties.Count > 1)
                    {
                        List<int> mod_ties = new List<int>();
                        for (int u = 0; u < ties.Count; u++)
                        {
                            if (tmp_idx[ties[u]] != hull[j - 1]) mod_ties.Add(u);
                        }

                        for (int u = 0; u < mod_ties.Count; u++)
                        {
                            hull.Add(tmp_idx[ties[mod_ties[u]]]);
                            j++;
                        }
                    }
                }
            }
            //////////////////////////////////////////////////////////////////////////////////////



            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //          %-- 2. Compute lb and ub for rects on hub
            double[] lbound = calc_lbound(hull);
            double[] ubound = calc_ubound(hull);
            //////////////////////////////////////////////////////////////////////////////////////


            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 3. Find indeces of hull who satisfy
            //%--    1st condition
            List<int> maybe_po = new List<int>();
            for (int i = 0; i < lbound.Length; i++)
            {
                if (lbound[i] - ubound[i] <= 0) maybe_po.Add(i);
            }
            //////////////////////////////////////////////////////////////////////////////////////




            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 4. Find indeces of hull who satisfy
            //%--    2nd condition
            int t_len = maybe_po.Count;

            List<int> po = new List<int>();
            if (minval != 0)
            {
                for (int i = 0; i < t_len; i++)
                {
                    if (((minval - fc[hull[maybe_po[i]]]) / Math.Abs(minval) +
                        szes[hull[maybe_po[i]]] * ubound[maybe_po[i]] / Math.Abs(minval))
                        >= par_ep)
                    {
                        po.Add(i);
                    }
                }
            }
            else
            {
                for (int i = 0; i < t_len; i++)
                {
                    if ((fc[hull[maybe_po[i]]] - szes[hull[maybe_po[i]]] * ubound[maybe_po[i]])
                       <= 0)
                    {
                        po.Add(i);
                    }
                }
            }

            List<int> final_pos = new List<int>();
            for (int i = 0; i < po.Count; i++)
            {
                final_pos.Add(hull[maybe_po[po[i]]]);
            }
            //////////////////////////////////////////////////////////////////////////////////////





            //____________________________________________________________________________________
            //////////////////////////////////////////////////////////////////////////////////////
            //%-- 6. Return dictionary with potential boxes
            for (int i = 0; i < final_pos.Count; i++)
            {
                S.Add(final_pos[i], szes[final_pos[i]]);
            }
            //////////////////////////////////////////////////////////////////////////////////////
        }




        /// <summary>
        /// Calculate the lbound used in determing potentially optimal hrectangles.
        /// </summary>
        double[] calc_lbound(List<int> hull)
        {
            int hull_length = hull.Count;
            double[] lbDir = new double[hull_length];

            int[,] hull_lengths = new int[lengths.Length, hull_length];
            for (int i = 0; i < lengths.Length; i++)
            {
                for (int u = 0; u < hull_length; u++)
                {
                    hull_lengths[i, u] = lengths[i][hull[u]];
                }
            }


            for (int i = 0; i < hull_length; i++)
            {
                List<int> tmp_rects = new List<int>();
                int[] sum_hull_lengths = new int[hull_lengths.GetLength(1)];
                for (int u = 0; u < hull_lengths.GetLength(1); u++)
                {
                    sum_hull_lengths[u] = 0;
                    for (int k = 0; k < hull_lengths.GetLength(0); k++)
                    {
                        sum_hull_lengths[u] += hull_lengths[k, u];
                    }
                }
                int sum_lengths = 0;
                for (int k = 0; k < lengths.Length; k++)
                {
                    sum_lengths += lengths[k][hull[i]];
                }

                for (int u = 0; u < sum_hull_lengths.Length; u++)
                {
                    if (sum_hull_lengths[u] > sum_lengths)
                    {
                        tmp_rects.Add(u);
                    }
                }


                if (tmp_rects.Count > 0)
                {
                    List<double> tmp_f = new List<double>();
                    List<double> tmp_szes = new List<double>();
                    List<double> tmp_lbs = new List<double>();
                    for (int u = 0; u < tmp_rects.Count; u++)
                    {
                        tmp_f.Add(fc[hull[tmp_rects[u]]]);
                        tmp_szes.Add(szes[hull[tmp_rects[u]]]);
                        tmp_lbs.Add((fc[hull[i]] - tmp_f[u]) / (szes[hull[i]] - tmp_szes[u]));
                    }
                    lbDir[i] = tmp_lbs.Max();
                }
                else
                {
                    lbDir[i] = -1.976e14;
                }

            }
            return lbDir;
        }


        /// <summary>
        /// Calculate the ubound used in determing potentially optimal hrectangles.
        /// </summary>
        double[] calc_ubound(List<int> hull)
        {
            int hull_length = hull.Count;
            double[] ubDir = new double[hull_length];

            int[,] hull_lengths = new int[lengths.Length, hull_length];
            for (int i = 0; i < lengths.Length; i++)
            {
                for (int u = 0; u < hull_length; u++)
                {
                    hull_lengths[i, u] = lengths[i][hull[u]];
                }
            }

            for (int i = 0; i < hull_length; i++)
            {
                List<int> tmp_rects = new List<int>();
                int[] sum_hull_lengths = new int[hull_lengths.GetLength(1)];
                for (int u = 0; u < hull_lengths.GetLength(1); u++)
                {
                    sum_hull_lengths[u] = 0;
                    for (int k = 0; k < hull_lengths.GetLength(0); k++)
                    {
                        sum_hull_lengths[u] += hull_lengths[k, u];
                    }
                }
                int sum_lengths = 0;
                for (int k = 0; k < lengths.Length; k++)
                {
                    sum_lengths += lengths[k][hull[i]];
                }

                for (int u = 0; u < sum_hull_lengths.Length; u++)
                {
                    if (sum_hull_lengths[u] < sum_lengths)
                    {
                        tmp_rects.Add(u);
                    }
                }


                if (tmp_rects.Count > 0)
                {
                    List<double> tmp_f = new List<double>();
                    List<double> tmp_szes = new List<double>();
                    List<double> tmp_ubs = new List<double>();
                    for (int u = 0; u < tmp_rects.Count; u++)
                    {
                        tmp_f.Add(fc[hull[tmp_rects[u]]]);
                        tmp_szes.Add(szes[hull[tmp_rects[u]]]);
                        tmp_ubs.Add((tmp_f[u] - fc[hull[i]]) / (tmp_szes[u] - szes[hull[i]]));
                    }
                    ubDir[i] = tmp_ubs.Min();
                }
                else
                {
                    ubDir[i] = 1.976e14;
                }

            }
            return ubDir;
        }




        private double[] NormalizeX(double[] _x)
        {
            double[] _xNorm = new double[_x.Length];
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
