using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using MetaheuristicsLibrary.SolversSO;

namespace Tester
{
    public class HillclimberFROG : ISolver
    {
        /// <summary>
        /// Variable vector of final solution.
        /// </summary>
        public double[] xopt { get; private set; }
        /// <summary>
        /// Cost of final solution.
        /// </summary>
        public double fxopt { get; private set; }

        public Dictionary<string, string> settings = new Dictionary<string, string>();

        public List<string> presets = new List<string>();

        public HillclimberFROG()
        {
            settings.Add("seed", Convert.ToString(1));
            settings.Add("stepsize", Convert.ToString(0.1));
            settings.Add("itermax", Convert.ToString(150));


            presets.Add("don't know");
        }


        public bool RunSolver(List<Variable> variables, Func<List<decimal>, double> evaluate, string preset, string installFolder)
        {

            int dvar = variables.Count;
            double[] lb = new double[dvar];
            double[] ub = new double[dvar];
            bool[] xint = new bool[dvar];
            for (int i = 0; i < dvar; i++)
            {
                lb[i] = Convert.ToDouble(variables[i].LowerB);
                ub[i] = Convert.ToDouble(variables[i].UpperB);
                xint[i] = variables[i].Integer;
            }
            Func<double[], double> eval = x =>
            {
                List<decimal> decis = new List<decimal>();
                foreach (double _x in x)
                {
                    decis.Add(Convert.ToDecimal(_x));
                }
                return evaluate(decis);
            };




            try
            {
                int seed = int.Parse(settings["seed"]);
                double stepsize = double.Parse(settings["stepsize"]);
                int itermax = int.Parse(settings["itermax"]);
                Hillclimber hc = new Hillclimber(lb, ub, xint, itermax, eval, seed, stepsize);
                hc.solve();
                xopt = hc.get_Xoptimum();
                fxopt = hc.get_fxoptimum();

                return true;
            }
            catch (Exception ex)
            {
                return false;
            }

        }


        /// <summary>
        /// Get the variable vector of the final solution.
        /// </summary>
        /// <returns>Variable vector.</returns>
        public double[] get_Xoptimum()
        {
            return this.xopt;
        }

        /// <summary>
        /// Get the cost value of the final solution.
        /// </summary>
        /// <returns>Cost value.</returns>
        public double get_fxoptimum()
        {
            return this.fxopt;
        }




        public string GetErrorMessage()
        {
            return "don't know";
        }

        public IEnumerable<string> GetPresetNames()
        {
            return presets;
        }
    }
}
