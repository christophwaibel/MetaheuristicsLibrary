using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Tester
{
    public static class ParetoFront
    {

        /// <summary>
        /// Non-dominated sorting. Calculates pareto front and all subsequent fronts according to objective function values.
        /// Source: A fast and elitist multiobjective genetic algorithm: NSGA-II. Deb et al (2002).
        /// </summary>
        /// <param name="x">List of variable vectors, i.e. the solutions.</param>
        /// <param name="fx">obj vals for each objective and solution candidate</param>
        /// <param name="Fronts">empty array</param>
        /// <param name="Rank">empty list</param>
        public static void S_MO_NonDominatedFront(ref List<double[]> x, ref List<double[]> fx, ref int?[] Rank, ref List<List<int>> Fronts) 
        {
            /*  Non-Dominated sort. 
            The initialized population is sorted based on non-domination. The fast
            sort algorithm [1] is described as below for each

            • for each individual p in main population P do the following
                – Initialize Sp = []. This set would contain all the individuals that is
                    being dominated by p.
                – Initialize np = 0. This would be the number of individuals that domi-
                    nate p.
                – for each individual q in P
                    * if p dominated q then
                        · add q to the set Sp i.e. Sp = Sp ? {q}
                    * else if q dominates p then
                        · increment the domination counter for p i.e. np = np + 1
                – if np = 0 i.e. no individuals dominate p then p belongs to the first
                    front; Set rank of individual p to one i.e prank = 1. Update the first
                    front set by adding p to front one i.e F1 = F1 ? {p}
            • This is carried out for all the individuals in main population P.
            • Initialize the front counter to one. i = 1
            • following is carried out while the ith front is nonempty i.e. Fi != []
                – Q = []. The set for storing the individuals for (i + 1)th front.
                – for each individual p in front Fi
                    * for each individual q in Sp (Sp is the set of individuals
                      dominated by p)
                        · nq = nq?1, decrement the domination count for individual q.
                        · if nq = 0 then none of the individuals in the subsequent
                          fronts would dominate q. Hence set qrank = i + 1. Update
                          the set Q with individual q i.e. Q = Q ? q.
                – Increment the front counter by one.
                – Now the set Q is the next front and hence Fi = Q.
    
            This algorithm is better than the original NSGA ([2]) since it utilize
            the informatoion about the set that an individual dominate (Sp) and
            number of individuals that dominate the individual (np).*/





            int mObj = fx[0].Count();



            //Rank = new int?[x.Count];

            int[] dominated_amount = new int[x.Count];     //number of individuals that dominate this individual
            List<int>[] dominate_who = new List<int>[x.Count];     //indices of individuals, this individual dominates

            int index_front = 0;
            //Fronts = new List<List<int>>();

            int[] idx = new int[x.Count];

            Fronts.Add(new List<int>());
            for (int i = 0; i < x.Count; i++)
            {
                idx[i] = i;         //storing indices for sorting later on

                dominated_amount[i] = 0;
                dominate_who[i] = new List<int>();

                for (int j = 0; j < x.Count; j++)
                {
                    int dom_less = 0;
                    int dom_equal = 0;
                    int dom_more = 0;
                    for (int m = 0; m < mObj; m++)
                    {
                        if (fx[i][m] < fx[j][m])
                        {
                            dom_less++;
                        }
                        else if (fx[i][m] == fx[j][m])
                        {
                            dom_equal++;
                        }
                        else
                        {
                            dom_more++;
                        }
                    }
                    if (dom_less == 0 && dom_equal != mObj)
                    {
                        dominated_amount[i]++;
                    }
                    else if (dom_more == 0 && dom_equal != mObj)
                    {
                        dominate_who[i].Add(j);
                    }
                }
                if (dominated_amount[i] == 0)
                {
                    Rank[i] = 0;
                    Fronts[index_front].Add(i);     //individuals on the first pareto front (rank 0)
                }
                else
                {
                    Rank[i] = null;
                }
            }



            //find the subsequent fronts
            while (Fronts[index_front].Count > 0)
            {
                List<int> Q = new List<int>();
                for (int i = 0; i < Fronts[index_front].Count; i++)
                {
                    if (dominate_who[Fronts[index_front][i]].Count > 0)      //if this individual dominates other
                    {
                        for (int j = 0; j < dominate_who[Fronts[index_front][i]].Count; j++)
                        {
                            dominated_amount[dominate_who[Fronts[index_front][i]][j]]--;
                            if (dominated_amount[dominate_who[Fronts[index_front][i]][j]] == 0)
                            {
                                Rank[dominate_who[Fronts[index_front][i]][j]] = index_front + 1;
                                Q.Add(dominate_who[Fronts[index_front][i]][j]);
                            }
                        }
                    }
                }
                index_front++;
                Fronts.Add(Q);
            }



            Array.Sort(Rank, idx);      //sorting indices according to parate fronts


            double[][] copyXTest = x.ToArray();
            double[][] copyObjVals = fx.ToArray();

            x = new List<double[]>();
            fx = new List<double[]>();

            for (int i = 0; i < idx.Length; i++)
            {
                x.Add(copyXTest[idx[i]]);
                fx.Add(copyObjVals[idx[i]]);
            }







        }


    }
}
