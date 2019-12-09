using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using MetaheuristicRepository.Solvers_SO;

namespace MetaheuristicsLibrary
{
    public static class HybridAlgorithms
    {
        public static Hybrid.settings PSO(int hybridMaxFncCalls, int hybridpopsize)
        {

            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };        //3 operators, so add two parameter vectors
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVE
            settings._Moving = new string[] { "M_CopyXtoXTest", "M_AddStep" };           //2 operators, so 2 parameter vector
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 1 });                            //always first element of lambda for a move type 1 operator must be 1 or 0, : replacing operator (1), or if it creates new x' (0) 

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { " ", "M__SigmaPSOsteps" };         //first one empty, because copyXtoX' doesnt need stepsize
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[4] { 0.4, 0.9, 0.9, 0.2 });           // 1: weight/inertia; 2: bias to own best; 3: bias to global best; 4: max v for first velocity at t=0

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest", "C_BestSolutionOfTrajectory" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AcceptEntireTestPop", "S_MOVE_AcceptIncrementalIndex" };  //since only one move op of type 1, we dont need move op selection mechanism
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });


            return settings;
        }

        public static Hybrid.settings SA(int hybridMaxFncCalls)
        {
            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //2. MOVE
            settings._Moving = new string[] { "M_PerturbGaussDiminProbs" };
            settings._parameters_t0.Add(new double[2] { 0, 0.5 });                          //1: adding (0). 2: probability p

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { "M__SigmaSimAnn" };
            settings._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)

            //3. EVALUATION
            settings._Evaluation = new string[] { };

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_SimAnn", "S_MOVE_AcceptIncrementalIndex" };
            settings._parameters_t0.Add(new double[1] { 10 });                            //multiplication of delta E
            settings._parameters_t0.Add(new double[] { });

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });


            return settings;
        }

        public static Hybrid.settings ES(int hybridMaxFncCalls, int hybridpopsize)
        {
            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            ////2. MOVE
            //settings._Moving = new string[] { "M_PerturbGaussDiminProbs","M_InterRecomb", "M_CopyXtoXTest"  };
            //settings._parameters_t0.Add(new double[2] { 0, 0.5 });                          //1: adding (0). 2: probability p
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize
            //settings._Moving_stepsize = new string[] { "M__SigmaSimAnn", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)
            //settings._parameters_t0.Add(new double[1] { 0.05 });

            ////2. MOVE
            //settings._Moving = new string[] { "M_MutationPohlheim", "M_InterRecomb", "M_CopyXtoXTest" };
            //settings._parameters_t0.Add(new double[3] { 0, 0.1, 16.0 });                          //1: adding (0). 2: probability p
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize
            //settings._Moving_stepsize = new string[] { "", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[] { });
            //settings._parameters_t0.Add(new double[1] { 0.1 });    

            //2. MOVEM_MutationGaussian
            settings._Moving = new string[] { "M_MutationGaussian", "M_InterRecomb", "M_CopyXtoXTest" };
            settings._parameters_t0.Add(new double[1] { 0 });                            //adding x' 
            settings._parameters_t0.Add(new double[2] { 0, 2 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            settings._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize  M__SigmaSchwefel1977
            settings._Moving_stepsize = new string[] { "M__SigmaSchwefel1977", "M__SigmaConstant" };
            settings._parameters_t0.Add(new double[2] { 0.02, 1 });                      //0: initial sigma,   1: mutation mode (0 or 1)
            settings._parameters_t0.Add(new double[1] { 0.05 });                           //initial sigma. here: d, hyper cube extension

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });

            //S_MOVE_TournamentRecomb
            //S_MOVE_SUSRecomb
            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AccepptMuBest", "S_MOVE_TournamentRecomb", "S_MOVE_TournamentMutate" };
            settings._parameters_t0.Add(new double[1] { hybridpopsize });                //accepts pop size
            settings._parameters_t0.Add(new double[1] { 5 });                          //how many couples (each couple makes one or two kid(s).)
            settings._parameters_t0.Add(new double[1] { 5 });                           //how many to mutate (each couple makes one kid.)

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });




            return settings;
        }

        public static Hybrid.settings ES2(int hybridMaxFncCalls, int hybridpopsize)
        {
            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            ////2. MOVE
            //settings._Moving = new string[] { "M_PerturbGaussDiminProbs","M_InterRecomb", "M_CopyXtoXTest"  };
            //settings._parameters_t0.Add(new double[2] { 0, 0.5 });                          //1: adding (0). 2: probability p
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize
            //settings._Moving_stepsize = new string[] { "M__SigmaSimAnn", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)
            //settings._parameters_t0.Add(new double[1] { 0.05 });

            //2. MOVE
            settings._Moving = new string[] { "M_MutationPohlheim", "M_InterRecomb", "M_CopyXtoXTest" };
            settings._parameters_t0.Add(new double[3] { 0, 0.1, 16.0 });                          //1: adding (0). 2: probability p
            settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            settings._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { "", "M__SigmaConstant" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 0.1 });

            ////2. MOVEM_MutationGaussian
            //settings._Moving = new string[] { "M_MutationGaussian", "M_InterRecomb", "M_CopyXtoXTest" };
            //settings._parameters_t0.Add(new double[1] { 0 });                            //adding x' 
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize  M__SigmaSchwefel1977
            //settings._Moving_stepsize = new string[] { "M__SigmaSchwefel1977", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[2] { 0.02, 1 });                      //0: initial sigma,   1: mutation mode (0 or 1)
            //settings._parameters_t0.Add(new double[1] { 0.05 });                           //initial sigma. here: d, hyper cube extension

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });

            //S_MOVE_TournamentRecomb
            //S_MOVE_SUSRecomb
            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AccepptMuBest", "S_MOVE_TournamentRecomb", "S_MOVE_TournamentMutate" };
            settings._parameters_t0.Add(new double[1] { hybridpopsize });                //accepts pop size
            settings._parameters_t0.Add(new double[1] { 1 });                          //how many couples (each couple makes one or two kid(s).)
            settings._parameters_t0.Add(new double[1] { 1 });                           //how many to mutate (each couple makes one kid.)

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });




            return settings;
        }

        public static Hybrid.settings ES3(int hybridMaxFncCalls, int hybridpopsize)
        {
            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            ////2. MOVE
            //settings._Moving = new string[] { "M_PerturbGaussDiminProbs","M_InterRecomb", "M_CopyXtoXTest"  };
            //settings._parameters_t0.Add(new double[2] { 0, 0.5 });                          //1: adding (0). 2: probability p
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize
            //settings._Moving_stepsize = new string[] { "M__SigmaSimAnn", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)
            //settings._parameters_t0.Add(new double[1] { 0.05 });

            //2. MOVE
            settings._Moving = new string[] { "M_MutationPohlheim", "M_InterRecomb", "M_CopyXtoXTest" };
            settings._parameters_t0.Add(new double[3] { 0, 0.1, 16.0 });                          //1: adding (0). 2: probability p
            settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            settings._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { "", "M__SigmaConstant" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 0.1 });

            ////2. MOVEM_MutationGaussian
            //settings._Moving = new string[] { "M_MutationGaussian", "M_InterRecomb", "M_CopyXtoXTest" };
            //settings._parameters_t0.Add(new double[1] { 0 });                            //adding x' 
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize  M__SigmaSchwefel1977
            //settings._Moving_stepsize = new string[] { "M__SigmaSchwefel1977", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[2] { 0.02, 1 });                      //0: initial sigma,   1: mutation mode (0 or 1)
            //settings._parameters_t0.Add(new double[1] { 0.05 });                           //initial sigma. here: d, hyper cube extension

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });

            //S_MOVE_TournamentRecomb
            //S_MOVE_SUSRecomb
            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AccepptMuBest", "S_MOVE_SUSRecomb", "S_MOVE_TournamentMutate" };
            settings._parameters_t0.Add(new double[1] { hybridpopsize });                //accepts pop size
            settings._parameters_t0.Add(new double[1] { 1 });                          //how many couples (each couple makes one or two kid(s).)
            settings._parameters_t0.Add(new double[1] { 1 });                           //how many to mutate (each couple makes one kid.)

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });




            return settings;
        }

        public static Hybrid.settings ESSA(int hybridMaxFncCalls, int hybridpopsize)
        {
            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVE
            settings._Moving = new string[] { "M_PerturbGaussDiminProbs", "M_InterRecomb", "M_CopyXtoXTest" };
            settings._parameters_t0.Add(new double[2] { 0, 0.5 });                          //1: adding (0). 2: probability p
            settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            settings._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { "M__SigmaSimAnn", "M__SigmaConstant" };
            settings._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)
            settings._parameters_t0.Add(new double[1] { 0.05 });

            ////2. MOVE
            //settings._Moving = new string[] { "M_MutationPohlheim", "M_InterRecomb", "M_CopyXtoXTest" };
            //settings._parameters_t0.Add(new double[3] { 0, 0.1, 16.0 });                          //1: adding (0). 2: probability p
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize
            //settings._Moving_stepsize = new string[] { "", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[] { });
            //settings._parameters_t0.Add(new double[1] { 0.1 });

            ////2. MOVEM_MutationGaussian
            //settings._Moving = new string[] { "M_MutationGaussian", "M_InterRecomb", "M_CopyXtoXTest" };
            //settings._parameters_t0.Add(new double[1] { 0 });                            //adding x' 
            //settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            //settings._parameters_t0.Add(new double[] { });

            ////2.b MOVE stepsize  M__SigmaSchwefel1977
            //settings._Moving_stepsize = new string[] { "M__SigmaSchwefel1977", "M__SigmaConstant" };
            //settings._parameters_t0.Add(new double[2] { 0.02, 1 });                      //0: initial sigma,   1: mutation mode (0 or 1)
            //settings._parameters_t0.Add(new double[1] { 0.05 });                           //initial sigma. here: d, hyper cube extension

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });


            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_SimAnnPop", "S_MOVE_TournamentRecomb", "S_MOVE_TournamentMutate" };
            settings._parameters_t0.Add(new double[2] { hybridpopsize, 10 });                            //multiplication of delta E
            settings._parameters_t0.Add(new double[1] { 1 });                          //how many couples (each couple makes one or two kid(s).)
            settings._parameters_t0.Add(new double[1] { 1 });
            //S_MOVE_TournamentRecomb
            //S_MOVE_SUSRecomb
            ////5. SELECTION
            //settings._Selection = new string[] { "S_MIGR_AccepptMuBest", "S_MOVE_TournamentRecomb", "S_MOVE_TournamentMutate" };
            //settings._parameters_t0.Add(new double[1] { hybridpopsize });                //accepts pop size
            //settings._parameters_t0.Add(new double[1] { 1 });                          //how many couples (each couple makes one or two kid(s).)
            //settings._parameters_t0.Add(new double[1] { 1 });                           //how many to mutate (each couple makes one kid.)

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });




            return settings;
        }

        public static Hybrid.settings PSOES(int hybridMaxFncCalls, int hybridpopsize)
        {
            //PSO - with recombination ES.. and  popsize 25

            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };        //3 operators, so add two parameter vectors
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVE
            settings._Moving = new string[] { "M_CopyXtoXTest", "M_AddStep", "M_InterRecomb" };           //2 operators, so 2 parameter vector
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 1 });                            //always first element of lambda for a move type 1 operator must be 1 or 0, : replacing operator (1), or if it creates new x' (0) 
            settings._parameters_t0.Add(new double[2] { 1, 1 });

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { " ", "M__SigmaPSOsteps", "M__SigmaConstant" };         //first one empty, because copyXtoX' doesnt need stepsize
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[4] { 0.8, 0.7, 0.5, 0.4 });           // 1: weight/inertia; 2: bias to own best; 3: bias to global best; 4: max v for first velocity at t=0
            settings._parameters_t0.Add(new double[1] { 0.05 });

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest", "C_BestSolutionOfTrajectory" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AcceptEntireTestPop", "S_MOVE_AcceptIncrementalIndex", "S_MOVE_TournamentRecomb" };  //since only one move op of type 1, we dont need move op selection mechanism
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 1 });

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });


            return settings;
        }

        public static Hybrid.settings ESNM(int hybridMaxFncCalls, int hybridpopsize)
        {
            //ES with reflection of worst

            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVEM_MutationGaussian
            settings._Moving = new string[] { "M_MutationGaussian", "M_InterRecomb", "M_ReflectWorst", "M_CopyXtoXTest" };
            settings._parameters_t0.Add(new double[1] { 0 });                            //adding x' 
            settings._parameters_t0.Add(new double[2] { 0, 1 });                       //0: adding x'     1: how many kids per couple, 1 or 2
            settings._parameters_t0.Add(new double[1] { 0 });
            settings._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize  M__SigmaSchwefel1977
            settings._Moving_stepsize = new string[] { "M__SigmaSchwefel1977", "M__SigmaConstant", "M__SigmaConstant" };
            settings._parameters_t0.Add(new double[2] { 0.02, 1 });                      //0: initial sigma,   1: mutation mode (0 or 1)
            settings._parameters_t0.Add(new double[1] { 0.05 });                           //initial sigma. here: d, hyper cube extension
            settings._parameters_t0.Add(new double[1] { 0.05 });                           //initial sigma. here: alpha for reflecting worst

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });

            //S_MOVE_TournamentRecomb
            //S_MOVE_SUSRecomb
            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AccepptMuBest", "S_MOVE_TournamentRecomb", "S_MOVE_TournamentMutate" };
            settings._parameters_t0.Add(new double[1] { hybridpopsize });                //accepts pop size
            settings._parameters_t0.Add(new double[1] { 1 });                          //how many couples (each couple makes one or two kid(s).)
            settings._parameters_t0.Add(new double[1] { 1 });                           //how many to mutate (each couple makes one kid.)

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });




            return settings;
        }


        public static Hybrid.settings PSONM(int hybridMaxFncCalls, int hybridpopsize)
        {

            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_OrthoSimplex" };        //3 operators, so add two parameter vectors
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[1] { 0.1 });                                  // uniform random needs no parameter.

            //2. MOVE
            settings._Moving = new string[] { "M_CopyXtoXTest", "M_AddStep", "M_ReflectWorst" };           //2 operators, so 2 parameter vector
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 1 });                            //always first element of lambda for a move type 1 operator must be 1 or 0, : replacing operator (1), or if it creates new x' (0) 
            settings._parameters_t0.Add(new double[1] { 1 });

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { " ", "M__SigmaPSOsteps", "M__SigmaConstant" };         //first one empty, because copyXtoX' doesnt need stepsize
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[4] { 0.8, 0.7, 0.5, 0.4 });           // 1: weight/inertia; 2: bias to own best; 3: bias to global best; 4: max v for first velocity at t=0
            settings._parameters_t0.Add(new double[1] { 0.1 });

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest", "C_BestSolutionOfTrajectory" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AcceptEntireTestPop", "S_MOVE_AcceptIncrementalIndex" };  //since only one move op of type 1, we dont need move op selection mechanism
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });


            return settings;
        }


        public static Hybrid.settings PSOSA1(int hybridMaxFncCalls, int hybridpopsize)
        {
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////// PSO + SA selection rule     ////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Hybrid.settings settingsPSOSA = new Hybrid.settings();
            settingsPSOSA._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settingsPSOSA._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };        //3 operators, so add two parameter vectors
            settingsPSOSA._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settingsPSOSA._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settingsPSOSA._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVE
            settingsPSOSA._Moving = new string[] { "M_CopyXtoXTest", "M_AddStep", " " };           //2 operators, so 2 parameter vector
            settingsPSOSA._parameters_t0.Add(new double[] { });
            settingsPSOSA._parameters_t0.Add(new double[1] { 1 });                            //always first element of lambda for a move type 1 operator must be 1 or 0, : replacing operator (1), or if it creates new x' (0) 
            settingsPSOSA._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize
            settingsPSOSA._Moving_stepsize = new string[] { " ", "M__SigmaPSOsteps", "M__SigmaSimAnn " };         //first one empty, because copyXtoX' doesnt need stepsize
            settingsPSOSA._parameters_t0.Add(new double[] { });
            settingsPSOSA._parameters_t0.Add(new double[4] { 0.8, 0.7, 0.5, 0.4 });           // 1: weight/inertia; 2: bias to own best; 3: bias to global best; 4: max v for first velocity at t=0
            settingsPSOSA._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)

            //3. EVALUATION
            settingsPSOSA._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settingsPSOSA._Collection = new string[] { "C_GlobalBest", "C_BestSolutionOfTrajectory" };
            settingsPSOSA._parameters_t0.Add(new double[] { });
            settingsPSOSA._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settingsPSOSA._Selection = new string[] { "S_MIGR_AcceptEntireTestPop", " ", "S_MIGR_SimAnn", "S_MOVE_AcceptIncrementalIndex" };  //since only one move op of type 1, we dont need move op selection mechanism
            settingsPSOSA._parameters_t0.Add(new double[] { });
            settingsPSOSA._parameters_t0.Add(new double[] { });
            settingsPSOSA._parameters_t0.Add(new double[1] { 10 });
            settingsPSOSA._parameters_t0.Add(new double[] { });

            //6. TERMINATION
            settingsPSOSA._Termination = new string[] { "T_MaxFunctionCalls" };
            settingsPSOSA._parameters_t0.Add(new double[1] { hybridMaxFncCalls });
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            return settingsPSOSA;
        }

        public static Hybrid.settings PSOSA2(int hybridMaxFncCalls, int hybridpopsize)
        {

            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            ///////////////////////////////////////////// PSO + SA 2 :selection and move   ///////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            Hybrid.settings settingsPSOSA2 = new Hybrid.settings();
            settingsPSOSA2._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settingsPSOSA2._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };        //3 operators, so add two parameter vectors
            settingsPSOSA2._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settingsPSOSA2._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settingsPSOSA2._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVE
            settingsPSOSA2._Moving = new string[] { "M_CopyXtoXTest", "M_AddStep", "M_PerturbGaussDiminProbs" };           //2 operators, so 2 parameter vector
            settingsPSOSA2._parameters_t0.Add(new double[] { });
            settingsPSOSA2._parameters_t0.Add(new double[1] { 1 });                            //always first element of lambda for a move type 1 operator must be 1 or 0, : replacing operator (1), or if it creates new x' (0) 
            settingsPSOSA2._parameters_t0.Add(new double[2] { 0, 0.5 });

            //2.b MOVE stepsize
            settingsPSOSA2._Moving_stepsize = new string[] { " ", "M__SigmaPSOsteps", "M__SigmaSimAnn" };         //first one empty, because copyXtoX' doesnt need stepsize
            settingsPSOSA2._parameters_t0.Add(new double[] { });
            settingsPSOSA2._parameters_t0.Add(new double[4] { 0.8, 0.7, 0.5, 0.4 });           // 1: weight/inertia; 2: bias to own best; 3: bias to global best; 4: max v for first velocity at t=0
            settingsPSOSA2._parameters_t0.Add(new double[8] { 1.0, 100.0, 100.0, 0.8, 15, 0, 100, 0 });      //beta, T0, T, c, max_cnt_success, curr_cnt_success, max_cnt_total, curr_cnt_total)

            //3. EVALUATION
            settingsPSOSA2._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settingsPSOSA2._Collection = new string[] { "C_GlobalBest", "C_BestSolutionOfTrajectory" };
            settingsPSOSA2._parameters_t0.Add(new double[] { });
            settingsPSOSA2._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settingsPSOSA2._Selection = new string[] { "S_MIGR_AcceptEntireTestPop", " ", "S_MIGR_SimAnn", "S_MOVE_AcceptIncrementalIndex" };  //since only one move op of type 1, we dont need move op selection mechanism
            settingsPSOSA2._parameters_t0.Add(new double[] { });
            settingsPSOSA2._parameters_t0.Add(new double[] { });
            settingsPSOSA2._parameters_t0.Add(new double[1] { 10 });
            settingsPSOSA2._parameters_t0.Add(new double[] { });

            //6. TERMINATION
            settingsPSOSA2._Termination = new string[] { "T_MaxFunctionCalls" };
            settingsPSOSA2._parameters_t0.Add(new double[1] { hybridMaxFncCalls });
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
            //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            return settingsPSOSA2;

        }



        public static Hybrid.settings PSOHYPER(int hybridMaxFncCalls, int hybridpopsize, double inertia, double weight1, double weight2, double sigma0)
        {

            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };        //3 operators, so add two parameter vectors
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVE
            settings._Moving = new string[] { "M_CopyXtoXTest", "M_AddStep" };           //2 operators, so 2 parameter vector
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[1] { 1 });                            //always first element of lambda for a move type 1 operator must be 1 or 0, : replacing operator (1), or if it creates new x' (0) 

            //2.b MOVE stepsize
            settings._Moving_stepsize = new string[] { " ", "M__SigmaPSOsteps" };         //first one empty, because copyXtoX' doesnt need stepsize
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[4] { inertia, weight1, weight2, sigma0 });           // 1: weight/inertia; 2: bias to own best; 3: bias to global best; 4: max v for first velocity at t=0

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest", "C_BestSolutionOfTrajectory" };
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AcceptEntireTestPop", "S_MOVE_AcceptIncrementalIndex" };  //since only one move op of type 1, we dont need move op selection mechanism
            settings._parameters_t0.Add(new double[] { });
            settings._parameters_t0.Add(new double[] { });

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });


            return settings;
        }

        public static Hybrid.settings ESHYPER(int hybridMaxFncCalls, int hybridpopsize, int kidspercouple, double sigma, double sigma_d, int mutmode, double dblcouples, double dblmutants)
        {
            int mutants = Convert.ToInt32(dblmutants * hybridpopsize);
            int couples = Convert.ToInt32(dblcouples * hybridpopsize / 2);
            if (mutants < 1) mutants = 1;
            if (couples < 1) couples = 1;

            Hybrid.settings settings = new Hybrid.settings();
            settings._parameters_t0 = new List<double[]>();

            //1. INITIALISATION 
            settings._Initialization = new string[] { "I_normalizeX", "I_popSize", "I_RndUniformPop" };
            settings._parameters_t0.Add(new double[] { });                                  //"I_normalizeX", no element needed here, invoking this operator already means true
            settings._parameters_t0.Add(new double[] { hybridpopsize });                    //"I_popSize", double is round to integer 
            settings._parameters_t0.Add(new double[] { });                                  // uniform random needs no parameter.

            //2. MOVEM_MutationGaussian
            settings._Moving = new string[] { "M_MutationGaussian", "M_InterRecomb", "M_CopyXtoXTest" };
            settings._parameters_t0.Add(new double[1] { 0 });                            //adding x' 
            settings._parameters_t0.Add(new double[2] { 0, kidspercouple });                       //0: adding x'     1: how many kids per couple, 1 or 2
            settings._parameters_t0.Add(new double[] { });

            //2.b MOVE stepsize  M__SigmaSchwefel1977
            settings._Moving_stepsize = new string[] { "M__SigmaSchwefel1977", "M__SigmaConstant" };
            settings._parameters_t0.Add(new double[2] { sigma, mutmode });                      //0: initial sigma,   1: mutation mode (0 or 1)
            settings._parameters_t0.Add(new double[1] { sigma_d });                           //initial sigma. here: d, hyper cube extension

            //3. EVALUATION
            settings._Evaluation = new string[] { };                                     //empty, so no parameter vector added

            //4. COLLECTION
            settings._Collection = new string[] { "C_GlobalBest" };
            settings._parameters_t0.Add(new double[] { });

            //S_MOVE_TournamentRecomb
            //S_MOVE_SUSRecomb
            //5. SELECTION
            settings._Selection = new string[] { "S_MIGR_AccepptMuBest", "S_MOVE_TournamentRecomb", "S_MOVE_TournamentMutate" };
            settings._parameters_t0.Add(new double[1] { hybridpopsize });                //accepts pop size
            settings._parameters_t0.Add(new double[1] { couples });                          //how many couples (each couple makes one or two kid(s).)
            settings._parameters_t0.Add(new double[1] { mutants });                           //how many to mutate (each couple makes one kid.)

            //6. TERMINATION
            settings._Termination = new string[] { "T_MaxFunctionCalls" };
            settings._parameters_t0.Add(new double[1] { hybridMaxFncCalls });




            return settings;
        }

    }
}
