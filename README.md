# Optimization_Runtime

Problem Structure:
                                                      Helper
                                                        |
                    |                                   |                                                   |
external C          |                                   |                                                   |
call simulation -------- NLP specific problem  --- specific NLP #--- generic NLP --- NLP_Solver_Interface  ---# {Ipopt, Uno, ...}
callbacks           |       (also class ?)              |                                                   |
                    |                                   |                                                   |
--------------------------------------------------------|---------------------------------------------------------
                                              generic Solver / Mode
                                                        |  
                                                        |  
                                                        #
                                    {Mesh Refinement, NN Training, Events, ...}