// for xml, bin, json run in 
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/collocation.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/ipopt_solver.h>
#include <nlp/instances/gdop/gdop.h>

#include "gdop_problem.h"

// TODO: rename all variables on this OpenModelica side, clearly we use stuff from OPT (choose camelCase or _!)
// TODO: wrap all this into a namespace OpenModelica or so

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    printf("Entry point [OPT] - _main_OptimitationRuntime\n\n");

    /* create info struct <-> same purpose as DATA* in OpenModeica */
    auto info = std::make_unique<InfoGDOP>(); // TODO: maybe refactor: s.t. info also holds DATA, threadData_t, and flags
    auto nlp_solver_flags = std::make_unique<NLPSolverFlags>(argc, argv);
    info->set_omc_flags(data, *nlp_solver_flags);
    nlp_solver_flags->set("IpoptDerivativeTest", "false"); // debug
    nlp_solver_flags->print();

    auto collocation = std::make_unique<Collocation>();
    auto mesh = std::make_unique<Mesh>(Mesh::create_equidistant_fixed_stages(info->tf, info->intervals, info->stages, *collocation));
    auto problem = std::make_unique<Problem>(create_gdop(data, threadData, *info, *mesh, *collocation));

    // TODO: add more strategies here
    auto simulation_result = simulate(data, threadData, *info, S_DASSL, info->intervals)    ;
    //auto initial_guess = create_constant_guess(data, threadData, *info);
    GDOP gdop(*problem, *collocation, *mesh, *simulation_result);

    IpoptSolver ipopt_solver(gdop, *nlp_solver_flags);
    ipopt_solver.optimize();

    simulation_result->to_csv("simulation_moo.csv");
    gdop.optimal_solution->to_csv("optimal_moo.csv");

    // TODO: add one more stage: trajectory -> emit (set states, inputs, parameters for each t_ij, then all functionDAE each time)

    return 0;
}
