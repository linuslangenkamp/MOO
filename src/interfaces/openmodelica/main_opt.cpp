// for xml, bin, json run in 
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/collocation.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/ipopt_solver.h>
#include <nlp/instances/gdop/gdop.h>

#include "gdop_problem.h"

// TODO: wrap all this into a namespace OpenModelica or so

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    printf("Entry point [OPT] - _main_OptimitationRuntime\n\n");

    /* create info struct <-> same purpose as DATA* in OpenModeica */
    auto info = std::make_unique<InfoGDOP>(); // TODO: maybe refactor: s.t. info also holds DATA, threadData_t, and flags
    auto nlp_solver_flags = std::make_unique<NLPSolverFlags>(argc, argv);
    nlp_solver_flags->set("Hessian", "LBFGS");
    nlp_solver_flags->set("Tolerance", "1e-10");
    nlp_solver_flags->set("CPUTime", "3600");
    nlp_solver_flags->set("LinearSolver", "MUMPS");
    nlp_solver_flags->print();
    // TODO: add flag to set this 1, degree
    // stages = atoi((char*)omc_flagValue[FLAG_OPTIMIZER_NP]); // but please rename this flag to FLAG_OPT_STAGES or so

    info->set_time_horizon(data, 3);
    auto collocation = std::make_unique<Collocation>();
    auto mesh = std::make_unique<Mesh>(Mesh::create_equidistant_fixed_stages(info->tf, info->intervals, info->stages, *collocation));
    auto problem = std::make_unique<Problem>(create_gdop(data, threadData, *info, *mesh, *collocation));
    auto initial_guess = std::make_unique<Trajectory>(create_constant_guess(data, threadData, *info)); // TODO: add proper strategies here

    printf("tf = %f, intervals = %d, stages = %d\n", info->tf, info->intervals, info->stages);


    GDOP gdop(*problem, *collocation, *mesh, *initial_guess);

    IpoptSolver ipopt_solver(gdop, *nlp_solver_flags);
    ipopt_solver.optimize();

    auto jac = info->exc_jac->C;
    print_jacobian_sparsity(jac, true, "C");
    HESSIAN_PATTERN* pattern = generateHessianPattern(jac);
    printHessianPattern(pattern);

    return 0;
}
