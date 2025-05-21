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

    auto jacobian = info->exc_jac->C;
    print_jacobian_sparsity(jacobian, true, "C");
    HESSIAN_PATTERN* pattern = __generateHessianPattern(jacobian);
    info->auto_free.attach(pattern, __freeHessianPattern);

    __printHessianPattern(pattern);
    
    modelica_real* lambda = (modelica_real*)malloc(5 * sizeof(modelica_real));
    lambda[0] = 1;
    lambda[1] = 1;
    lambda[2] = 1;
    lambda[3] = 1;
    lambda[4] = 1;
    modelica_real* hes = (modelica_real*)calloc(6, sizeof(modelica_real));
    __evalHessianForwardDifferences(data, threadData, pattern, 1e-8, lambda, hes);

    print_real_var_names_values(data);
    printf("\n[Std Finite Differences]\n");
    for (int i = 0; i < 6; i++) {
        printf("hes[%d] = %.15g\n", i, hes[i]);
    }
    
    HessianFiniteDiffArgs args = {
        .data = data,
        .threadData = threadData,
        .hes_pattern = pattern,
        .lambda = lambda
    };
    
    ExtrapolationData* extrData = __initExtrapolationData(pattern->lnnz, 10);
    info->auto_free.attach(extrData, __freeExtrapolationData);
    //  1e-3, 3, 2!
    __richardsonExtrapolation(extrData, __forwardDiffHessianWrapper, &args, 1e-3, 3, 2, 1, hes);
    printf("\n[Extrapolation]\n");
    for (int i = 0; i < 6; i++) { printf("hes[%d] = %.15g\n", i, hes[i]); }

    free(lambda);
    free(hes);
    return 0;
}

