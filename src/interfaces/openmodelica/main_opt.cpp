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
    printf("Entry point - _main_OptimitationRuntime\n\n");

    /* create info struct <-> same purpose as DATA* in OpenModeica */
    auto info = std::make_unique<InfoGDOP>();
    auto nlp_solver_flags = std::make_unique<NLPSolverFlags>(argc, argv);
    nlp_solver_flags->set("Hessian", "LBFGS");
    nlp_solver_flags->set("Tolerance", "1e-8");
    nlp_solver_flags->set("CPUTime", "3600");
    nlp_solver_flags->set("LinearSolver", "MA57");
    nlp_solver_flags->print();
    // TODO: add flag to set this degree
    // stages = atoi((char*)omc_flagValue[FLAG_OPTIMIZER_NP]); // but please rename this flag to FLAG_OPT_STAGES or so

    // TODO: fix potential start / stop time offset, e.g. startTime = 1 => in callback make offset +=1 for t
    int stages = 3;
    F64 tf = data->simulationInfo->stopTime - data->simulationInfo->startTime;
    int intervals = (int)(round(tf/data->simulationInfo->stepSize));
    auto fLGR = std::make_unique<Collocation>();
    auto mesh = std::make_unique<Mesh>(Mesh::create_equidistant_fixed_stages(tf, intervals, stages, *fLGR));
    auto problem = std::make_unique<Problem>(create_gdop(data, threadData, *info, *mesh));

    Trajectory initial_guess({0, tf}, {{1, 0}, {0, 1}}, {{0, 5}}, {}, InterpolationMethod::LINEAR);
    GDOP gdop(*problem, *fLGR, *mesh, initial_guess);


    IpoptSolver ipopt_solver(gdop, *nlp_solver_flags);
    ipopt_solver.optimize();

    return 0;
}


/*
print_jacobian_sparsity(om_ptr->info.exc_jac->A, true, "A");
print_jacobian_sparsity(om_ptr->info.exc_jac->B, true, "B");
print_jacobian_sparsity(om_ptr->info.exc_jac->C, true, "C");
print_jacobian_sparsity(om_ptr->info.exc_jac->D, true, "D");

problem->full->print_jacobian_sparsity_pattern();
problem->boundary->print_jacobian_sparsity_pattern();
print_bounds_fixed_vector(problem->u_bounds);*/
