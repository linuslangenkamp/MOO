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
// TODO: add a simple optional plotting lib (open source)

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    printf("Entry point [OPT] - _main_OptimitationRuntime\n\n");

    /* create info struct <-> same purpose as DATA* in OpenModeica */
    auto info = std::make_unique<InfoGDOP>(data, threadData, argc, argv); // TODO: maybe refactor: s.t. info also holds DATA, threadData_t, and flags
    auto nlp_solver_flags = std::make_unique<NLPSolverFlags>(argc, argv);
    info->set_omc_flags(*nlp_solver_flags);
    nlp_solver_flags->print();

    auto collocation = std::make_unique<Collocation>();
    auto mesh = std::make_unique<Mesh>(Mesh::create_equidistant_fixed_stages(info->tf, info->intervals, info->stages, *collocation));
    auto problem = std::make_unique<Problem>(create_gdop(*info, *mesh, *collocation));

    // TODO: add more strategies here
    auto const_trajectories = create_constant_guess(*info);
    auto const_controls = std::make_unique<ControlTrajectory>(const_trajectories->copy_extract_controls());
    auto simulated_initial_guess = simulate(*info, S_DASSL, info->intervals, *const_controls);
    simulated_initial_guess->print();

    // TODO: i dont really like this pattern, think about it
    auto gdop = std::make_unique<GDOP>(*problem, *collocation, *mesh, *simulated_initial_guess);
    auto scaling = std::make_unique<NominalScaling>(create_gdop_om_nominal_scaling(*gdop, *info));
    gdop->set_scaling(std::move(scaling));

    IpoptSolver ipopt_solver(*gdop, *nlp_solver_flags);
    ipopt_solver.optimize();

    // TODO: add verification step in Solver - costates or resimulate with simulation runtime
    auto optimal_control = gdop->optimal_solution->copy_extract_controls();
    auto simulated_optimum = simulate(*info, S_DASSL, info->intervals, optimal_control);
    simulated_optimum->print();

    emit_trajectory_om(*gdop->optimal_solution, *info);

    return 0;
}
