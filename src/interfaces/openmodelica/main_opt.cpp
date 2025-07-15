// for xml, bin, json run in
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/collocation.h>
#include <base/log.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/ipopt_solver.h>
#include <nlp/instances/gdop/gdop.h>

#include "gdop_problem.h"

using namespace OpenModelica;

// TODO: rename all variables on this OpenModelica side, clearly we use stuff from OPT (choose camelCase or _!)
// TODO: add a simple optional plotting lib (open source)

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    LOG_PREFIX('*', "Entry point [OPT] - _main_OptimitationRuntime\n");

    /* create info struct <-> same purpose as DATA* in OpenModeica */
    auto info = InfoGDOP(data, threadData, argc, argv);
    auto nlp_solver_flags = NLP::NLPSolverFlags(argc, argv);
    info.set_omc_flags(nlp_solver_flags);
    nlp_solver_flags.print();

    auto collocation = Collocation();
    auto mesh = Mesh::create_equidistant_fixed_stages(info.tf, info.intervals, info.stages, collocation);
    auto problem = create_gdop(info, mesh, collocation);

    auto strategies = std::make_unique<GDOP::Strategies>(default_strategies(info, S_DASSL));
    auto gdop = GDOP::GDOP(problem, collocation, mesh, std::move(strategies));

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_flags);
    ipopt_solver.optimize();

    // TODO: add verification step in Solver - costates or resimulate with simulation runtime
    auto optimal_control = gdop.optimal_solution->copy_extract_controls();
    auto simulated_optimum = gdop.strategies->simulate(gdop, optimal_control, info.intervals, 0.0, info.tf, gdop.get_curr_x_x0());
    simulated_optimum->print();
    gdop.strategies->emitter = std::make_shared<GDOP::CSVEmitter>(GDOP::CSVEmitter("opt_out.csv"));
    gdop.strategies->emit(gdop, *gdop.optimal_solution);

    gdop.strategies->verify(gdop, *gdop.optimal_solution);

    return 0;
}
