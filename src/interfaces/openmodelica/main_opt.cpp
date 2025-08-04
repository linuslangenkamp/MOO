// for xml, bin, json run in
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/fLGR.h>
#include <base/log.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/ipopt_solver.h>
#include <nlp/instances/gdop/gdop_orchestrator.h>

#include "gdop_problem.h"

using namespace OpenModelica;

// TODO: rename all variables on this OpenModelica side, clearly we use stuff from moo (choose camelCase or _!)
// TODO: add a simple optional plotting lib (open source)

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimizationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    LOG_PREFIX('*', "Entry point [OPT] - _main_OptimizationRuntime\n");

    // disable omc logs
    //memset(omc_useStream, 0, OMC_SIM_LOG_MAX * sizeof(int));

    /* create info struct <-> same purpose as DATA* in OpenModelica */
    auto info = InfoGDOP(data, threadData, argc, argv);
    auto nlp_solver_settings = NLP::NLPSolverSettings(argc, argv);
    info.set_omc_flags(nlp_solver_settings);
    nlp_solver_settings.print();
    auto mesh = Mesh::create_equidistant_fixed_stages(info.tf, info.intervals, info.stages);
    auto problem = create_gdop(info, mesh);

    // auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    auto strategies = std::make_unique<GDOP::Strategies>(default_strategies(info, S_GBODE));
    auto gdop = GDOP::GDOP(problem, mesh);

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_settings);

    auto orchestrator = GDOP::MeshRefinementOrchestrator(gdop, std::move(strategies), ipopt_solver);

    orchestrator.optimize();

    return 0;
}
