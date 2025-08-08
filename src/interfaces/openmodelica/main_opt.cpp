// for xml, bin, json run in
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>
#include <interfaces/gdopt/generated.h>

#include <base/fLGR.h>
#include <base/log.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/solver.h>
#include <nlp/instances/gdop/orchestrator.h>

#include <simulation/radau/test.h>

#include "problem.h"

using namespace OpenModelica;

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimizationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    LOG_PREFIX('*', "Entry point [MOO] _main_OptimizationRuntime()\n");

    // disable omc logs
    //memset(omc_useStream, 0, OMC_SIM_LOG_MAX * sizeof(int));

    /* create info struct <-> same purpose as DATA* in OpenModelica */
    auto info = InfoGDOP(data, threadData, argc, argv);
    auto nlp_solver_settings = NLP::NLPSolverSettings(argc, argv);
    info.set_omc_flags(nlp_solver_settings);

    nlp_solver_settings.print();
    auto mesh = Mesh::create_equidistant_fixed_stages(info.tf, info.intervals, info.stages);
    auto problem = create_gdop(info, *mesh);

    // auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    auto strategies = std::make_unique<GDOP::Strategies>(default_strategies(info));
    auto gdop = GDOP::GDOP(problem);

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_settings);

    auto orchestrator = GDOP::MeshRefinementOrchestrator(gdop, std::move(strategies), ipopt_solver);

    orchestrator.optimize();

    run_model(argc, argv);

    radau_wrapper_test();

    return 0;
}
