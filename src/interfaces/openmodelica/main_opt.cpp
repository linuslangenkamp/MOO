// for xml, bin, json run in
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/collocation.h>
#include <base/log.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/ipopt_solver.h>
#include <nlp/instances/gdop/gdop_orchestrator.h>

#include "gdop_problem.h"

using namespace OpenModelica;

// TODO: rename all variables on this OpenModelica side, clearly we use stuff from OPT (choose camelCase or _!)
// TODO: add a simple optional plotting lib (open source)

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    LOG_PREFIX('*', "Entry point [OPT] - _main_OptimitationRuntime\n");

    // disable omc logs
    memset(omc_useStream, 0, OMC_SIM_LOG_MAX * sizeof(int));

    /* create info struct <-> same purpose as DATA* in OpenModelica */
    auto info = InfoGDOP(data, threadData, argc, argv);
    auto nlp_solver_flags = NLP::NLPSolverFlags(argc, argv);
    info.set_omc_flags(nlp_solver_flags);
    nlp_solver_flags.print();
    auto collocation = Collocation();
    auto mesh = Mesh::create_equidistant_fixed_stages(info.tf, info.intervals, info.stages, collocation);
    auto problem = create_gdop(info, mesh, collocation);

    // auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    auto strategies = std::make_unique<GDOP::Strategies>(default_strategies(info, S_GBODE));
    auto gdop = GDOP::GDOP(problem, collocation, mesh);

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_flags);

    auto orchestrator = GDOP::MeshRefinementOrchestrator(gdop, std::move(strategies), ipopt_solver);

    orchestrator.optimize();

    /*
    gdop.mesh = mesh;
    gdop.problem.resize_buffers();
    ipopt_solver.optimize();*/
    return 0;
}
