// for xml, bin, json run in 
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/collocation.h>
#include <base/mesh.h>

#include <nlp/solvers/ipopt/ipopt_solver.h>
#include <nlp/instances/gdop/gdop.h>

#include "gdop_problem.h"

// TODO: warp all this into a namespace OpenModelica or so

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    printf("Entry point - _main_OptimitationRuntime\n\n");

    /* create info struct <-> same purpose as DATA* in OpenModeica */
    InfoGDOP info;

    // TODO: add flag to set this degree
    // stages = atoi((char*)omc_flagValue[FLAG_OPTIMIZER_NP]); // but please rename this flag to FLAG_OPT_STAGES or so

    // TODO: fix potential start / stop time offset, e.g. startTime = 1 => in callback make offset +=1 for t
    int stages = 3;
    F64 tf = data->simulationInfo->stopTime - data->simulationInfo->startTime;
    int intervals = (int)(round(tf/data->simulationInfo->stepSize));
    auto fLGR = std::make_unique<Collocation>();
    auto mesh = std::make_unique<Mesh>(Mesh::create_equidistant_fixed_stages(tf, intervals, stages, *fLGR));
    auto problem = create_gdop(data, threadData, info, *mesh);

    Trajectory initial_guess({0, tf}, {{1, 0}}, {{0, 0}}, {}, InterpolationMethod::LINEAR);
    GDOP gdop(problem, *fLGR, *mesh, initial_guess);
    
    std::unordered_map<std::string, std::string> nlp_solver_flags;
    nlp_solver_flags["Hessian"] = "LBFGS";

    IpoptSolver ipopt_solver(gdop, nlp_solver_flags);
    ipopt_solver.optimize();

    return 0;
}
