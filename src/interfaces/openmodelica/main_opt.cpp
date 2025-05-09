// for xml, bin, json run in 
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

#include <iostream>

#include <interfaces/openmodelica/main_opt.h>

#include <base/collocation.h>
#include <base/mesh.h>

#include "gdop_problem_impl.h"

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime(int argc, char** argv, DATA* data, threadData_t* threadData) {
    printf("Entry point - _main_OptimitationRuntime\n\n");

    // TODO: remove shared_ptr / unique_ptr for NLP + Mesh + Collocation? Do we really need them?

    // TODO: add flag to set this degree
    // TODO: fix potential start / stop time offset, e.g. startTime = 1 => in callback make offset +=1 for t
    // stages = atoi((char*)omc_flagValue[FLAG_OPTIMIZER_NP]); // but please rename this flag to FLAG_OPT_STAGES or so
    int stages = 3;
    F64 tf = data->simulationInfo->stopTime - data->simulationInfo->startTime;
    int intervals = (int)(round(tf/data->simulationInfo->stepSize));
    auto fLGR = std::make_unique<Collocation>();
    auto mesh = std::make_unique<Mesh>(Mesh::create_equidistant_fixed_stages(tf, intervals, stages, *fLGR));
    auto problem = create_gdop_om(data, threadData, *mesh);

    // problem
    // fs, bs, set sparsity pattern and data
    // create solver => run()
    

    return 0;
}
