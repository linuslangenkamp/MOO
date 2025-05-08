#include <iostream>
#include <interfaces/openmodelica/main_opt.h>

// for xml, bin, json run in 
// Projects/Optimization/src/modelica/example$ ../../../build/src/modelica/example/include_test_cmk

/* entry point to the optimization runtime from OpenModelica generated code
 * this dir, i.e. interfaces/openmodelica, defines the glue code (Mesh, Problem, Flags, CallSimulation) between the runtime and the simulation code */
int _main_OptimitationRuntime() {
    printf("Entry Point to Optimization Runtime\n");
    //printf("START, STOP, STEPS, STEPSIZE: %f, %f, %lu, %f\n", data->simulationInfo->startTime, data->simulationInfo->stopTime, data->simulationInfo->stepSize);

    // create initial mesh, integrator, problem
    // fs, bs, set sparsity pattern and data
    // create solver => run()

    return 0;
}
