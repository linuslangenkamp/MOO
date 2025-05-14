#ifndef OPT_OM_EXTENSIONS_H
#define OPT_OM_EXTENSIONS_H

/** 
 * has all the missing structures and functions that should be implemented in OpenModelica SimulationRuntime
 * we collect them for now
 */
#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "debug_om.h"

/* simple extension to evalJacobian of SimulationRuntime */
void __evalJacobian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac);

#endif // OPT_OM_EXTENSIONS_H
