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

/* numerical Hessian using foward differences on OpenModelica Jacobian */
void __evalForwardDifferencesHessian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian,
                                     modelica_real h, modelica_real* jac, modelica_real* hes);

/* numerical Hessian using extrapolation on foward differences and OpenModelica Jacobian */
void __evalNumericalHessianExtrapolation(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian,
                                         modelica_real h0, int steps, modelica_real* jac, modelica_real** hes);


#endif // OPT_OM_EXTENSIONS_H
