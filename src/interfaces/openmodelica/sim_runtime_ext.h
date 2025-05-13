#ifndef OPT_OM_EXTENSIONS
#define OPT_OM_EXTENSIONS

/* has all the missing structures and functions that should be implemented in OpenModelica SimulationRuntime
 * we collect them for now
 */
#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "debug_om.h"
#include "helper.h"

/* buffer fromat of evalJacobian */
enum class new_JACOBIAN_OUTPUT_FORMAT
{
  new_JAC_OUTPUT_DENSE = 0, /* array of size #cols * #rows; fills at dense indices */
  new_JAC_OUTPUT_CSC   = 1  /* array of size nnz; fills at corresponding CSC indices */
};

/* simple extension to evalJacobian of simrt */
void new_evalJacobian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac, new_JACOBIAN_OUTPUT_FORMAT bufferFormat);

#endif // OPT_OM_EXTENSIONS
