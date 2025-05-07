/* Initialization */
#include "include_test_model.h"
#include "include_test_11mix.h"
#include "include_test_12jac.h"
#if defined(__cplusplus)
extern "C" {
#endif

void include_test_functionInitialEquations_0(DATA *data, threadData_t *threadData);

/*
equation index: 1
type: SIMPLE_ASSIGN
x = $START.x
*/
void include_test_eqFunction_1(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,1};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x STATE(1,x) */) = (data->modelData->realVarsData[0] /* x STATE(1,x) */).attribute .start;
  TRACE_POP
}
extern void include_test_eqFunction_3(DATA *data, threadData_t *threadData);

OMC_DISABLE_OPT
void include_test_functionInitialEquations_0(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  include_test_eqFunction_1(data, threadData);
  include_test_eqFunction_3(data, threadData);
  TRACE_POP
}

int include_test_functionInitialEquations(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  data->simulationInfo->discreteCall = 1;
  include_test_functionInitialEquations_0(data, threadData);
  data->simulationInfo->discreteCall = 0;
  
  TRACE_POP
  return 0;
}

/* No include_test_functionInitialEquations_lambda0 function */

int include_test_functionRemovedInitialEquations(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int *equationIndexes = NULL;
  double res = 0.0;

  
  TRACE_POP
  return 0;
}


#if defined(__cplusplus)
}
#endif

