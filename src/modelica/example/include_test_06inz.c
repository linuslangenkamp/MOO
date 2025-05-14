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
x1 = $START.x1
*/
void include_test_eqFunction_1(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,1};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) = (data->modelData->realVarsData[0] /* x1 STATE(1) */).attribute .start;
  TRACE_POP
}
extern void include_test_eqFunction_8(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_9(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_10(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_11(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_12(DATA *data, threadData_t *threadData);


/*
equation index: 7
type: SIMPLE_ASSIGN
x2 = $START.x2
*/
void include_test_eqFunction_7(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,7};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */) = (data->modelData->realVarsData[1] /* x2 STATE(1) */).attribute .start;
  TRACE_POP
}
OMC_DISABLE_OPT
void include_test_functionInitialEquations_0(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  include_test_eqFunction_1(data, threadData);
  include_test_eqFunction_8(data, threadData);
  include_test_eqFunction_9(data, threadData);
  include_test_eqFunction_10(data, threadData);
  include_test_eqFunction_11(data, threadData);
  include_test_eqFunction_12(data, threadData);
  include_test_eqFunction_7(data, threadData);
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

