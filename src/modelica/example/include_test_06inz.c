/* Initialization */
#include "include_test_model.h"
#include "include_test_11mix.h"
#include "include_test_12jac.h"
#if defined(__cplusplus)
extern "C" {
#endif

void include_test_functionInitialEquations_0(DATA *data, threadData_t *threadData);
extern void include_test_eqFunction_13(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_14(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_15(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_16(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_17(DATA *data, threadData_t *threadData);


/*
equation index: 6
type: SIMPLE_ASSIGN
x1 = $START.x1
*/
void include_test_eqFunction_6(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,6};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) = (data->modelData->realVarsData[0] /* x1 STATE(1) */).attribute .start;
  TRACE_POP
}

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
extern void include_test_eqFunction_18(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_19(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_20(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_21(DATA *data, threadData_t *threadData);

extern void include_test_eqFunction_22(DATA *data, threadData_t *threadData);

OMC_DISABLE_OPT
void include_test_functionInitialEquations_0(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  include_test_eqFunction_13(data, threadData);
  include_test_eqFunction_14(data, threadData);
  include_test_eqFunction_15(data, threadData);
  include_test_eqFunction_16(data, threadData);
  include_test_eqFunction_17(data, threadData);
  include_test_eqFunction_6(data, threadData);
  include_test_eqFunction_7(data, threadData);
  include_test_eqFunction_18(data, threadData);
  include_test_eqFunction_19(data, threadData);
  include_test_eqFunction_20(data, threadData);
  include_test_eqFunction_21(data, threadData);
  include_test_eqFunction_22(data, threadData);
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

