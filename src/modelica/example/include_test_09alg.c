/* Algebraic */
#include "include_test_model.h"

#ifdef __cplusplus
extern "C" {
#endif

/* forwarded equations */
extern void include_test_eqFunction_12(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_13(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_16(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_17(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_19(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_20(DATA* data, threadData_t *threadData);

static void functionAlg_system0(DATA *data, threadData_t *threadData)
{
  int id;

  static void (*const eqFunctions[6])(DATA*, threadData_t*) = {
    include_test_eqFunction_12,
    include_test_eqFunction_13,
    include_test_eqFunction_16,
    include_test_eqFunction_17,
    include_test_eqFunction_19,
    include_test_eqFunction_20
  };
  
  static const int eqIndices[6] = {
    12,
    13,
    16,
    17,
    19,
    20
  };
  
  for (id = 0; id < 6; id++) {
    eqFunctions[id](data, threadData);
    threadData->lastEquationSolved = eqIndices[id];
  }
}
/* for continuous time variables */
int include_test_functionAlgebraics(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

#if !defined(OMC_MINIMAL_RUNTIME)
  if (measure_time_flag) rt_tick(SIM_TIMER_ALGEBRAICS);
#endif
  data->simulationInfo->callStatistics.functionAlgebraics++;

  include_test_function_savePreSynchronous(data, threadData);
  
  functionAlg_system0(data, threadData);

#if !defined(OMC_MINIMAL_RUNTIME)
  if (measure_time_flag) rt_accumulate(SIM_TIMER_ALGEBRAICS);
#endif

  TRACE_POP
  return 0;
}

#ifdef __cplusplus
}
#endif
