/* Main Simulation File */

#if defined(__cplusplus)
extern "C" {
#endif

#include "include_test_model.h"
#include "simulation/solver/events.h"

/* FIXME these defines are ugly and hard to read, why not use direct function pointers instead? */
#define prefixedName_performSimulation include_test_performSimulation
#define prefixedName_updateContinuousSystem include_test_updateContinuousSystem
#include <simulation/solver/perform_simulation.c.inc>

#define prefixedName_performQSSSimulation include_test_performQSSSimulation
#include <simulation/solver/perform_qss_simulation.c.inc>


/* dummy VARINFO and FILEINFO */
const VAR_INFO dummyVAR_INFO = omc_dummyVarInfo;

int include_test_input_function(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */) = data->simulationInfo->inputVars[0];
  
  TRACE_POP
  return 0;
}

int include_test_input_function_init(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  data->simulationInfo->inputVars[0] = data->modelData->realVarsData[14].attribute.start;
  
  TRACE_POP
  return 0;
}

int include_test_input_function_updateStartValues(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  data->modelData->realVarsData[14].attribute.start = data->simulationInfo->inputVars[0];
  
  TRACE_POP
  return 0;
}

int include_test_inputNames(DATA *data, char ** names){
  TRACE_PUSH

  names[0] = (char *) data->modelData->realVarsData[14].info.name;
  
  TRACE_POP
  return 0;
}

int include_test_data_function(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  TRACE_POP
  return 0;
}

int include_test_dataReconciliationInputNames(DATA *data, char ** names){
  TRACE_PUSH

  
  TRACE_POP
  return 0;
}

int include_test_dataReconciliationUnmeasuredVariables(DATA *data, char ** names)
{
  TRACE_PUSH

  
  TRACE_POP
  return 0;
}

int include_test_output_function(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  data->simulationInfo->outputVars[0] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[4]] /* $OMC$objectLagrangeTerm variable */);
  data->simulationInfo->outputVars[1] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[5]] /* $OMC$objectMayerTerm variable */);
  data->simulationInfo->outputVars[2] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[15]] /* $con$CONSTR OPT_CONSTR */);
  data->simulationInfo->outputVars[3] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[16]] /* $finalCon$FINALCONSTR OPT_FCONSTR */);
  data->simulationInfo->outputVars[4] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */);
  data->simulationInfo->outputVars[5] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* cost_l variable */);
  data->simulationInfo->outputVars[6] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* cost_m variable */);
  
  TRACE_POP
  return 0;
}

int include_test_setc_function(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  
  TRACE_POP
  return 0;
}

int include_test_setb_function(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  
  TRACE_POP
  return 0;
}


/*
equation index: 17
type: SIMPLE_ASSIGN
cost_l = u - x2
*/
void include_test_eqFunction_17(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,17};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* cost_l variable */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */);
  TRACE_POP
}
/*
equation index: 18
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm = cost_l
*/
void include_test_eqFunction_18(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,18};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[4]] /* $OMC$objectLagrangeTerm variable */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* cost_l variable */);
  TRACE_POP
}
/*
equation index: 19
type: SIMPLE_ASSIGN
FINALCONSTR = x1 * x2
*/
void include_test_eqFunction_19(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,19};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */) = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}
/*
equation index: 20
type: SIMPLE_ASSIGN
$finalCon$FINALCONSTR = FINALCONSTR
*/
void include_test_eqFunction_20(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,20};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[16]] /* $finalCon$FINALCONSTR OPT_FCONSTR */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */);
  TRACE_POP
}
/*
equation index: 21
type: SIMPLE_ASSIGN
$con$CONSTR = x1 + u
*/
void include_test_eqFunction_21(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,21};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[15]] /* $con$CONSTR OPT_CONSTR */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) + (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  TRACE_POP
}
/*
equation index: 22
type: SIMPLE_ASSIGN
cost_m = -x1
*/
void include_test_eqFunction_22(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,22};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* cost_m variable */) = (-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}
/*
equation index: 23
type: SIMPLE_ASSIGN
$OMC$objectMayerTerm = cost_m
*/
void include_test_eqFunction_23(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,23};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[5]] /* $OMC$objectMayerTerm variable */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* cost_m variable */);
  TRACE_POP
}
/*
equation index: 24
type: SIMPLE_ASSIGN
k1 = exp(8.86 - 10215.37842190016 / u)
*/
void include_test_eqFunction_24(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,24};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */) = exp(8.86 - (DIVISION_SIM(10215.37842190016,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u",equationIndexes)));
  TRACE_POP
}
/*
equation index: 25
type: SIMPLE_ASSIGN
k2 = exp(24.25 - 18820.450885668277 / u)
*/
void include_test_eqFunction_25(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,25};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k2 variable */) = exp(24.25 - (DIVISION_SIM(18820.450885668277,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u",equationIndexes)));
  TRACE_POP
}
/*
equation index: 26
type: SIMPLE_ASSIGN
k3 = exp(23.67 - 17008.856682769725 / u)
*/
void include_test_eqFunction_26(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,26};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) = exp(23.67 - (DIVISION_SIM(17008.856682769725,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u",equationIndexes)));
  TRACE_POP
}
/*
equation index: 27
type: SIMPLE_ASSIGN
$DER.x2 = k1 * x1 + k3 * FINALCONSTR - k2 * x2
*/
void include_test_eqFunction_27(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,27};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[3]] /* der(x2) STATE_DER */) = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */)) - (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k2 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)));
  TRACE_POP
}
/*
equation index: 28
type: SIMPLE_ASSIGN
k4 = exp(18.75 - 14190.821256038647 / u)
*/
void include_test_eqFunction_28(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,28};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */) = exp(18.75 - (DIVISION_SIM(14190.821256038647,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u",equationIndexes)));
  TRACE_POP
}
/*
equation index: 29
type: SIMPLE_ASSIGN
k5 = exp(20.7 - 15599.838969404187 / u)
*/
void include_test_eqFunction_29(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,29};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */) = exp(20.7 - (DIVISION_SIM(15599.838969404187,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u",equationIndexes)));
  TRACE_POP
}
/*
equation index: 30
type: SIMPLE_ASSIGN
$DER.x1 = (((-k4) - k5 - k3) * x2 - k1) * x1
*/
void include_test_eqFunction_30(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,30};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[2]] /* der(x1) STATE_DER */) = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

OMC_DISABLE_OPT
int include_test_functionDAE(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  int equationIndexes[1] = {0};
#if !defined(OMC_MINIMAL_RUNTIME)
  if (measure_time_flag) rt_tick(SIM_TIMER_DAE);
#endif

  data->simulationInfo->needToIterate = 0;
  data->simulationInfo->discreteCall = 1;
  include_test_functionLocalKnownVars(data, threadData);
  include_test_eqFunction_17(data, threadData);

  include_test_eqFunction_18(data, threadData);

  include_test_eqFunction_19(data, threadData);

  include_test_eqFunction_20(data, threadData);

  include_test_eqFunction_21(data, threadData);

  include_test_eqFunction_22(data, threadData);

  include_test_eqFunction_23(data, threadData);

  include_test_eqFunction_24(data, threadData);

  include_test_eqFunction_25(data, threadData);

  include_test_eqFunction_26(data, threadData);

  include_test_eqFunction_27(data, threadData);

  include_test_eqFunction_28(data, threadData);

  include_test_eqFunction_29(data, threadData);

  include_test_eqFunction_30(data, threadData);
  data->simulationInfo->discreteCall = 0;
  
#if !defined(OMC_MINIMAL_RUNTIME)
  if (measure_time_flag) rt_accumulate(SIM_TIMER_DAE);
#endif
  TRACE_POP
  return 0;
}


int include_test_functionLocalKnownVars(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  
  TRACE_POP
  return 0;
}

/* forwarded equations */
extern void include_test_eqFunction_19(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_24(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_25(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_26(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_27(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_28(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_29(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_30(DATA* data, threadData_t *threadData);

static void functionODE_system0(DATA *data, threadData_t *threadData)
{
  int id;

  static void (*const eqFunctions[8])(DATA*, threadData_t*) = {
    include_test_eqFunction_19,
    include_test_eqFunction_24,
    include_test_eqFunction_25,
    include_test_eqFunction_26,
    include_test_eqFunction_27,
    include_test_eqFunction_28,
    include_test_eqFunction_29,
    include_test_eqFunction_30
  };
  
  static const int eqIndices[8] = {
    19,
    24,
    25,
    26,
    27,
    28,
    29,
    30
  };
  
  for (id = 0; id < 8; id++) {
    eqFunctions[id](data, threadData);
    threadData->lastEquationSolved = eqIndices[id];
  }
}

int include_test_functionODE(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
#if !defined(OMC_MINIMAL_RUNTIME)
  if (measure_time_flag) rt_tick(SIM_TIMER_FUNCTION_ODE);
#endif

  
  data->simulationInfo->callStatistics.functionODE++;
  
  include_test_functionLocalKnownVars(data, threadData);
  functionODE_system0(data, threadData);

#if !defined(OMC_MINIMAL_RUNTIME)
  if (measure_time_flag) rt_accumulate(SIM_TIMER_FUNCTION_ODE);
#endif

  TRACE_POP
  return 0;
}

void include_test_computeVarIndices(size_t* realIndex, size_t* integerIndex, size_t* booleanIndex, size_t* stringIndex)
{
  TRACE_PUSH

  size_t i_real = 0;
  size_t i_integer = 0;
  size_t i_boolean = 0;
  size_t i_string = 0;

  realIndex[0] = 0;
  integerIndex[0] = 0;
  booleanIndex[0] = 0;
  stringIndex[0] = 0;

  /* stateVars */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* x1 STATE(1) */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* x2 STATE(1) */
  
  /* derivativeVars */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* der(x1) STATE_DER */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* der(x2) STATE_DER */
  
  /* algVars */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* $OMC$objectLagrangeTerm variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* $OMC$objectMayerTerm variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* FINALCONSTR variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* cost_l variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* cost_m variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* k1 variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* k2 variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* k3 variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* k4 variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* k5 variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* u variable */
  
  /* discreteAlgVars */
  
  /* realOptimizeConstraintsVars */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* $con$CONSTR OPT_CONSTR */
  
  /* realOptimizeFinalConstraintsVars */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* $finalCon$FINALCONSTR OPT_FCONSTR */
  
  
  /* intAlgVars */
  
  /* boolAlgVars */
  
  /* stringAlgVars */
  
  TRACE_POP
}

/* forward the main in the simulation runtime */
extern int _main_SimulationRuntime(int argc, char**argv, DATA *data, threadData_t *threadData);

#include "include_test_12jac.h"
#include "include_test_13opt.h"

struct OpenModelicaGeneratedFunctionCallbacks include_test_callback = {
   (int (*)(DATA *, threadData_t *, void *)) include_test_performSimulation,    /* performSimulation */
   (int (*)(DATA *, threadData_t *, void *)) include_test_performQSSSimulation,    /* performQSSSimulation */
   include_test_updateContinuousSystem,    /* updateContinuousSystem */
   include_test_callExternalObjectDestructors,    /* callExternalObjectDestructors */
   NULL,    /* initialNonLinearSystem */
   NULL,    /* initialLinearSystem */
   NULL,    /* initialMixedSystem */
   #if !defined(OMC_NO_STATESELECTION)
   include_test_initializeStateSets,
   #else
   NULL,
   #endif    /* initializeStateSets */
   include_test_initializeDAEmodeData,
   include_test_computeVarIndices,
   include_test_functionODE,
   include_test_functionAlgebraics,
   include_test_functionDAE,
   include_test_functionLocalKnownVars,
   include_test_input_function,
   include_test_input_function_init,
   include_test_input_function_updateStartValues,
   include_test_data_function,
   include_test_output_function,
   include_test_setc_function,
   include_test_setb_function,
   include_test_function_storeDelayed,
   include_test_function_storeSpatialDistribution,
   include_test_function_initSpatialDistribution,
   include_test_updateBoundVariableAttributes,
   include_test_functionInitialEquations,
   1, /* useHomotopy - 0: local homotopy (equidistant lambda), 1: global homotopy (equidistant lambda), 2: new global homotopy approach (adaptive lambda), 3: new local homotopy approach (adaptive lambda)*/
   NULL,
   include_test_functionRemovedInitialEquations,
   include_test_updateBoundParameters,
   include_test_checkForAsserts,
   include_test_function_ZeroCrossingsEquations,
   include_test_function_ZeroCrossings,
   include_test_function_updateRelations,
   include_test_zeroCrossingDescription,
   include_test_relationDescription,
   include_test_function_initSample,
   include_test_INDEX_JAC_A,
   include_test_INDEX_JAC_B,
   include_test_INDEX_JAC_C,
   include_test_INDEX_JAC_D,
   include_test_INDEX_JAC_F,
   include_test_INDEX_JAC_H,
   include_test_initialAnalyticJacobianA,
   include_test_initialAnalyticJacobianB,
   include_test_initialAnalyticJacobianC,
   include_test_initialAnalyticJacobianD,
   include_test_initialAnalyticJacobianF,
   include_test_initialAnalyticJacobianH,
   include_test_functionJacA_column,
   include_test_functionJacB_column,
   include_test_functionJacC_column,
   include_test_functionJacD_column,
   include_test_functionJacF_column,
   include_test_functionJacH_column,
   include_test_linear_model_frame,
   include_test_linear_model_datarecovery_frame,
   include_test_mayer,
   include_test_lagrange,
   include_test_getInputVarIndices,
   include_test_pickUpBoundsForInputsInOptimization,
   include_test_setInputData,
   include_test_getTimeGrid,
   include_test_symbolicInlineSystem,
   include_test_function_initSynchronous,
   include_test_function_updateSynchronous,
   include_test_function_equationsSynchronous,
   include_test_inputNames,
   include_test_dataReconciliationInputNames,
   include_test_dataReconciliationUnmeasuredVariables,
   NULL,
   NULL,
   NULL,
   NULL,
   -1,
   NULL,
   NULL,
   -1
};

#define _OMC_LIT_RESOURCE_0_name_data "include_test"
#define _OMC_LIT_RESOURCE_0_dir_data "/home/linus/Projects/Optimization/src/modelica/example"
static const MMC_DEFSTRINGLIT(_OMC_LIT_RESOURCE_0_name,12,_OMC_LIT_RESOURCE_0_name_data);
static const MMC_DEFSTRINGLIT(_OMC_LIT_RESOURCE_0_dir,54,_OMC_LIT_RESOURCE_0_dir_data);

static const MMC_DEFSTRUCTLIT(_OMC_LIT_RESOURCES,2,MMC_ARRAY_TAG) {MMC_REFSTRINGLIT(_OMC_LIT_RESOURCE_0_name), MMC_REFSTRINGLIT(_OMC_LIT_RESOURCE_0_dir)}};
void include_test_setupDataStruc(DATA *data, threadData_t *threadData)
{
  assertStreamPrint(threadData,0!=data, "Error while initialize Data");
  threadData->localRoots[LOCAL_ROOT_SIMULATION_DATA] = data;
  data->callback = &include_test_callback;
  OpenModelica_updateUriMapping(threadData, MMC_REFSTRUCTLIT(_OMC_LIT_RESOURCES));
  data->modelData->modelName = "include_test";
  data->modelData->modelFilePrefix = "include_test";
  data->modelData->resultFileName = NULL;
  data->modelData->modelDir = "/home/linus/Projects/Optimization/src/modelica/example";
  data->modelData->modelGUID = "{dfd06e5b-41a4-493e-9423-9e51909e32a5}";
  #if defined(OPENMODELICA_XML_FROM_FILE_AT_RUNTIME)
  data->modelData->initXMLData = NULL;
  data->modelData->modelDataXml.infoXMLData = NULL;
  #else
  #if defined(_MSC_VER) /* handle joke compilers */
  {
  /* for MSVC we encode a string like char x[] = {'a', 'b', 'c', '\0'} */
  /* because the string constant limit is 65535 bytes */
  static const char contents_init[] =
    #include "include_test_init.c"
    ;
  static const char contents_info[] =
    #include "include_test_info.c"
    ;
    data->modelData->initXMLData = contents_init;
    data->modelData->modelDataXml.infoXMLData = contents_info;
  }
  #else /* handle real compilers */
  data->modelData->initXMLData =
  #include "include_test_init.c"
    ;
  data->modelData->modelDataXml.infoXMLData =
  #include "include_test_info.c"
    ;
  #endif /* defined(_MSC_VER) */
  #endif /* defined(OPENMODELICA_XML_FROM_FILE_AT_RUNTIME) */
  data->modelData->modelDataXml.fileName = "include_test_info.json";
  data->modelData->resourcesDir = NULL;
  data->modelData->runTestsuite = 0;
  data->modelData->nStates = 2;
  data->modelData->nVariablesRealArray = 17;
  data->modelData->nDiscreteReal = 0;
  data->modelData->nVariablesIntegerArray = 0;
  data->modelData->nVariablesBooleanArray = 0;
  data->modelData->nVariablesStringArray = 0;
  data->modelData->nParametersReal = 0;
  data->modelData->nParametersInteger = 0;
  data->modelData->nParametersBoolean = 0;
  data->modelData->nParametersString = 0;
  data->modelData->nInputVars = 1;
  data->modelData->nOutputVars = 7;
  data->modelData->nAliasReal = 1;
  data->modelData->nAliasInteger = 0;
  data->modelData->nAliasBoolean = 0;
  data->modelData->nAliasString = 0;
  data->modelData->nZeroCrossings = 0;
  data->modelData->nSamples = 0;
  data->modelData->nRelations = 0;
  data->modelData->nMathEvents = 0;
  data->modelData->nExtObjs = 0;
  data->modelData->modelDataXml.modelInfoXmlLength = 0;
  data->modelData->modelDataXml.nFunctions = 0;
  data->modelData->modelDataXml.nProfileBlocks = 0;
  data->modelData->modelDataXml.nEquations = 85;
  data->modelData->nMixedSystems = 0;
  data->modelData->nLinearSystems = 0;
  data->modelData->nNonLinearSystems = 0;
  data->modelData->nStateSets = 0;
  data->modelData->nJacobians = 6;
  data->modelData->nOptimizeConstraints = 1;
  data->modelData->nOptimizeFinalConstraints = 1;
  data->modelData->nDelayExpressions = 0;
  data->modelData->nBaseClocks = 0;
  data->modelData->nSpatialDistributions = 0;
  data->modelData->nSensitivityVars = 0;
  data->modelData->nSensitivityParamVars = 0;
  data->modelData->nSetcVars = 0;
  data->modelData->ndataReconVars = 0;
  data->modelData->nSetbVars = 0;
  data->modelData->nRelatedBoundaryConditions = 0;
  data->modelData->linearizationDumpLanguage = OMC_LINEARIZE_DUMP_LANGUAGE_MODELICA;
}

static int rml_execution_failed()
{
  fflush(NULL);
  fprintf(stderr, "Execution failed!\n");
  fflush(NULL);
  return 1;
}


#if defined(__MINGW32__) || defined(_MSC_VER)

#if !defined(_UNICODE)
#define _UNICODE
#endif
#if !defined(UNICODE)
#define UNICODE
#endif

#include <windows.h>
char** omc_fixWindowsArgv(int argc, wchar_t **wargv)
{
  char** newargv;
  /* Support for non-ASCII characters
  * Read the unicode command line arguments and translate it to char*
  */
  newargv = (char**)malloc(argc*sizeof(char*));
  for (int i = 0; i < argc; i++) {
    newargv[i] = omc_wchar_to_multibyte_str(wargv[i]);
  }
  return newargv;
}

#define OMC_MAIN wmain
#define OMC_CHAR wchar_t
#define OMC_EXPORT __declspec(dllexport) extern

#else
#define omc_fixWindowsArgv(N, A) (A)
#define OMC_MAIN main
#define OMC_CHAR char
#define OMC_EXPORT extern
#endif

#if defined(threadData)
#undef threadData
#endif
/* call the simulation runtime main from our main! */
#if defined(OMC_DLL_MAIN_DEFINE)
OMC_EXPORT int omcDllMain(int argc, OMC_CHAR **argv)
#else
int OMC_MAIN(int argc, OMC_CHAR** argv)
#endif
{
  char** newargv = omc_fixWindowsArgv(argc, argv);
  /*
    Set the error functions to be used for simulation.
    The default value for them is 'functions' version. Change it here to 'simulation' versions
  */
  omc_assert = omc_assert_simulation;
  omc_assert_withEquationIndexes = omc_assert_simulation_withEquationIndexes;

  omc_assert_warning_withEquationIndexes = omc_assert_warning_simulation_withEquationIndexes;
  omc_assert_warning = omc_assert_warning_simulation;
  omc_terminate = omc_terminate_simulation;
  omc_throw = omc_throw_simulation;

  int res;
  DATA data;
  MODEL_DATA modelData;
  SIMULATION_INFO simInfo;
  data.modelData = &modelData;
  data.simulationInfo = &simInfo;
  measure_time_flag = 0;
  compiledInDAEMode = 0;
  compiledWithSymSolver = 0;
  MMC_INIT(0);
  omc_alloc_interface.init();
  {
    MMC_TRY_TOP()
  
    MMC_TRY_STACK()
  
    include_test_setupDataStruc(&data, threadData);
    res = _main_initRuntimeAndSimulation(argc, newargv, &data, threadData);
    if(res == 0) {
      res = _main_SimulationRuntime(argc, newargv, &data, threadData);
    }
    
    MMC_ELSE()
    rml_execution_failed();
    fprintf(stderr, "Stack overflow detected and was not caught.\nSend us a bug report at https://trac.openmodelica.org/OpenModelica/newticket\n    Include the following trace:\n");
    printStacktraceMessages();
    fflush(NULL);
    return 1;
    MMC_CATCH_STACK()
    
    MMC_CATCH_TOP(return rml_execution_failed());
  }

  fflush(NULL);
  return res;
}

#ifdef __cplusplus
}
#endif


