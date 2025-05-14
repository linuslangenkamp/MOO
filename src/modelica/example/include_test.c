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

  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */) = data->simulationInfo->inputVars[0];
  
  TRACE_POP
  return 0;
}

int include_test_input_function_init(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  data->simulationInfo->inputVars[0] = data->modelData->realVarsData[6].attribute.start;
  
  TRACE_POP
  return 0;
}

int include_test_input_function_updateStartValues(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  data->modelData->realVarsData[6].attribute.start = data->simulationInfo->inputVars[0];
  
  TRACE_POP
  return 0;
}

int include_test_inputNames(DATA *data, char ** names){
  TRACE_PUSH

  names[0] = (char *) data->modelData->realVarsData[6].info.name;
  
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
  data->simulationInfo->outputVars[1] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* $con$g OPT_CONSTR */);
  data->simulationInfo->outputVars[2] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[5]] /* cost_l variable */);
  
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
equation index: 8
type: SIMPLE_ASSIGN
$DER.x1 = ((-0.5) * u ^ 2.0 - u) * x1
*/
void include_test_eqFunction_8(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,8};
  modelica_real tmp0;
  tmp0 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */);
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[2]] /* der(x1) STATE_DER */) = ((-0.5) * ((tmp0 * tmp0)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}
/*
equation index: 9
type: SIMPLE_ASSIGN
$con$g = x1 + u
*/
void include_test_eqFunction_9(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,9};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* $con$g OPT_CONSTR */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) + (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */);
  TRACE_POP
}
/*
equation index: 10
type: SIMPLE_ASSIGN
cost_l = (-u) * x1
*/
void include_test_eqFunction_10(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,10};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[5]] /* cost_l variable */) = ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */))) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}
/*
equation index: 11
type: SIMPLE_ASSIGN
$DER.x2 = -cost_l
*/
void include_test_eqFunction_11(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,11};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[3]] /* der(x2) STATE_DER */) = (-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[5]] /* cost_l variable */));
  TRACE_POP
}
/*
equation index: 12
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm = cost_l
*/
void include_test_eqFunction_12(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,12};
  (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[4]] /* $OMC$objectLagrangeTerm variable */) = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[5]] /* cost_l variable */);
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
  include_test_eqFunction_8(data, threadData);

  include_test_eqFunction_9(data, threadData);

  include_test_eqFunction_10(data, threadData);

  include_test_eqFunction_11(data, threadData);

  include_test_eqFunction_12(data, threadData);
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
extern void include_test_eqFunction_8(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_10(DATA* data, threadData_t *threadData);
extern void include_test_eqFunction_11(DATA* data, threadData_t *threadData);

static void functionODE_system0(DATA *data, threadData_t *threadData)
{
  int id;

  static void (*const eqFunctions[3])(DATA*, threadData_t*) = {
    include_test_eqFunction_8,
    include_test_eqFunction_10,
    include_test_eqFunction_11
  };
  
  static const int eqIndices[3] = {
    8,
    10,
    11
  };
  
  for (id = 0; id < 3; id++) {
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
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* cost_l variable */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* u variable */
  
  /* discreteAlgVars */
  
  /* realOptimizeConstraintsVars */
  realIndex[i_real+1] = realIndex[i_real] + ((modelica_integer) 1); i_real++; /* $con$g OPT_CONSTR */
  
  /* realOptimizeFinalConstraintsVars */
  
  
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
  include_test_getInputVarIndicesInOptimization,
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
  data->modelData->modelGUID = "{1dcc975d-b1cf-41a4-abc5-d2766a3931d3}";
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
  data->modelData->nVariablesRealArray = 8;
  data->modelData->nDiscreteReal = 0;
  data->modelData->nVariablesIntegerArray = 0;
  data->modelData->nVariablesBooleanArray = 0;
  data->modelData->nVariablesStringArray = 0;
  data->modelData->nParametersReal = 0;
  data->modelData->nParametersInteger = 0;
  data->modelData->nParametersBoolean = 0;
  data->modelData->nParametersString = 0;
  data->modelData->nInputVars = 1;
  data->modelData->nOutputVars = 3;
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
  data->modelData->modelDataXml.nEquations = 27;
  data->modelData->nMixedSystems = 0;
  data->modelData->nLinearSystems = 0;
  data->modelData->nNonLinearSystems = 0;
  data->modelData->nStateSets = 0;
  data->modelData->nJacobians = 6;
  data->modelData->nOptimizeConstraints = 1;
  data->modelData->nOptimizeFinalConstraints = 0;
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


