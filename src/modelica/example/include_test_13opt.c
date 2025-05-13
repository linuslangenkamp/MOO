/* Optimization */
#include "include_test_model.h"
#include "include_test_12jac.h"
#if defined(__cplusplus)
extern "C" {
#endif
/* objectiveFunction */
int include_test_mayer(DATA* data, modelica_real** res, short * index_Dres)
{
  return  -1;
}
/* objectiveIntegrand */
int include_test_lagrange(DATA* data, modelica_real** res, short * index_DresB, short *index_DresC)
{
  *res =  &(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[2]] /* $OMC$objectLagrangeTerm variable */);
  *index_DresB = 1;
  *index_DresC = 1;
  return 0;
  return -1;
}

/* fill buffer with optimization input indices */
int include_test_getInputVarIndicesInOptimization(DATA* data, int* input_var_indices)
{
  input_var_indices[0] = 4 /* u */;
  return 0;
}

/* opt vars  */
int include_test_pickUpBoundsForInputsInOptimization(DATA* data, modelica_real* min, modelica_real* max, modelica_real*nominal, modelica_boolean *useNominal, char ** name, modelica_real * start, modelica_real* startTimeOpt)
{
  min[0] = data->modelData->realVarsData[4].attribute /* u */.min;
  max[0] = data->modelData->realVarsData[4].attribute /* u */.max;
  nominal[0] = data->modelData->realVarsData[4].attribute /* u */.nominal;
  useNominal[0] = data->modelData->realVarsData[4].attribute /* u */.useNominal;
  name[0] =(char *) data->modelData->realVarsData[4].info /* u */.name;
  start[0] = data->modelData->realVarsData[4].attribute /* u */.start;
  *startTimeOpt = data->simulationInfo->startTime - 1.0;
  return 0;
}

int include_test_setInputData(DATA *data, const modelica_boolean file)
{
  TRACE_PUSH
  if(file){
  }
  data->simulationInfo->inputVars[0] = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[4]] /* u variable */);
  TRACE_POP
  return 0;
}
int include_test_getTimeGrid(DATA *data, modelica_integer * nsi, modelica_real**t){
  *nsi=(-1 );
  *t = (modelica_real*) malloc((*nsi+1)*sizeof(modelica_real));
  return 0;
}
#if defined(__cplusplus)
}
#endif