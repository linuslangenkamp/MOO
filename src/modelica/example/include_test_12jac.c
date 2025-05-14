/* Jacobians 6 */
#include "include_test_model.h"
#include "include_test_12jac.h"
#include "simulation/jacobian_util.h"
#include "util/omc_file.h"
int include_test_functionJacH_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  TRACE_POP
  return 0;
}
int include_test_functionJacF_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  TRACE_POP
  return 0;
}
int include_test_functionJacD_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 22
type: SIMPLE_ASSIGN
cost_l.$pDERC.dummyVarC = (-u) * x1.SeedC - u.SeedC * x1
*/
void include_test_eqFunction_22(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,22};
  jacobian->tmpVars[0] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */))) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) - ((jacobian->seedVars[2] /* u.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 23
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERC.dummyVarC = cost_l.$pDERC.dummyVarC
*/
void include_test_eqFunction_23(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,23};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->tmpVars[0] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */;
  TRACE_POP
}

/*
equation index: 24
type: SIMPLE_ASSIGN
$DER.x2.$pDERC.dummyVarC = -$OMC$objectLagrangeTerm.$pDERC.dummyVarC
*/
void include_test_eqFunction_24(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,24};
  jacobian->resultVars[1] /* der(x2.$pDERC.dummyVarC) JACOBIAN_VAR */ = (-jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERC.dummyVarC JACOBIAN_VAR */);
  TRACE_POP
}

/*
equation index: 25
type: SIMPLE_ASSIGN
$con$g.$pDERC.dummyVarC = x1.SeedC + u.SeedC
*/
void include_test_eqFunction_25(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,25};
  jacobian->resultVars[3] /* $con$g.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedC SEED_VAR */ + jacobian->seedVars[2] /* u.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 26
type: SIMPLE_ASSIGN
$DER.x1.$pDERC.dummyVarC = ((-0.5) * u ^ 2.0 - u) * x1.SeedC + ((-u) * u.SeedC - u.SeedC) * x1
*/
void include_test_eqFunction_26(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,26};
  modelica_real tmp0;
  tmp0 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */);
  jacobian->resultVars[0] /* der(x1.$pDERC.dummyVarC) JACOBIAN_VAR */ = ((-0.5) * ((tmp0 * tmp0)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */))) * (jacobian->seedVars[2] /* u.SeedC SEED_VAR */) - jacobian->seedVars[2] /* u.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

OMC_DISABLE_OPT
int include_test_functionJacC_constantEqns(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_C;
  
  TRACE_POP
  return 0;
}

int include_test_functionJacC_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_C;
  include_test_eqFunction_22(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_23(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_24(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_25(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_26(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 17
type: SIMPLE_ASSIGN
cost_l.$pDERB.dummyVarB = (-u) * x1.SeedB - u.SeedB * x1
*/
void include_test_eqFunction_17(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,17};
  jacobian->tmpVars[0] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */))) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) - ((jacobian->seedVars[2] /* u.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 18
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERB.dummyVarB = cost_l.$pDERB.dummyVarB
*/
void include_test_eqFunction_18(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,18};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->tmpVars[0] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */;
  TRACE_POP
}

/*
equation index: 19
type: SIMPLE_ASSIGN
$DER.x2.$pDERB.dummyVarB = -$OMC$objectLagrangeTerm.$pDERB.dummyVarB
*/
void include_test_eqFunction_19(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,19};
  jacobian->resultVars[1] /* der(x2.$pDERB.dummyVarB) JACOBIAN_VAR */ = (-jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERB.dummyVarB JACOBIAN_VAR */);
  TRACE_POP
}

/*
equation index: 20
type: SIMPLE_ASSIGN
$con$g.$pDERB.dummyVarB = x1.SeedB + u.SeedB
*/
void include_test_eqFunction_20(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,20};
  jacobian->resultVars[3] /* $con$g.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedB SEED_VAR */ + jacobian->seedVars[2] /* u.SeedB SEED_VAR */;
  TRACE_POP
}

/*
equation index: 21
type: SIMPLE_ASSIGN
$DER.x1.$pDERB.dummyVarB = ((-0.5) * u ^ 2.0 - u) * x1.SeedB + ((-u) * u.SeedB - u.SeedB) * x1
*/
void include_test_eqFunction_21(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,21};
  modelica_real tmp1;
  tmp1 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */);
  jacobian->resultVars[0] /* der(x1.$pDERB.dummyVarB) JACOBIAN_VAR */ = ((-0.5) * ((tmp1 * tmp1)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */))) * (jacobian->seedVars[2] /* u.SeedB SEED_VAR */) - jacobian->seedVars[2] /* u.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

OMC_DISABLE_OPT
int include_test_functionJacB_constantEqns(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_B;
  
  TRACE_POP
  return 0;
}

int include_test_functionJacB_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_B;
  include_test_eqFunction_17(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_18(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_19(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_20(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_21(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 14
type: SIMPLE_ASSIGN
cost_l.$pDERA.dummyVarA = (-u) * x1.SeedA
*/
void include_test_eqFunction_14(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,14};
  jacobian->tmpVars[1] /* cost_l.$pDERA.dummyVarA JACOBIAN_TMP_VAR */ = ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */))) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */);
  TRACE_POP
}

/*
equation index: 15
type: SIMPLE_ASSIGN
$DER.x2.$pDERA.dummyVarA = -cost_l.$pDERA.dummyVarA
*/
void include_test_eqFunction_15(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,15};
  jacobian->resultVars[1] /* der(x2.$pDERA.dummyVarA) JACOBIAN_VAR */ = (-jacobian->tmpVars[1] /* cost_l.$pDERA.dummyVarA JACOBIAN_TMP_VAR */);
  TRACE_POP
}

/*
equation index: 16
type: SIMPLE_ASSIGN
$DER.x1.$pDERA.dummyVarA = ((-0.5) * u ^ 2.0 - u) * x1.SeedA
*/
void include_test_eqFunction_16(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,16};
  modelica_real tmp2;
  tmp2 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */);
  jacobian->resultVars[0] /* der(x1.$pDERA.dummyVarA) JACOBIAN_VAR */ = ((-0.5) * ((tmp2 * tmp2)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */);
  TRACE_POP
}

OMC_DISABLE_OPT
int include_test_functionJacA_constantEqns(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_A;
  
  TRACE_POP
  return 0;
}

int include_test_functionJacA_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_A;
  include_test_eqFunction_14(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_15(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_16(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}

int include_test_initialAnalyticJacobianH(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  TRACE_POP
  jacobian->availability = JACOBIAN_NOT_AVAILABLE;
  return 1;
}
int include_test_initialAnalyticJacobianF(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  TRACE_POP
  jacobian->availability = JACOBIAN_NOT_AVAILABLE;
  return 1;
}
int include_test_initialAnalyticJacobianD(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  TRACE_POP
  jacobian->availability = JACOBIAN_NOT_AVAILABLE;
  return 1;
}
OMC_DISABLE_OPT
int include_test_initialAnalyticJacobianC(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  size_t count;

  FILE* pFile = openSparsePatternFile(data, threadData, "include_test_JacC.bin");
  
  initJacobian(jacobian, 3, 4, 5, include_test_functionJacC_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 8, 2);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 8, pFile, FALSE);
  if (count != 8) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 8, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 3);
  /* color 2 with 2 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 2, 2, 3);
  
  omc_fclose(pFile);
  
  TRACE_POP
  return 0;
}
OMC_DISABLE_OPT
int include_test_initialAnalyticJacobianB(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  size_t count;

  FILE* pFile = openSparsePatternFile(data, threadData, "include_test_JacB.bin");
  
  initJacobian(jacobian, 3, 4, 5, include_test_functionJacB_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 8, 2);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 8, pFile, FALSE);
  if (count != 8) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 8, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 3);
  /* color 2 with 2 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 2, 2, 3);
  
  omc_fclose(pFile);
  
  TRACE_POP
  return 0;
}
OMC_DISABLE_OPT
int include_test_initialAnalyticJacobianA(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  size_t count;

  FILE* pFile = openSparsePatternFile(data, threadData, "include_test_JacA.bin");
  
  initJacobian(jacobian, 2, 2, 5, include_test_functionJacA_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(2, 2, 1);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 2+1, pFile, FALSE);
  if (count != 2+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 2+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 2, pFile, FALSE);
  if (count != 2) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 2, count);
  }
  
  /* write color array */
  /* color 1 with 2 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 2, 2);
  
  omc_fclose(pFile);
  
  TRACE_POP
  return 0;
}



