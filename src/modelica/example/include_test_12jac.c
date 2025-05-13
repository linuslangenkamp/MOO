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
equation index: 12
type: SIMPLE_ASSIGN
$DER.x.$pDERC.dummyVarC = x.SeedC + u.SeedC
*/
void include_test_eqFunction_12(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,12};
  jacobian->resultVars[0] /* der(x.$pDERC.dummyVarC) JACOBIAN_VAR */ = jacobian->seedVars[0] /* x.SeedC SEED_VAR */ + jacobian->seedVars[1] /* u.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 13
type: SIMPLE_ASSIGN
cost_l.$pDERC.dummyVarC = 2.0 * (x * x.SeedC + u * u.SeedC)
*/
void include_test_eqFunction_13(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,13};
  jacobian->tmpVars[0] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (2.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x STATE(1) */)) * (jacobian->seedVars[0] /* x.SeedC SEED_VAR */) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[4]] /* u variable */)) * (jacobian->seedVars[1] /* u.SeedC SEED_VAR */));
  TRACE_POP
}

/*
equation index: 14
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERC.dummyVarC = cost_l.$pDERC.dummyVarC
*/
void include_test_eqFunction_14(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,14};
  jacobian->resultVars[1] /* $OMC$objectLagrangeTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->tmpVars[0] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */;
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
  include_test_eqFunction_12(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_13(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_14(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 9
type: SIMPLE_ASSIGN
$DER.x.$pDERB.dummyVarB = x.SeedB + u.SeedB
*/
void include_test_eqFunction_9(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,9};
  jacobian->resultVars[0] /* der(x.$pDERB.dummyVarB) JACOBIAN_VAR */ = jacobian->seedVars[0] /* x.SeedB SEED_VAR */ + jacobian->seedVars[1] /* u.SeedB SEED_VAR */;
  TRACE_POP
}

/*
equation index: 10
type: SIMPLE_ASSIGN
cost_l.$pDERB.dummyVarB = 2.0 * (x * x.SeedB + u * u.SeedB)
*/
void include_test_eqFunction_10(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,10};
  jacobian->tmpVars[0] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (2.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x STATE(1) */)) * (jacobian->seedVars[0] /* x.SeedB SEED_VAR */) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[4]] /* u variable */)) * (jacobian->seedVars[1] /* u.SeedB SEED_VAR */));
  TRACE_POP
}

/*
equation index: 11
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERB.dummyVarB = cost_l.$pDERB.dummyVarB
*/
void include_test_eqFunction_11(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,11};
  jacobian->resultVars[1] /* $OMC$objectLagrangeTerm.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->tmpVars[0] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */;
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
  include_test_eqFunction_9(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_10(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_11(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 8
type: SIMPLE_ASSIGN
$DER.x.$pDERA.dummyVarA = x.SeedA
*/
void include_test_eqFunction_8(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,8};
  jacobian->resultVars[0] /* der(x.$pDERA.dummyVarA) JACOBIAN_VAR */ = jacobian->seedVars[0] /* x.SeedA SEED_VAR */;
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
  include_test_eqFunction_8(data, threadData, jacobian, parentJacobian);
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
  
  initJacobian(jacobian, 2, 2, 3, include_test_functionJacC_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(2, 4, 2);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 2+1, pFile, FALSE);
  if (count != 2+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 2+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 4, pFile, FALSE);
  if (count != 4) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 4, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 2);
  /* color 2 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 2, 1, 2);
  
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
  
  initJacobian(jacobian, 2, 2, 3, include_test_functionJacB_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(2, 4, 2);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 2+1, pFile, FALSE);
  if (count != 2+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 2+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 4, pFile, FALSE);
  if (count != 4) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 4, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 2);
  /* color 2 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 2, 1, 2);
  
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
  
  initJacobian(jacobian, 1, 1, 3, include_test_functionJacA_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(1, 1, 1);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 1+1, pFile, FALSE);
  if (count != 1+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 1+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 1, pFile, FALSE);
  if (count != 1) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 1, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 1);
  
  omc_fclose(pFile);
  
  TRACE_POP
  return 0;
}



