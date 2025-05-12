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
/* constant equations */
/* dynamic equations */

/*
equation index: 71
type: SIMPLE_ASSIGN
$cse11 = exp(8.86 - 10215.37842190016 / u)
*/
void include_test_eqFunction_71(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,71};
  jacobian->tmpVars[4] /* $cse11 JACOBIAN_TMP_VAR */ = exp(8.86 - (DIVISION(10215.37842190016,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 72
type: SIMPLE_ASSIGN
$cse12 = exp(24.25 - 18820.450885668277 / u)
*/
void include_test_eqFunction_72(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,72};
  jacobian->tmpVars[3] /* $cse12 JACOBIAN_TMP_VAR */ = exp(24.25 - (DIVISION(18820.450885668277,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 73
type: SIMPLE_ASSIGN
$cse13 = exp(23.67 - 17008.856682769725 / u)
*/
void include_test_eqFunction_73(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,73};
  jacobian->tmpVars[2] /* $cse13 JACOBIAN_TMP_VAR */ = exp(23.67 - (DIVISION(17008.856682769725,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 74
type: SIMPLE_ASSIGN
$cse14 = exp(18.75 - 14190.821256038647 / u)
*/
void include_test_eqFunction_74(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,74};
  jacobian->tmpVars[1] /* $cse14 JACOBIAN_TMP_VAR */ = exp(18.75 - (DIVISION(14190.821256038647,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 75
type: SIMPLE_ASSIGN
$cse15 = exp(20.7 - 15599.838969404187 / u)
*/
void include_test_eqFunction_75(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,75};
  jacobian->tmpVars[0] /* $cse15 JACOBIAN_TMP_VAR */ = exp(20.7 - (DIVISION(15599.838969404187,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 76
type: SIMPLE_ASSIGN
k5.$pDERD.dummyVarD = 15599.838969404187 * $cse15 * u.SeedD / u ^ 2.0
*/
void include_test_eqFunction_76(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,76};
  modelica_real tmp0;
  tmp0 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[16] /* k5.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = (15599.838969404187) * ((jacobian->tmpVars[0] /* $cse15 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedD SEED_VAR */,(tmp0 * tmp0),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 77
type: SIMPLE_ASSIGN
k4.$pDERD.dummyVarD = 14190.821256038647 * $cse14 * u.SeedD / u ^ 2.0
*/
void include_test_eqFunction_77(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 6;
  const int equationIndexes[2] = {1,77};
  modelica_real tmp1;
  tmp1 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[15] /* k4.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = (14190.821256038647) * ((jacobian->tmpVars[1] /* $cse14 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedD SEED_VAR */,(tmp1 * tmp1),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 78
type: SIMPLE_ASSIGN
k3.$pDERD.dummyVarD = 17008.856682769725 * $cse13 * u.SeedD / u ^ 2.0
*/
void include_test_eqFunction_78(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 7;
  const int equationIndexes[2] = {1,78};
  modelica_real tmp2;
  tmp2 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[14] /* k3.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = (17008.856682769725) * ((jacobian->tmpVars[2] /* $cse13 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedD SEED_VAR */,(tmp2 * tmp2),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 79
type: SIMPLE_ASSIGN
k2.$pDERD.dummyVarD = 18820.450885668277 * $cse12 * u.SeedD / u ^ 2.0
*/
void include_test_eqFunction_79(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 8;
  const int equationIndexes[2] = {1,79};
  modelica_real tmp3;
  tmp3 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[13] /* k2.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = (18820.450885668277) * ((jacobian->tmpVars[3] /* $cse12 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedD SEED_VAR */,(tmp3 * tmp3),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 80
type: SIMPLE_ASSIGN
k1.$pDERD.dummyVarD = 10215.37842190016 * $cse11 * u.SeedD / u ^ 2.0
*/
void include_test_eqFunction_80(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 9;
  const int equationIndexes[2] = {1,80};
  modelica_real tmp4;
  tmp4 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[12] /* k1.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = (10215.37842190016) * ((jacobian->tmpVars[4] /* $cse11 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedD SEED_VAR */,(tmp4 * tmp4),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 81
type: SIMPLE_ASSIGN
$DER.x1.$pDERD.dummyVarD = (((-k5) - k3 - k4) * x2 - k1) * x1.SeedD + (((-k5) - k3 - k4) * x2.SeedD + ((-k5.$pDERD.dummyVarD) - k3.$pDERD.dummyVarD - k4.$pDERD.dummyVarD) * x2 - k1.$pDERD.dummyVarD) * x1
*/
void include_test_eqFunction_81(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 10;
  const int equationIndexes[2] = {1,81};
  jacobian->tmpVars[5] /* der(x1.$pDERD.dummyVarD) JACOBIAN_TMP_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedD SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) * (jacobian->seedVars[1] /* x2.SeedD SEED_VAR */) + ((-jacobian->tmpVars[16] /* k5.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) - jacobian->tmpVars[14] /* k3.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ - jacobian->tmpVars[15] /* k4.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - jacobian->tmpVars[12] /* k1.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

/*
equation index: 82
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERD.dummyVarD = x1 * x2.SeedD + x1.SeedD * x2
*/
void include_test_eqFunction_82(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 11;
  const int equationIndexes[2] = {1,82};
  jacobian->tmpVars[9] /* FINALCONSTR.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedD SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedD SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 83
type: SIMPLE_ASSIGN
$finalCon$FINALCONSTR.$pDERD.dummyVarD = FINALCONSTR.$pDERD.dummyVarD
*/
void include_test_eqFunction_83(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 12;
  const int equationIndexes[2] = {1,83};
  jacobian->resultVars[0] /* $finalCon$FINALCONSTR.$pDERD.dummyVarD JACOBIAN_VAR */ = jacobian->tmpVars[9] /* FINALCONSTR.$pDERD.dummyVarD JACOBIAN_TMP_VAR */;
  TRACE_POP
}

/*
equation index: 84
type: SIMPLE_ASSIGN
$DER.x2.$pDERD.dummyVarD = k1 * x1.SeedD + k1.$pDERD.dummyVarD * x1 + k3 * FINALCONSTR.$pDERD.dummyVarD + k3.$pDERD.dummyVarD * FINALCONSTR + (-k2) * x2.SeedD - k2.$pDERD.dummyVarD * x2
*/
void include_test_eqFunction_84(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 13;
  const int equationIndexes[2] = {1,84};
  jacobian->tmpVars[6] /* der(x2.$pDERD.dummyVarD) JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedD SEED_VAR */) + (jacobian->tmpVars[12] /* k1.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) * (jacobian->tmpVars[9] /* FINALCONSTR.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) + (jacobian->tmpVars[14] /* k3.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */)) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k2 variable */))) * (jacobian->seedVars[1] /* x2.SeedD SEED_VAR */) - ((jacobian->tmpVars[13] /* k2.$pDERD.dummyVarD JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)));
  TRACE_POP
}

OMC_DISABLE_OPT
int include_test_functionJacD_constantEqns(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_D;
  
  TRACE_POP
  return 0;
}

int include_test_functionJacD_column(DATA* data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH

  int index = include_test_INDEX_JAC_D;
  include_test_eqFunction_71(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_72(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_73(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_74(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_75(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_76(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_77(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_78(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_79(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_80(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_81(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_82(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_83(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_84(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 53
type: SIMPLE_ASSIGN
$cse6 = exp(8.86 - 10215.37842190016 / u)
*/
void include_test_eqFunction_53(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,53};
  jacobian->tmpVars[4] /* $cse6 JACOBIAN_TMP_VAR */ = exp(8.86 - (DIVISION(10215.37842190016,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 54
type: SIMPLE_ASSIGN
$cse7 = exp(24.25 - 18820.450885668277 / u)
*/
void include_test_eqFunction_54(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,54};
  jacobian->tmpVars[3] /* $cse7 JACOBIAN_TMP_VAR */ = exp(24.25 - (DIVISION(18820.450885668277,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 55
type: SIMPLE_ASSIGN
$cse8 = exp(23.67 - 17008.856682769725 / u)
*/
void include_test_eqFunction_55(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,55};
  jacobian->tmpVars[2] /* $cse8 JACOBIAN_TMP_VAR */ = exp(23.67 - (DIVISION(17008.856682769725,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 56
type: SIMPLE_ASSIGN
$cse9 = exp(18.75 - 14190.821256038647 / u)
*/
void include_test_eqFunction_56(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,56};
  jacobian->tmpVars[1] /* $cse9 JACOBIAN_TMP_VAR */ = exp(18.75 - (DIVISION(14190.821256038647,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 57
type: SIMPLE_ASSIGN
$cse10 = exp(20.7 - 15599.838969404187 / u)
*/
void include_test_eqFunction_57(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,57};
  jacobian->tmpVars[0] /* $cse10 JACOBIAN_TMP_VAR */ = exp(20.7 - (DIVISION(15599.838969404187,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 58
type: SIMPLE_ASSIGN
$OMC$objectMayerTerm.$pDERC.dummyVarC = -x1.SeedC
*/
void include_test_eqFunction_58(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,58};
  jacobian->resultVars[3] /* $OMC$objectMayerTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = (-jacobian->seedVars[0] /* x1.SeedC SEED_VAR */);
  TRACE_POP
}

/*
equation index: 59
type: SIMPLE_ASSIGN
cost_m.$pDERC.dummyVarC = -x1.SeedC
*/
void include_test_eqFunction_59(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 6;
  const int equationIndexes[2] = {1,59};
  jacobian->tmpVars[7] /* cost_m.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (-jacobian->seedVars[0] /* x1.SeedC SEED_VAR */);
  TRACE_POP
}

/*
equation index: 60
type: SIMPLE_ASSIGN
k5.$pDERC.dummyVarC = 15599.838969404187 * $cse10 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_60(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 7;
  const int equationIndexes[2] = {1,60};
  modelica_real tmp5;
  tmp5 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[12] /* k5.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (15599.838969404187) * ((jacobian->tmpVars[0] /* $cse10 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp5 * tmp5),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 61
type: SIMPLE_ASSIGN
k4.$pDERC.dummyVarC = 14190.821256038647 * $cse9 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_61(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 8;
  const int equationIndexes[2] = {1,61};
  modelica_real tmp6;
  tmp6 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[11] /* k4.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (14190.821256038647) * ((jacobian->tmpVars[1] /* $cse9 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp6 * tmp6),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 62
type: SIMPLE_ASSIGN
k3.$pDERC.dummyVarC = 17008.856682769725 * $cse8 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_62(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 9;
  const int equationIndexes[2] = {1,62};
  modelica_real tmp7;
  tmp7 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[10] /* k3.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (17008.856682769725) * ((jacobian->tmpVars[2] /* $cse8 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp7 * tmp7),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 63
type: SIMPLE_ASSIGN
k2.$pDERC.dummyVarC = 18820.450885668277 * $cse7 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_63(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 10;
  const int equationIndexes[2] = {1,63};
  modelica_real tmp8;
  tmp8 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[9] /* k2.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (18820.450885668277) * ((jacobian->tmpVars[3] /* $cse7 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp8 * tmp8),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 64
type: SIMPLE_ASSIGN
k1.$pDERC.dummyVarC = 10215.37842190016 * $cse6 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_64(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 11;
  const int equationIndexes[2] = {1,64};
  modelica_real tmp9;
  tmp9 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[8] /* k1.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (10215.37842190016) * ((jacobian->tmpVars[4] /* $cse6 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp9 * tmp9),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 65
type: SIMPLE_ASSIGN
$DER.x1.$pDERC.dummyVarC = (((-k5) - k3 - k4) * x2 - k1) * x1.SeedC + (((-k5) - k3 - k4) * x2.SeedC + ((-k5.$pDERC.dummyVarC) - k3.$pDERC.dummyVarC - k4.$pDERC.dummyVarC) * x2 - k1.$pDERC.dummyVarC) * x1
*/
void include_test_eqFunction_65(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 12;
  const int equationIndexes[2] = {1,65};
  jacobian->resultVars[0] /* der(x1.$pDERC.dummyVarC) JACOBIAN_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) + ((-jacobian->tmpVars[12] /* k5.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) - jacobian->tmpVars[10] /* k3.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ - jacobian->tmpVars[11] /* k4.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - jacobian->tmpVars[8] /* k1.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

/*
equation index: 66
type: SIMPLE_ASSIGN
$con$CONSTR.$pDERC.dummyVarC = x1.SeedC + u.SeedC
*/
void include_test_eqFunction_66(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 13;
  const int equationIndexes[2] = {1,66};
  jacobian->resultVars[4] /* $con$CONSTR.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedC SEED_VAR */ + jacobian->seedVars[2] /* u.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 67
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERC.dummyVarC = x1 * x2.SeedC + x1.SeedC * x2
*/
void include_test_eqFunction_67(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 14;
  const int equationIndexes[2] = {1,67};
  jacobian->tmpVars[5] /* FINALCONSTR.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 68
type: SIMPLE_ASSIGN
$DER.x2.$pDERC.dummyVarC = k1 * x1.SeedC + k1.$pDERC.dummyVarC * x1 + k3 * FINALCONSTR.$pDERC.dummyVarC + k3.$pDERC.dummyVarC * FINALCONSTR + (-k2) * x2.SeedC - k2.$pDERC.dummyVarC * x2
*/
void include_test_eqFunction_68(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 15;
  const int equationIndexes[2] = {1,68};
  jacobian->resultVars[1] /* der(x2.$pDERC.dummyVarC) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) + (jacobian->tmpVars[8] /* k1.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) * (jacobian->tmpVars[5] /* FINALCONSTR.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) + (jacobian->tmpVars[10] /* k3.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */)) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k2 variable */))) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) - ((jacobian->tmpVars[9] /* k2.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 69
type: SIMPLE_ASSIGN
cost_l.$pDERC.dummyVarC = u.SeedC - x2.SeedC
*/
void include_test_eqFunction_69(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 16;
  const int equationIndexes[2] = {1,69};
  jacobian->tmpVars[6] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = jacobian->seedVars[2] /* u.SeedC SEED_VAR */ - jacobian->seedVars[1] /* x2.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 70
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERC.dummyVarC = cost_l.$pDERC.dummyVarC
*/
void include_test_eqFunction_70(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 17;
  const int equationIndexes[2] = {1,70};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->tmpVars[6] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */;
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
  include_test_eqFunction_53(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_54(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_55(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_56(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_57(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_58(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_59(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_60(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_61(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_62(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_63(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_64(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_65(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_66(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_67(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_68(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_69(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_70(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 37
type: SIMPLE_ASSIGN
$cse1 = exp(8.86 - 10215.37842190016 / u)
*/
void include_test_eqFunction_37(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,37};
  jacobian->tmpVars[4] /* $cse1 JACOBIAN_TMP_VAR */ = exp(8.86 - (DIVISION(10215.37842190016,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 38
type: SIMPLE_ASSIGN
$cse2 = exp(24.25 - 18820.450885668277 / u)
*/
void include_test_eqFunction_38(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,38};
  jacobian->tmpVars[3] /* $cse2 JACOBIAN_TMP_VAR */ = exp(24.25 - (DIVISION(18820.450885668277,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 39
type: SIMPLE_ASSIGN
$cse3 = exp(23.67 - 17008.856682769725 / u)
*/
void include_test_eqFunction_39(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,39};
  jacobian->tmpVars[2] /* $cse3 JACOBIAN_TMP_VAR */ = exp(23.67 - (DIVISION(17008.856682769725,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 40
type: SIMPLE_ASSIGN
$cse4 = exp(18.75 - 14190.821256038647 / u)
*/
void include_test_eqFunction_40(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,40};
  jacobian->tmpVars[1] /* $cse4 JACOBIAN_TMP_VAR */ = exp(18.75 - (DIVISION(14190.821256038647,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 41
type: SIMPLE_ASSIGN
$cse5 = exp(20.7 - 15599.838969404187 / u)
*/
void include_test_eqFunction_41(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,41};
  jacobian->tmpVars[0] /* $cse5 JACOBIAN_TMP_VAR */ = exp(20.7 - (DIVISION(15599.838969404187,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 42
type: SIMPLE_ASSIGN
k5.$pDERB.dummyVarB = 15599.838969404187 * $cse5 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_42(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,42};
  modelica_real tmp10;
  tmp10 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[13] /* k5.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (15599.838969404187) * ((jacobian->tmpVars[0] /* $cse5 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp10 * tmp10),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 43
type: SIMPLE_ASSIGN
k4.$pDERB.dummyVarB = 14190.821256038647 * $cse4 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_43(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 6;
  const int equationIndexes[2] = {1,43};
  modelica_real tmp11;
  tmp11 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[12] /* k4.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (14190.821256038647) * ((jacobian->tmpVars[1] /* $cse4 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp11 * tmp11),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 44
type: SIMPLE_ASSIGN
k3.$pDERB.dummyVarB = 17008.856682769725 * $cse3 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_44(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 7;
  const int equationIndexes[2] = {1,44};
  modelica_real tmp12;
  tmp12 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[11] /* k3.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (17008.856682769725) * ((jacobian->tmpVars[2] /* $cse3 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp12 * tmp12),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 45
type: SIMPLE_ASSIGN
k2.$pDERB.dummyVarB = 18820.450885668277 * $cse2 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_45(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 8;
  const int equationIndexes[2] = {1,45};
  modelica_real tmp13;
  tmp13 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[10] /* k2.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (18820.450885668277) * ((jacobian->tmpVars[3] /* $cse2 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp13 * tmp13),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 46
type: SIMPLE_ASSIGN
k1.$pDERB.dummyVarB = 10215.37842190016 * $cse1 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_46(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 9;
  const int equationIndexes[2] = {1,46};
  modelica_real tmp14;
  tmp14 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* u variable */);
  jacobian->tmpVars[9] /* k1.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (10215.37842190016) * ((jacobian->tmpVars[4] /* $cse1 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp14 * tmp14),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 47
type: SIMPLE_ASSIGN
$DER.x1.$pDERB.dummyVarB = (((-k5) - k3 - k4) * x2 - k1) * x1.SeedB + (((-k5) - k3 - k4) * x2.SeedB + ((-k5.$pDERB.dummyVarB) - k3.$pDERB.dummyVarB - k4.$pDERB.dummyVarB) * x2 - k1.$pDERB.dummyVarB) * x1
*/
void include_test_eqFunction_47(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 10;
  const int equationIndexes[2] = {1,47};
  jacobian->resultVars[0] /* der(x1.$pDERB.dummyVarB) JACOBIAN_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */)) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) + ((-jacobian->tmpVars[13] /* k5.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) - jacobian->tmpVars[11] /* k3.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ - jacobian->tmpVars[12] /* k4.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - jacobian->tmpVars[9] /* k1.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

/*
equation index: 48
type: SIMPLE_ASSIGN
$con$CONSTR.$pDERB.dummyVarB = x1.SeedB + u.SeedB
*/
void include_test_eqFunction_48(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 11;
  const int equationIndexes[2] = {1,48};
  jacobian->resultVars[3] /* $con$CONSTR.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedB SEED_VAR */ + jacobian->seedVars[2] /* u.SeedB SEED_VAR */;
  TRACE_POP
}

/*
equation index: 49
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERB.dummyVarB = x1 * x2.SeedB + x1.SeedB * x2
*/
void include_test_eqFunction_49(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 12;
  const int equationIndexes[2] = {1,49};
  jacobian->tmpVars[6] /* FINALCONSTR.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 50
type: SIMPLE_ASSIGN
$DER.x2.$pDERB.dummyVarB = k1 * x1.SeedB + k1.$pDERB.dummyVarB * x1 + k3 * FINALCONSTR.$pDERB.dummyVarB + k3.$pDERB.dummyVarB * FINALCONSTR + (-k2) * x2.SeedB - k2.$pDERB.dummyVarB * x2
*/
void include_test_eqFunction_50(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 13;
  const int equationIndexes[2] = {1,50};
  jacobian->resultVars[1] /* der(x2.$pDERB.dummyVarB) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) + (jacobian->tmpVars[9] /* k1.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) * (jacobian->tmpVars[6] /* FINALCONSTR.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) + (jacobian->tmpVars[11] /* k3.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */)) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k2 variable */))) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) - ((jacobian->tmpVars[10] /* k2.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 51
type: SIMPLE_ASSIGN
cost_l.$pDERB.dummyVarB = u.SeedB - x2.SeedB
*/
void include_test_eqFunction_51(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 14;
  const int equationIndexes[2] = {1,51};
  jacobian->tmpVars[7] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = jacobian->seedVars[2] /* u.SeedB SEED_VAR */ - jacobian->seedVars[1] /* x2.SeedB SEED_VAR */;
  TRACE_POP
}

/*
equation index: 52
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERB.dummyVarB = cost_l.$pDERB.dummyVarB
*/
void include_test_eqFunction_52(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 15;
  const int equationIndexes[2] = {1,52};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->tmpVars[7] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */;
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
  include_test_eqFunction_37(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_38(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_39(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_40(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_41(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_42(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_43(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_44(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_45(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_46(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_47(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_48(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_49(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_50(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_51(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_52(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 34
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERA.dummyVarA = x1 * x2.SeedA + x1.SeedA * x2
*/
void include_test_eqFunction_34(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,34};
  jacobian->tmpVars[2] /* FINALCONSTR.$pDERA.dummyVarA JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedA SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 35
type: SIMPLE_ASSIGN
$DER.x2.$pDERA.dummyVarA = k1 * x1.SeedA + k3 * FINALCONSTR.$pDERA.dummyVarA - k2 * x2.SeedA
*/
void include_test_eqFunction_35(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,35};
  jacobian->resultVars[1] /* der(x2.$pDERA.dummyVarA) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) * (jacobian->tmpVars[2] /* FINALCONSTR.$pDERA.dummyVarA JACOBIAN_TMP_VAR */) - (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k2 variable */)) * (jacobian->seedVars[1] /* x2.SeedA SEED_VAR */));
  TRACE_POP
}

/*
equation index: 36
type: SIMPLE_ASSIGN
$DER.x1.$pDERA.dummyVarA = (((-k3) - k4 - k5) * x2 - k1) * x1.SeedA + ((-k3) - k4 - k5) * x2.SeedA * x1
*/
void include_test_eqFunction_36(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,36};
  jacobian->resultVars[0] /* der(x1.$pDERA.dummyVarA) JACOBIAN_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k3 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k4 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* k5 variable */)) * ((jacobian->seedVars[1] /* x2.SeedA SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)));
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
  include_test_eqFunction_34(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_35(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_36(data, threadData, jacobian, parentJacobian);
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
OMC_DISABLE_OPT
int include_test_initialAnalyticJacobianD(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  size_t count;

  FILE* pFile = openSparsePatternFile(data, threadData, "include_test_JacD.bin");
  
  initJacobian(jacobian, 3, 1, 19, include_test_functionJacD_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 2, 2);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 2, pFile, FALSE);
  if (count != 2) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 2, count);
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
int include_test_initialAnalyticJacobianC(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  size_t count;

  FILE* pFile = openSparsePatternFile(data, threadData, "include_test_JacC.bin");
  
  initJacobian(jacobian, 3, 5, 19, include_test_functionJacC_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 11, 3);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 11, pFile, FALSE);
  if (count != 11) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 11, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 3);
  /* color 2 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 2, 1, 3);
  /* color 3 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 3, 1, 3);
  
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
  
  initJacobian(jacobian, 3, 4, 19, include_test_functionJacB_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 10, 3);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 10, pFile, FALSE);
  if (count != 10) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 10, count);
  }
  
  /* write color array */
  /* color 1 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 1, 1, 3);
  /* color 2 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 2, 1, 3);
  /* color 3 with 1 columns */
  readSparsePatternColor(threadData, pFile, jacobian->sparsePattern->colorCols, 3, 1, 3);
  
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
  
  initJacobian(jacobian, 2, 2, 14, include_test_functionJacA_column, NULL, NULL);
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



