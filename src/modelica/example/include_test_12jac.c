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
equation index: 41
type: SIMPLE_ASSIGN
$cse6 = exp(20.7 - 15599.838969404187 / u)
*/
void include_test_eqFunction_41(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,41};
  jacobian->tmpVars[4] /* $cse6 JACOBIAN_TMP_VAR */ = exp(20.7 - (DIVISION(15599.838969404187,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 42
type: SIMPLE_ASSIGN
$cse7 = exp(18.75 - 14190.821256038647 / u)
*/
void include_test_eqFunction_42(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,42};
  jacobian->tmpVars[3] /* $cse7 JACOBIAN_TMP_VAR */ = exp(18.75 - (DIVISION(14190.821256038647,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 43
type: SIMPLE_ASSIGN
$cse8 = exp(23.67 - 17008.856682769725 / u)
*/
void include_test_eqFunction_43(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,43};
  jacobian->tmpVars[2] /* $cse8 JACOBIAN_TMP_VAR */ = exp(23.67 - (DIVISION(17008.856682769725,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 44
type: SIMPLE_ASSIGN
$cse9 = exp(24.25 - 18820.450885668277 / u)
*/
void include_test_eqFunction_44(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,44};
  jacobian->tmpVars[1] /* $cse9 JACOBIAN_TMP_VAR */ = exp(24.25 - (DIVISION(18820.450885668277,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 45
type: SIMPLE_ASSIGN
$cse10 = exp(8.86 - 10215.37842190016 / u)
*/
void include_test_eqFunction_45(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,45};
  jacobian->tmpVars[0] /* $cse10 JACOBIAN_TMP_VAR */ = exp(8.86 - (DIVISION(10215.37842190016,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 46
type: SIMPLE_ASSIGN
$OMC$objectMayerTerm.$pDERC.dummyVarC = x1.SeedC
*/
void include_test_eqFunction_46(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,46};
  jacobian->resultVars[3] /* $OMC$objectMayerTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 47
type: SIMPLE_ASSIGN
cost_m.$pDERC.dummyVarC = x1.SeedC
*/
void include_test_eqFunction_47(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 6;
  const int equationIndexes[2] = {1,47};
  jacobian->tmpVars[6] /* cost_m.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = jacobian->seedVars[0] /* x1.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 48
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERC.dummyVarC = -x2.SeedC
*/
void include_test_eqFunction_48(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 7;
  const int equationIndexes[2] = {1,48};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = (-jacobian->seedVars[1] /* x2.SeedC SEED_VAR */);
  TRACE_POP
}

/*
equation index: 49
type: SIMPLE_ASSIGN
cost_l.$pDERC.dummyVarC = -x2.SeedC
*/
void include_test_eqFunction_49(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 8;
  const int equationIndexes[2] = {1,49};
  jacobian->tmpVars[5] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (-jacobian->seedVars[1] /* x2.SeedC SEED_VAR */);
  TRACE_POP
}

/*
equation index: 50
type: SIMPLE_ASSIGN
k1.$pDERC.dummyVarC = 10215.37842190016 * $cse10 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_50(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 9;
  const int equationIndexes[2] = {1,50};
  modelica_real tmp0;
  tmp0 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[7] /* k1.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (10215.37842190016) * ((jacobian->tmpVars[0] /* $cse10 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp0 * tmp0),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 51
type: SIMPLE_ASSIGN
k2.$pDERC.dummyVarC = 18820.450885668277 * $cse9 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_51(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 10;
  const int equationIndexes[2] = {1,51};
  modelica_real tmp1;
  tmp1 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[8] /* k2.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (18820.450885668277) * ((jacobian->tmpVars[1] /* $cse9 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp1 * tmp1),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 52
type: SIMPLE_ASSIGN
k3.$pDERC.dummyVarC = 17008.856682769725 * $cse8 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_52(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 11;
  const int equationIndexes[2] = {1,52};
  modelica_real tmp2;
  tmp2 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[9] /* k3.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (17008.856682769725) * ((jacobian->tmpVars[2] /* $cse8 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp2 * tmp2),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 53
type: SIMPLE_ASSIGN
k4.$pDERC.dummyVarC = 14190.821256038647 * $cse7 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_53(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 12;
  const int equationIndexes[2] = {1,53};
  modelica_real tmp3;
  tmp3 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[10] /* k4.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (14190.821256038647) * ((jacobian->tmpVars[3] /* $cse7 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp3 * tmp3),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 54
type: SIMPLE_ASSIGN
k5.$pDERC.dummyVarC = 15599.838969404187 * $cse6 * u.SeedC / u ^ 2.0
*/
void include_test_eqFunction_54(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 13;
  const int equationIndexes[2] = {1,54};
  modelica_real tmp4;
  tmp4 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[11] /* k5.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (15599.838969404187) * ((jacobian->tmpVars[4] /* $cse6 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedC SEED_VAR */,(tmp4 * tmp4),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 55
type: SIMPLE_ASSIGN
$DER.x2.$pDERC.dummyVarC = k1 * x1.SeedC + k1.$pDERC.dummyVarC * x1 + k3 * (x1 * x2.SeedC + x1.SeedC * x2) + k3.$pDERC.dummyVarC * x1 * x2 + (-k2) * x2.SeedC - k2.$pDERC.dummyVarC * x2
*/
void include_test_eqFunction_55(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 14;
  const int equationIndexes[2] = {1,55};
  jacobian->resultVars[1] /* der(x2.$pDERC.dummyVarC) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) + (jacobian->tmpVars[7] /* k1.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */)) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) + (jacobian->tmpVars[9] /* k3.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k2 variable */))) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) - ((jacobian->tmpVars[8] /* k2.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 56
type: SIMPLE_ASSIGN
$DER.x1.$pDERC.dummyVarC = (((-k5) - k3 - k4) * x2 - k1) * x1.SeedC + (((-k5) - k3 - k4) * x2.SeedC + ((-k5.$pDERC.dummyVarC) - k3.$pDERC.dummyVarC - k4.$pDERC.dummyVarC) * x2 - k1.$pDERC.dummyVarC) * x1
*/
void include_test_eqFunction_56(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 15;
  const int equationIndexes[2] = {1,56};
  jacobian->resultVars[0] /* der(x1.$pDERC.dummyVarC) JACOBIAN_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k4 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k4 variable */)) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) + ((-jacobian->tmpVars[11] /* k5.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) - jacobian->tmpVars[9] /* k3.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ - jacobian->tmpVars[10] /* k4.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - jacobian->tmpVars[7] /* k1.$pDERC.dummyVarC JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
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
  include_test_eqFunction_53(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_54(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_55(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_56(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 27
type: SIMPLE_ASSIGN
$cse1 = exp(20.7 - 15599.838969404187 / u)
*/
void include_test_eqFunction_27(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,27};
  jacobian->tmpVars[4] /* $cse1 JACOBIAN_TMP_VAR */ = exp(20.7 - (DIVISION(15599.838969404187,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 28
type: SIMPLE_ASSIGN
$cse2 = exp(18.75 - 14190.821256038647 / u)
*/
void include_test_eqFunction_28(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,28};
  jacobian->tmpVars[3] /* $cse2 JACOBIAN_TMP_VAR */ = exp(18.75 - (DIVISION(14190.821256038647,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 29
type: SIMPLE_ASSIGN
$cse3 = exp(23.67 - 17008.856682769725 / u)
*/
void include_test_eqFunction_29(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,29};
  jacobian->tmpVars[2] /* $cse3 JACOBIAN_TMP_VAR */ = exp(23.67 - (DIVISION(17008.856682769725,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 30
type: SIMPLE_ASSIGN
$cse4 = exp(24.25 - 18820.450885668277 / u)
*/
void include_test_eqFunction_30(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,30};
  jacobian->tmpVars[1] /* $cse4 JACOBIAN_TMP_VAR */ = exp(24.25 - (DIVISION(18820.450885668277,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 31
type: SIMPLE_ASSIGN
$cse5 = exp(8.86 - 10215.37842190016 / u)
*/
void include_test_eqFunction_31(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,31};
  jacobian->tmpVars[0] /* $cse5 JACOBIAN_TMP_VAR */ = exp(8.86 - (DIVISION(10215.37842190016,(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */),"u")));
  TRACE_POP
}

/*
equation index: 32
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERB.dummyVarB = -x2.SeedB
*/
void include_test_eqFunction_32(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,32};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERB.dummyVarB JACOBIAN_VAR */ = (-jacobian->seedVars[1] /* x2.SeedB SEED_VAR */);
  TRACE_POP
}

/*
equation index: 33
type: SIMPLE_ASSIGN
cost_l.$pDERB.dummyVarB = -x2.SeedB
*/
void include_test_eqFunction_33(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 6;
  const int equationIndexes[2] = {1,33};
  jacobian->tmpVars[6] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (-jacobian->seedVars[1] /* x2.SeedB SEED_VAR */);
  TRACE_POP
}

/*
equation index: 34
type: SIMPLE_ASSIGN
k1.$pDERB.dummyVarB = 10215.37842190016 * $cse5 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_34(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 7;
  const int equationIndexes[2] = {1,34};
  modelica_real tmp5;
  tmp5 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[8] /* k1.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (10215.37842190016) * ((jacobian->tmpVars[0] /* $cse5 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp5 * tmp5),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 35
type: SIMPLE_ASSIGN
k2.$pDERB.dummyVarB = 18820.450885668277 * $cse4 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_35(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 8;
  const int equationIndexes[2] = {1,35};
  modelica_real tmp6;
  tmp6 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[9] /* k2.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (18820.450885668277) * ((jacobian->tmpVars[1] /* $cse4 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp6 * tmp6),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 36
type: SIMPLE_ASSIGN
k3.$pDERB.dummyVarB = 17008.856682769725 * $cse3 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_36(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 9;
  const int equationIndexes[2] = {1,36};
  modelica_real tmp7;
  tmp7 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[10] /* k3.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (17008.856682769725) * ((jacobian->tmpVars[2] /* $cse3 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp7 * tmp7),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 37
type: SIMPLE_ASSIGN
k4.$pDERB.dummyVarB = 14190.821256038647 * $cse2 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_37(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 10;
  const int equationIndexes[2] = {1,37};
  modelica_real tmp8;
  tmp8 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[11] /* k4.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (14190.821256038647) * ((jacobian->tmpVars[3] /* $cse2 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp8 * tmp8),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 38
type: SIMPLE_ASSIGN
k5.$pDERB.dummyVarB = 15599.838969404187 * $cse1 * u.SeedB / u ^ 2.0
*/
void include_test_eqFunction_38(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 11;
  const int equationIndexes[2] = {1,38};
  modelica_real tmp9;
  tmp9 = (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[13]] /* u variable */);
  jacobian->tmpVars[12] /* k5.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = (15599.838969404187) * ((jacobian->tmpVars[4] /* $cse1 JACOBIAN_TMP_VAR */) * (DIVISION(jacobian->seedVars[2] /* u.SeedB SEED_VAR */,(tmp9 * tmp9),"u ^ 2.0")));
  TRACE_POP
}

/*
equation index: 39
type: SIMPLE_ASSIGN
$DER.x2.$pDERB.dummyVarB = k1 * x1.SeedB + k1.$pDERB.dummyVarB * x1 + k3 * (x1 * x2.SeedB + x1.SeedB * x2) + k3.$pDERB.dummyVarB * x1 * x2 + (-k2) * x2.SeedB - k2.$pDERB.dummyVarB * x2
*/
void include_test_eqFunction_39(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 12;
  const int equationIndexes[2] = {1,39};
  jacobian->resultVars[1] /* der(x2.$pDERB.dummyVarB) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) + (jacobian->tmpVars[8] /* k1.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */)) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) + (jacobian->tmpVars[10] /* k3.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k2 variable */))) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) - ((jacobian->tmpVars[9] /* k2.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 40
type: SIMPLE_ASSIGN
$DER.x1.$pDERB.dummyVarB = (((-k5) - k3 - k4) * x2 - k1) * x1.SeedB + (((-k5) - k3 - k4) * x2.SeedB + ((-k5.$pDERB.dummyVarB) - k3.$pDERB.dummyVarB - k4.$pDERB.dummyVarB) * x2 - k1.$pDERB.dummyVarB) * x1
*/
void include_test_eqFunction_40(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 13;
  const int equationIndexes[2] = {1,40};
  jacobian->resultVars[0] /* der(x1.$pDERB.dummyVarB) JACOBIAN_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k4 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) + (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k5 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k4 variable */)) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) + ((-jacobian->tmpVars[12] /* k5.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) - jacobian->tmpVars[10] /* k3.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ - jacobian->tmpVars[11] /* k4.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - jacobian->tmpVars[8] /* k1.$pDERB.dummyVarB JACOBIAN_TMP_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
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
  include_test_eqFunction_27(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_28(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_29(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_30(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_31(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_32(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_33(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_34(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_35(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_36(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_37(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_38(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_39(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_40(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 25
type: SIMPLE_ASSIGN
$DER.x1.$pDERA.dummyVarA = (((-k3) - k4 - k5) * x2 - k1) * x1.SeedA + ((-k3) - k4 - k5) * x2.SeedA * x1
*/
void include_test_eqFunction_25(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,25};
  jacobian->resultVars[0] /* der(x1.$pDERA.dummyVarA) JACOBIAN_VAR */ = (((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k4 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k5 variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) + ((-(data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */)) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[11]] /* k4 variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[12]] /* k5 variable */)) * ((jacobian->seedVars[1] /* x2.SeedA SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)));
  TRACE_POP
}

/*
equation index: 26
type: SIMPLE_ASSIGN
$DER.x2.$pDERA.dummyVarA = k1 * x1.SeedA + k3 * (x1 * x2.SeedA + x1.SeedA * x2) - k2 * x2.SeedA
*/
void include_test_eqFunction_26(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,26};
  jacobian->resultVars[1] /* der(x2.$pDERA.dummyVarA) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[8]] /* k1 variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) + ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[10]] /* k3 variable */)) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedA SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* k2 variable */)) * (jacobian->seedVars[1] /* x2.SeedA SEED_VAR */));
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
  include_test_eqFunction_25(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_26(data, threadData, jacobian, parentJacobian);
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
  
  initJacobian(jacobian, 3, 4, 16, include_test_functionJacC_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 9, 3);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 9, pFile, FALSE);
  if (count != 9) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 9, count);
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
  
  initJacobian(jacobian, 3, 3, 16, include_test_functionJacB_column, NULL, NULL);
  jacobian->sparsePattern = allocSparsePattern(3, 7, 3);
  jacobian->availability = JACOBIAN_AVAILABLE;
  
  /* read lead index of compressed sparse column */
  count = omc_fread(jacobian->sparsePattern->leadindex, sizeof(unsigned int), 3+1, pFile, FALSE);
  if (count != 3+1) {
    throwStreamPrint(threadData, "Error while reading lead index list of sparsity pattern. Expected %d, got %zu", 3+1, count);
  }
  
  /* read sparse index */
  count = omc_fread(jacobian->sparsePattern->index, sizeof(unsigned int), 7, pFile, FALSE);
  if (count != 7) {
    throwStreamPrint(threadData, "Error while reading row index list of sparsity pattern. Expected %d, got %zu", 7, count);
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
  
  initJacobian(jacobian, 2, 2, 11, include_test_functionJacA_column, NULL, NULL);
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



