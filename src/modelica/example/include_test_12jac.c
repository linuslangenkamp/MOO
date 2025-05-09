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
equation index: 41
type: SIMPLE_ASSIGN
$DER.x1.$pDERD.dummyVarD = ((-3.0) * u * x2 - u) * x1.SeedD + ((-3.0) * u * x2.SeedD + (-3.0) * u.SeedD * x2 - u.SeedD) * x1
*/
void include_test_eqFunction_41(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,41};
  jacobian->tmpVars[0] /* der(x1.$pDERD.dummyVarD) JACOBIAN_TMP_VAR */ = ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedD SEED_VAR */) + ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[1] /* x2.SeedD SEED_VAR */)) + (-3.0) * ((jacobian->seedVars[2] /* u.SeedD SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - jacobian->seedVars[2] /* u.SeedD SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

/*
equation index: 42
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERD.dummyVarD = x1 * x2.SeedD + x1.SeedD * x2
*/
void include_test_eqFunction_42(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,42};
  jacobian->tmpVars[4] /* FINALCONSTR.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedD SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedD SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 43
type: SIMPLE_ASSIGN
$finalCon$FINALCONSTR.$pDERD.dummyVarD = FINALCONSTR.$pDERD.dummyVarD
*/
void include_test_eqFunction_43(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,43};
  jacobian->resultVars[0] /* $finalCon$FINALCONSTR.$pDERD.dummyVarD JACOBIAN_VAR */ = jacobian->tmpVars[4] /* FINALCONSTR.$pDERD.dummyVarD JACOBIAN_TMP_VAR */;
  TRACE_POP
}

/*
equation index: 44
type: SIMPLE_ASSIGN
$DER.x2.$pDERD.dummyVarD = u * (x1.SeedD + FINALCONSTR.$pDERD.dummyVarD - x2.SeedD) + u.SeedD * (x1 + FINALCONSTR - x2)
*/
void include_test_eqFunction_44(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,44};
  jacobian->tmpVars[1] /* der(x2.$pDERD.dummyVarD) JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedD SEED_VAR */ + jacobian->tmpVars[4] /* FINALCONSTR.$pDERD.dummyVarD JACOBIAN_TMP_VAR */ - jacobian->seedVars[1] /* x2.SeedD SEED_VAR */) + (jacobian->seedVars[2] /* u.SeedD SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) + (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
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
  include_test_eqFunction_41(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_42(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_43(data, threadData, jacobian, parentJacobian);
  include_test_eqFunction_44(data, threadData, jacobian, parentJacobian);
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 33
type: SIMPLE_ASSIGN
$OMC$objectMayerTerm.$pDERC.dummyVarC = -x1.SeedC
*/
void include_test_eqFunction_33(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,33};
  jacobian->resultVars[3] /* $OMC$objectMayerTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = (-jacobian->seedVars[0] /* x1.SeedC SEED_VAR */);
  TRACE_POP
}

/*
equation index: 34
type: SIMPLE_ASSIGN
cost_m.$pDERC.dummyVarC = -x1.SeedC
*/
void include_test_eqFunction_34(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,34};
  jacobian->tmpVars[2] /* cost_m.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = (-jacobian->seedVars[0] /* x1.SeedC SEED_VAR */);
  TRACE_POP
}

/*
equation index: 35
type: SIMPLE_ASSIGN
$DER.x1.$pDERC.dummyVarC = ((-3.0) * u * x2 - u) * x1.SeedC + ((-3.0) * u * x2.SeedC + (-3.0) * u.SeedC * x2 - u.SeedC) * x1
*/
void include_test_eqFunction_35(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,35};
  jacobian->resultVars[0] /* der(x1.$pDERC.dummyVarC) JACOBIAN_VAR */ = ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) + ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */)) + (-3.0) * ((jacobian->seedVars[2] /* u.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - jacobian->seedVars[2] /* u.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

/*
equation index: 36
type: SIMPLE_ASSIGN
$con$CONSTR.$pDERC.dummyVarC = x1.SeedC + u.SeedC
*/
void include_test_eqFunction_36(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,36};
  jacobian->resultVars[4] /* $con$CONSTR.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedC SEED_VAR */ + jacobian->seedVars[2] /* u.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 37
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERC.dummyVarC = x1 * x2.SeedC + x1.SeedC * x2
*/
void include_test_eqFunction_37(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,37};
  jacobian->tmpVars[0] /* FINALCONSTR.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 38
type: SIMPLE_ASSIGN
$DER.x2.$pDERC.dummyVarC = u * (x1.SeedC + FINALCONSTR.$pDERC.dummyVarC - x2.SeedC) + u.SeedC * (x1 + FINALCONSTR - x2)
*/
void include_test_eqFunction_38(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,38};
  jacobian->resultVars[1] /* der(x2.$pDERC.dummyVarC) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedC SEED_VAR */ + jacobian->tmpVars[0] /* FINALCONSTR.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ - jacobian->seedVars[1] /* x2.SeedC SEED_VAR */) + (jacobian->seedVars[2] /* u.SeedC SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) + (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 39
type: SIMPLE_ASSIGN
cost_l.$pDERC.dummyVarC = u.SeedC - x2.SeedC
*/
void include_test_eqFunction_39(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 6;
  const int equationIndexes[2] = {1,39};
  jacobian->tmpVars[1] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */ = jacobian->seedVars[2] /* u.SeedC SEED_VAR */ - jacobian->seedVars[1] /* x2.SeedC SEED_VAR */;
  TRACE_POP
}

/*
equation index: 40
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERC.dummyVarC = cost_l.$pDERC.dummyVarC
*/
void include_test_eqFunction_40(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 7;
  const int equationIndexes[2] = {1,40};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERC.dummyVarC JACOBIAN_VAR */ = jacobian->tmpVars[1] /* cost_l.$pDERC.dummyVarC JACOBIAN_TMP_VAR */;
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
equation index: 27
type: SIMPLE_ASSIGN
$DER.x1.$pDERB.dummyVarB = ((-3.0) * u * x2 - u) * x1.SeedB + ((-3.0) * u * x2.SeedB + (-3.0) * u.SeedB * x2 - u.SeedB) * x1
*/
void include_test_eqFunction_27(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,27};
  jacobian->resultVars[0] /* der(x1.$pDERB.dummyVarB) JACOBIAN_VAR */ = ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) + ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */)) + (-3.0) * ((jacobian->seedVars[2] /* u.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - jacobian->seedVars[2] /* u.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */));
  TRACE_POP
}

/*
equation index: 28
type: SIMPLE_ASSIGN
$con$CONSTR.$pDERB.dummyVarB = x1.SeedB + u.SeedB
*/
void include_test_eqFunction_28(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,28};
  jacobian->resultVars[3] /* $con$CONSTR.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->seedVars[0] /* x1.SeedB SEED_VAR */ + jacobian->seedVars[2] /* u.SeedB SEED_VAR */;
  TRACE_POP
}

/*
equation index: 29
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERB.dummyVarB = x1 * x2.SeedB + x1.SeedB * x2
*/
void include_test_eqFunction_29(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,29};
  jacobian->tmpVars[1] /* FINALCONSTR.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 30
type: SIMPLE_ASSIGN
$DER.x2.$pDERB.dummyVarB = u * (x1.SeedB + FINALCONSTR.$pDERB.dummyVarB - x2.SeedB) + u.SeedB * (x1 + FINALCONSTR - x2)
*/
void include_test_eqFunction_30(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 3;
  const int equationIndexes[2] = {1,30};
  jacobian->resultVars[1] /* der(x2.$pDERB.dummyVarB) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedB SEED_VAR */ + jacobian->tmpVars[1] /* FINALCONSTR.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ - jacobian->seedVars[1] /* x2.SeedB SEED_VAR */) + (jacobian->seedVars[2] /* u.SeedB SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */) + (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 31
type: SIMPLE_ASSIGN
cost_l.$pDERB.dummyVarB = u.SeedB - x2.SeedB
*/
void include_test_eqFunction_31(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 4;
  const int equationIndexes[2] = {1,31};
  jacobian->tmpVars[2] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */ = jacobian->seedVars[2] /* u.SeedB SEED_VAR */ - jacobian->seedVars[1] /* x2.SeedB SEED_VAR */;
  TRACE_POP
}

/*
equation index: 32
type: SIMPLE_ASSIGN
$OMC$objectLagrangeTerm.$pDERB.dummyVarB = cost_l.$pDERB.dummyVarB
*/
void include_test_eqFunction_32(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 5;
  const int equationIndexes[2] = {1,32};
  jacobian->resultVars[2] /* $OMC$objectLagrangeTerm.$pDERB.dummyVarB JACOBIAN_VAR */ = jacobian->tmpVars[2] /* cost_l.$pDERB.dummyVarB JACOBIAN_TMP_VAR */;
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
  TRACE_POP
  return 0;
}
/* constant equations */
/* dynamic equations */

/*
equation index: 24
type: SIMPLE_ASSIGN
$DER.x1.$pDERA.dummyVarA = ((-3.0) * u * x2 - u) * x1.SeedA + (-3.0) * u * x2.SeedA * x1
*/
void include_test_eqFunction_24(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 0;
  const int equationIndexes[2] = {1,24};
  jacobian->resultVars[0] /* der(x1.$pDERA.dummyVarA) JACOBIAN_VAR */ = ((-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */))) - (data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) + (-3.0) * (((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * ((jacobian->seedVars[1] /* x2.SeedA SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */))));
  TRACE_POP
}

/*
equation index: 25
type: SIMPLE_ASSIGN
FINALCONSTR.$pDERA.dummyVarA = x1 * x2.SeedA + x1.SeedA * x2
*/
void include_test_eqFunction_25(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 1;
  const int equationIndexes[2] = {1,25};
  jacobian->tmpVars[2] /* FINALCONSTR.$pDERA.dummyVarA JACOBIAN_TMP_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[0]] /* x1 STATE(1) */)) * (jacobian->seedVars[1] /* x2.SeedA SEED_VAR */) + (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */) * ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[1]] /* x2 STATE(1) */));
  TRACE_POP
}

/*
equation index: 26
type: SIMPLE_ASSIGN
$DER.x2.$pDERA.dummyVarA = u * (x1.SeedA + FINALCONSTR.$pDERA.dummyVarA - x2.SeedA)
*/
void include_test_eqFunction_26(DATA *data, threadData_t *threadData, JACOBIAN *jacobian, JACOBIAN *parentJacobian)
{
  TRACE_PUSH
  const int baseClockIndex = 0;
  const int subClockIndex = 2;
  const int equationIndexes[2] = {1,26};
  jacobian->resultVars[1] /* der(x2.$pDERA.dummyVarA) JACOBIAN_VAR */ = ((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[9]] /* u variable */)) * (jacobian->seedVars[0] /* x1.SeedA SEED_VAR */ + jacobian->tmpVars[2] /* FINALCONSTR.$pDERA.dummyVarA JACOBIAN_TMP_VAR */ - jacobian->seedVars[1] /* x2.SeedA SEED_VAR */);
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
  include_test_eqFunction_24(data, threadData, jacobian, parentJacobian);
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
OMC_DISABLE_OPT
int include_test_initialAnalyticJacobianD(DATA* data, threadData_t *threadData, JACOBIAN *jacobian)
{
  TRACE_PUSH
  size_t count;

  FILE* pFile = openSparsePatternFile(data, threadData, "include_test_JacD.bin");
  
  initJacobian(jacobian, 3, 1, 9, include_test_functionJacD_column, NULL, NULL);
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
  
  initJacobian(jacobian, 3, 5, 9, include_test_functionJacC_column, NULL, NULL);
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
  
  initJacobian(jacobian, 3, 4, 9, include_test_functionJacB_column, NULL, NULL);
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
  
  initJacobian(jacobian, 2, 2, 9, include_test_functionJacA_column, NULL, NULL);
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



