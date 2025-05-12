/* Asserts */
#include "include_test_model.h"
#if defined(__cplusplus)
extern "C" {
#endif


/*
equation index: 31
type: ALGORITHM

  assert($finalCon$FINALCONSTR >= -1.25 and $finalCon$FINALCONSTR <= 1.25, "Variable violating min/max constraint: -1.25 <= $finalCon$FINALCONSTR <= 1.25, has value: " + String($finalCon$FINALCONSTR, "g"));
*/
void include_test_eqFunction_31(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,31};
  modelica_boolean tmp0;
  modelica_boolean tmp1;
  static const MMC_DEFSTRINGLIT(tmp2,90,"Variable violating min/max constraint: -1.25 <= $finalCon$FINALCONSTR <= 1.25, has value: ");
  modelica_string tmp3;
  modelica_metatype tmpMeta4;
  static int tmp5 = 0;
  if(!tmp5)
  {
    tmp0 = GreaterEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[16]] /* $finalCon$FINALCONSTR OPT_FCONSTR */),-1.25);
    tmp1 = LessEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[16]] /* $finalCon$FINALCONSTR OPT_FCONSTR */),1.25);
    if(!(tmp0 && tmp1))
    {
      tmp3 = modelica_real_to_modelica_string_format((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[16]] /* $finalCon$FINALCONSTR OPT_FCONSTR */), (modelica_string) mmc_strings_len1[103]);
      tmpMeta4 = stringAppend(MMC_REFSTRINGLIT(tmp2),tmp3);
      {
        const char* assert_cond = "($finalCon$FINALCONSTR >= -1.25 and $finalCon$FINALCONSTR <= 1.25)";
        if (data->simulationInfo->noThrowAsserts) {
          FILE_INFO info = {"",0,0,0,0,0};
          infoStreamPrintWithEquationIndexes(OMC_LOG_ASSERT, info, 0, equationIndexes, "The following assertion has been violated %sat time %f\n(%s) --> \"%s\"", initial() ? "during initialization " : "", data->localData[0]->timeValue, assert_cond, MMC_STRINGDATA(tmpMeta4));
        } else {
          FILE_INFO info = {"",0,0,0,0,0};
          omc_assert_warning_withEquationIndexes(info, equationIndexes, "The following assertion has been violated %sat time %f\n(%s) --> \"%s\"", initial() ? "during initialization " : "", data->localData[0]->timeValue, assert_cond, MMC_STRINGDATA(tmpMeta4));
        }
      }
      tmp5 = 1;
    }
  }
  TRACE_POP
}

/*
equation index: 32
type: ALGORITHM

  assert($con$CONSTR >= -1.25 and $con$CONSTR <= 1000.0, "Variable violating min/max constraint: -1.25 <= $con$CONSTR <= 1000.0, has value: " + String($con$CONSTR, "g"));
*/
void include_test_eqFunction_32(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,32};
  modelica_boolean tmp6;
  modelica_boolean tmp7;
  static const MMC_DEFSTRINGLIT(tmp8,82,"Variable violating min/max constraint: -1.25 <= $con$CONSTR <= 1000.0, has value: ");
  modelica_string tmp9;
  modelica_metatype tmpMeta10;
  static int tmp11 = 0;
  if(!tmp11)
  {
    tmp6 = GreaterEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[15]] /* $con$CONSTR OPT_CONSTR */),-1.25);
    tmp7 = LessEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[15]] /* $con$CONSTR OPT_CONSTR */),1000.0);
    if(!(tmp6 && tmp7))
    {
      tmp9 = modelica_real_to_modelica_string_format((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[15]] /* $con$CONSTR OPT_CONSTR */), (modelica_string) mmc_strings_len1[103]);
      tmpMeta10 = stringAppend(MMC_REFSTRINGLIT(tmp8),tmp9);
      {
        const char* assert_cond = "($con$CONSTR >= -1.25 and $con$CONSTR <= 1000.0)";
        if (data->simulationInfo->noThrowAsserts) {
          FILE_INFO info = {"",0,0,0,0,0};
          infoStreamPrintWithEquationIndexes(OMC_LOG_ASSERT, info, 0, equationIndexes, "The following assertion has been violated %sat time %f\n(%s) --> \"%s\"", initial() ? "during initialization " : "", data->localData[0]->timeValue, assert_cond, MMC_STRINGDATA(tmpMeta10));
        } else {
          FILE_INFO info = {"",0,0,0,0,0};
          omc_assert_warning_withEquationIndexes(info, equationIndexes, "The following assertion has been violated %sat time %f\n(%s) --> \"%s\"", initial() ? "during initialization " : "", data->localData[0]->timeValue, assert_cond, MMC_STRINGDATA(tmpMeta10));
        }
      }
      tmp11 = 1;
    }
  }
  TRACE_POP
}

/*
equation index: 33
type: ALGORITHM

  assert(FINALCONSTR >= -1.25 and FINALCONSTR <= 1.25, "Variable violating min/max constraint: -1.25 <= FINALCONSTR <= 1.25, has value: " + String(FINALCONSTR, "g"));
*/
void include_test_eqFunction_33(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,33};
  modelica_boolean tmp12;
  modelica_boolean tmp13;
  static const MMC_DEFSTRINGLIT(tmp14,80,"Variable violating min/max constraint: -1.25 <= FINALCONSTR <= 1.25, has value: ");
  modelica_string tmp15;
  modelica_metatype tmpMeta16;
  static int tmp17 = 0;
  if(!tmp17)
  {
    tmp12 = GreaterEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */),-1.25);
    tmp13 = LessEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */),1.25);
    if(!(tmp12 && tmp13))
    {
      tmp15 = modelica_real_to_modelica_string_format((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[6]] /* FINALCONSTR variable */), (modelica_string) mmc_strings_len1[103]);
      tmpMeta16 = stringAppend(MMC_REFSTRINGLIT(tmp14),tmp15);
      {
        const char* assert_cond = "(FINALCONSTR >= -1.25 and FINALCONSTR <= 1.25)";
        if (data->simulationInfo->noThrowAsserts) {
          FILE_INFO info = {"/home/linus/Projects/Optimization/src/modelica/example/include_test.mo",11,3,11,92,0};
          infoStreamPrintWithEquationIndexes(OMC_LOG_ASSERT, info, 0, equationIndexes, "The following assertion has been violated %sat time %f\n(%s) --> \"%s\"", initial() ? "during initialization " : "", data->localData[0]->timeValue, assert_cond, MMC_STRINGDATA(tmpMeta16));
        } else {
          FILE_INFO info = {"/home/linus/Projects/Optimization/src/modelica/example/include_test.mo",11,3,11,92,0};
          omc_assert_warning_withEquationIndexes(info, equationIndexes, "The following assertion has been violated %sat time %f\n(%s) --> \"%s\"", initial() ? "during initialization " : "", data->localData[0]->timeValue, assert_cond, MMC_STRINGDATA(tmpMeta16));
        }
      }
      tmp17 = 1;
    }
  }
  TRACE_POP
}
/* function to check assert after a step is done */
OMC_DISABLE_OPT
int include_test_checkForAsserts(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  include_test_eqFunction_31(data, threadData);

  include_test_eqFunction_32(data, threadData);

  include_test_eqFunction_33(data, threadData);
  
  TRACE_POP
  return 0;
}

#if defined(__cplusplus)
}
#endif

