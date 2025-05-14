/* Asserts */
#include "include_test_model.h"
#if defined(__cplusplus)
extern "C" {
#endif


/*
equation index: 13
type: ALGORITHM

  assert($con$g >= 0.0 and $con$g <= 1.5, "Variable violating min/max constraint: 0.0 <= $con$g <= 1.5, has value: " + String($con$g, "g"));
*/
void include_test_eqFunction_13(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,13};
  modelica_boolean tmp0;
  modelica_boolean tmp1;
  static const MMC_DEFSTRINGLIT(tmp2,72,"Variable violating min/max constraint: 0.0 <= $con$g <= 1.5, has value: ");
  modelica_string tmp3;
  modelica_metatype tmpMeta4;
  static int tmp5 = 0;
  if(!tmp5)
  {
    tmp0 = GreaterEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* $con$g OPT_CONSTR */),0.0);
    tmp1 = LessEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* $con$g OPT_CONSTR */),1.5);
    if(!(tmp0 && tmp1))
    {
      tmp3 = modelica_real_to_modelica_string_format((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[7]] /* $con$g OPT_CONSTR */), (modelica_string) mmc_strings_len1[103]);
      tmpMeta4 = stringAppend(MMC_REFSTRINGLIT(tmp2),tmp3);
      {
        const char* assert_cond = "($con$g >= 0.0 and $con$g <= 1.5)";
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
/* function to check assert after a step is done */
OMC_DISABLE_OPT
int include_test_checkForAsserts(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH

  include_test_eqFunction_13(data, threadData);
  
  TRACE_POP
  return 0;
}

#if defined(__cplusplus)
}
#endif

