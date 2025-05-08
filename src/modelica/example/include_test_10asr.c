/* Asserts */
#include "include_test_model.h"
#if defined(__cplusplus)
extern "C" {
#endif


/*
equation index: 27
type: ALGORITHM

  assert($con$CONSTR >= -1.25 and $con$CONSTR <= 1.25, "Variable violating min/max constraint: -1.25 <= $con$CONSTR <= 1.25, has value: " + String($con$CONSTR, "g"));
*/
void include_test_eqFunction_27(DATA *data, threadData_t *threadData)
{
  TRACE_PUSH
  const int equationIndexes[2] = {1,27};
  modelica_boolean tmp0;
  modelica_boolean tmp1;
  static const MMC_DEFSTRINGLIT(tmp2,80,"Variable violating min/max constraint: -1.25 <= $con$CONSTR <= 1.25, has value: ");
  modelica_string tmp3;
  modelica_metatype tmpMeta4;
  static int tmp5 = 0;
  if(!tmp5)
  {
    tmp0 = GreaterEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* $con$CONSTR OPT_CONSTR */),-1.25);
    tmp1 = LessEq((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* $con$CONSTR OPT_CONSTR */),1.25);
    if(!(tmp0 && tmp1))
    {
      tmp3 = modelica_real_to_modelica_string_format((data->localData[0]->realVars[data->simulationInfo->realVarsIndex[14]] /* $con$CONSTR OPT_CONSTR */), (modelica_string) mmc_strings_len1[103]);
      tmpMeta4 = stringAppend(MMC_REFSTRINGLIT(tmp2),tmp3);
      {
        const char* assert_cond = "($con$CONSTR >= -1.25 and $con$CONSTR <= 1.25)";
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

  include_test_eqFunction_27(data, threadData);
  
  TRACE_POP
  return 0;
}

#if defined(__cplusplus)
}
#endif

