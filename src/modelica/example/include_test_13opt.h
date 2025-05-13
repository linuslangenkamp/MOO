#if defined(__cplusplus)
  extern "C" {
#endif
  int include_test_mayer(DATA* data, modelica_real** res, short*);
  int include_test_lagrange(DATA* data, modelica_real** res, short *, short *);
  int include_test_getInputVarIndicesInOptimization(DATA* data, int* input_var_indices);
  int include_test_pickUpBoundsForInputsInOptimization(DATA* data, modelica_real* min, modelica_real* max, modelica_real*nominal, modelica_boolean *useNominal, char ** name, modelica_real * start, modelica_real * startTimeOpt);
  int include_test_setInputData(DATA *data, const modelica_boolean file);
  int include_test_getTimeGrid(DATA *data, modelica_integer * nsi, modelica_real**t);
#if defined(__cplusplus)
}
#endif