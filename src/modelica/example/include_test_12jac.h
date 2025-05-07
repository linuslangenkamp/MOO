/* Jacobians */
static const REAL_ATTRIBUTE dummyREAL_ATTRIBUTE = omc_dummyRealAttribute;

#if defined(__cplusplus)
extern "C" {
#endif

/* Jacobian Variables */
#define include_test_INDEX_JAC_H 0
int include_test_functionJacH_column(DATA* data, threadData_t *threadData, JACOBIAN *thisJacobian, JACOBIAN *parentJacobian);
int include_test_initialAnalyticJacobianH(DATA* data, threadData_t *threadData, JACOBIAN *jacobian);


#define include_test_INDEX_JAC_F 1
int include_test_functionJacF_column(DATA* data, threadData_t *threadData, JACOBIAN *thisJacobian, JACOBIAN *parentJacobian);
int include_test_initialAnalyticJacobianF(DATA* data, threadData_t *threadData, JACOBIAN *jacobian);


#define include_test_INDEX_JAC_D 2
int include_test_functionJacD_column(DATA* data, threadData_t *threadData, JACOBIAN *thisJacobian, JACOBIAN *parentJacobian);
int include_test_initialAnalyticJacobianD(DATA* data, threadData_t *threadData, JACOBIAN *jacobian);


#define include_test_INDEX_JAC_C 3
int include_test_functionJacC_column(DATA* data, threadData_t *threadData, JACOBIAN *thisJacobian, JACOBIAN *parentJacobian);
int include_test_initialAnalyticJacobianC(DATA* data, threadData_t *threadData, JACOBIAN *jacobian);


#define include_test_INDEX_JAC_B 4
int include_test_functionJacB_column(DATA* data, threadData_t *threadData, JACOBIAN *thisJacobian, JACOBIAN *parentJacobian);
int include_test_initialAnalyticJacobianB(DATA* data, threadData_t *threadData, JACOBIAN *jacobian);


#define include_test_INDEX_JAC_A 5
int include_test_functionJacA_column(DATA* data, threadData_t *threadData, JACOBIAN *thisJacobian, JACOBIAN *parentJacobian);
int include_test_initialAnalyticJacobianA(DATA* data, threadData_t *threadData, JACOBIAN *jacobian);

#if defined(__cplusplus)
}
#endif

