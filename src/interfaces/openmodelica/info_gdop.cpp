#include "info_gdop.h"

void InfoGDOP::set_time_horizon(DATA* data, int collocation) {
    start_time = data->simulationInfo->startTime;
    stop_time = data->simulationInfo->stopTime;
    tf = stop_time - start_time;
    intervals = (int)(round(tf/data->simulationInfo->stepSize));
    stages = collocation;
}

ExchangeJacobians::ExchangeJacobians(DATA* data, threadData_t* threadData, InfoGDOP& info) :
    /* set OpenModelica Jacobian ptrs, allocate memory, initilization of A, B, C, D */
    A(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A])),
    B(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_B])),
    C(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_C])),
    D(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_D])),
    A_exists((bool)(data->callback->initialAnalyticJacobianA(data, threadData, A) == 0)),
    B_exists((bool)(data->callback->initialAnalyticJacobianB(data, threadData, B) == 0)),
    C_exists((bool)(data->callback->initialAnalyticJacobianC(data, threadData, C) == 0)),
    D_exists((bool)(data->callback->initialAnalyticJacobianD(data, threadData, D) == 0)),

    /* create COO sparsity and CSC(OM) <-> COO(OPT, reordered) mappings */
    A_coo(Exchange_COO_CSC::from_csc((int*)A->sparsePattern->leadindex, (int*)A->sparsePattern->index,
                                     (int)A->sizeCols, (int)A->sparsePattern->numberOfNonZeros)),
    B_coo(Exchange_COO_CSC::from_csc((int*)B->sparsePattern->leadindex, (int*)B->sparsePattern->index,
                                     (int)B->sizeCols, (int)B->sparsePattern->numberOfNonZeros,
                                     info.lagrange_exists ? info.x_size : -1)),
    C_coo(Exchange_COO_CSC::from_csc((int*)C->sparsePattern->leadindex, (int*)C->sparsePattern->index,
                                     (int)C->sizeCols, (int)C->sparsePattern->numberOfNonZeros,
                                     info.mayer_exists ? info.x_size + (int)(info.lagrange_exists) : -1)),
    D_coo(info.r_size != 0 ? Exchange_COO_CSC::from_csc((int*)D->sparsePattern->leadindex, (int*)D->sparsePattern->index,
                                     (int)D->sizeCols, (int)D->sparsePattern->numberOfNonZeros,
                                     -1, info.mayer_exists ? C_coo.row_nnz(0) : 0) : Exchange_COO_CSC()),
    
    /* create optional buffers, use when in-place buffers are no option */
    A_buffer(FixedVector<modelica_real>(A_coo.nnz)),
    B_buffer(FixedVector<modelica_real>(B_coo.nnz)),
    C_buffer(FixedVector<modelica_real>(C_coo.nnz)),
    D_buffer(FixedVector<modelica_real>(D_coo.nnz)) {
}

ExchangeHessians::ExchangeHessians(DATA* data, threadData_t* threadData, InfoGDOP& info) :
    A(__generateHessianPattern(info.exc_jac->A)),
    B(__generateHessianPattern(info.exc_jac->B)),
    C(__generateHessianPattern(info.exc_jac->C)),
    D(__generateHessianPattern(info.exc_jac->D)),

    A_exists(info.exc_jac->A_exists),
    B_exists(info.exc_jac->B_exists),
    C_exists(info.exc_jac->C_exists),
    D_exists(info.exc_jac->D_exists),

    A_extr(!A_exists ? nullptr : __initExtrapolationData(A->lnnz, 5)),
    B_extr(!B_exists ? nullptr : __initExtrapolationData(B->lnnz, 5)),
    C_extr(!C_exists ? nullptr : __initExtrapolationData(C->lnnz, 5)),
    D_extr(!D_exists ? nullptr : __initExtrapolationData(D->lnnz, 5)),

    A_buffer(FixedVector<modelica_real>(!A_exists ? 0 : A->lnnz)),
    B_buffer(FixedVector<modelica_real>(!B_exists ? 0 : B->lnnz)),
    C_buffer(FixedVector<modelica_real>(!C_exists ? 0 : C->lnnz)),
    D_buffer(FixedVector<modelica_real>(!D_exists ? 0 : D->lnnz)),

    A_lambda(FixedVector<modelica_real>(!A_exists ? 0 : A->numFuncs)),
    B_lambda(FixedVector<modelica_real>(!B_exists ? 0 : B->numFuncs)),
    C_lambda(FixedVector<modelica_real>(!C_exists ? 0 : C->numFuncs)),
    D_lambda(FixedVector<modelica_real>(!D_exists ? 0 : D->numFuncs)),

    A_args{data, threadData, A, A_lambda.raw()},
    B_args{data, threadData, B, B_lambda.raw()},
    C_args{data, threadData, C, C_lambda.raw()},
    D_args{data, threadData, D, D_lambda.raw()} {
    /* attach global, heap allocated C structs to auto free */
    info.auto_free.attach({A, B, C, D}, __freeHessianPattern);
    info.auto_free.attach({A_extr, B_extr, C_extr, D_extr}, __freeExtrapolationData);
}
