#include "info_gdop.h"

namespace OpenModelica {

InfoGDOP::InfoGDOP(DATA* data, threadData_t* threadData, int argc, char** argv) :
                   data(data), threadData(threadData), argc(argc), argv(argv) {}

void InfoGDOP::set_time_horizon(int steps) {
    model_start_time = data->simulationInfo->startTime;
    model_stop_time = data->simulationInfo->stopTime;
    tf = model_stop_time - model_start_time;
    intervals = (int)(round(tf/data->simulationInfo->stepSize));
    stages = steps;
}

void InfoGDOP::set_omc_flags(NLP::NLPSolverSettings& nlp_solver_settings) {
    char* cflags = (char*)omc_flagValue[FLAG_OPTIMIZER_NP];
    set_time_horizon(cflags ? atoi(cflags) : 3);

    nlp_solver_settings.set(NLP::Option::Tolerance, data->simulationInfo->tolerance);

    // Linear solver
    cflags = (char*)omc_flagValue[FLAG_LS_IPOPT];
    if (cflags) {
        std::string opt(cflags);
        std::string lower;
        std::transform(opt.begin(), opt.end(), std::back_inserter(lower), ::tolower);

        using LS = NLP::LinearSolverOption;
        if (lower == "mumps") {
            nlp_solver_settings.set(NLP::Option::LinearSolver, LS::MUMPS);
        } else if (lower == "ma27") {
            nlp_solver_settings.set(NLP::Option::LinearSolver, LS::MA27);
        } else if (lower == "ma57") {
            nlp_solver_settings.set(NLP::Option::LinearSolver, LS::MA57);
        } else if (lower == "ma77") {
            nlp_solver_settings.set(NLP::Option::LinearSolver, LS::MA77);
        } else if (lower == "ma86") {
            nlp_solver_settings.set(NLP::Option::LinearSolver, LS::MA86);
        } else if (lower == "ma97") {
            nlp_solver_settings.set(NLP::Option::LinearSolver, LS::MA97);
        } else {
            LOG_WARNING("Unsupported linear solver option: %s", cflags);
        }
    }

    // Maximum iterations
    cflags = (char*)omc_flagValue[FLAG_IPOPT_MAX_ITER];
    if (cflags) {
        try {
            nlp_solver_settings.set(NLP::Option::Iterations, std::stoi(cflags));
        } catch (...) {
            LOG_WARNING("Invalid integer for Iterations: %s", cflags);
        }
    }

    // Hessian option
    cflags = (char*)omc_flagValue[FLAG_IPOPT_HESSE];
    if (cflags) {
        std::string opt(cflags);
        std::string lower;
        std::transform(opt.begin(), opt.end(), std::back_inserter(lower), ::tolower);

        using H = NLP::HessianOption;
        if (lower == "bfgs" || lower == "lbfgs") {
            nlp_solver_settings.set(NLP::Option::Hessian, H::LBFGS);
        } else if (lower == "const" || lower == "qp") {
            nlp_solver_settings.set(NLP::Option::Hessian, H::CONST);
        } else if (lower == "exact") {
            nlp_solver_settings.set(NLP::Option::Hessian, H::Exact);
        } else {
            LOG_WARNING("Unsupported Hessian option: %s (use LBFGS, QP, or Exact)", cflags);
        }
    }
}

ExchangeJacobians::ExchangeJacobians(InfoGDOP& info) :
    /* set OpenModelica Jacobian ptrs, allocate memory, initilization of A, B, C, D */
    A(&(info.data->simulationInfo->analyticJacobians[info.data->callback->INDEX_JAC_A])),
    B(&(info.data->simulationInfo->analyticJacobians[info.data->callback->INDEX_JAC_B])),
    C(&(info.data->simulationInfo->analyticJacobians[info.data->callback->INDEX_JAC_C])),
    D(&(info.data->simulationInfo->analyticJacobians[info.data->callback->INDEX_JAC_D])),
    A_exists((bool)(info.data->callback->initialAnalyticJacobianA(info.data, info.threadData, A) == 0)),
    B_exists((bool)(info.data->callback->initialAnalyticJacobianB(info.data, info.threadData, B) == 0)),
    C_exists((bool)(info.data->callback->initialAnalyticJacobianC(info.data, info.threadData, C) == 0)),
    D_exists((bool)(info.data->callback->initialAnalyticJacobianD(info.data, info.threadData, D) == 0)),

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

ExchangeHessians::ExchangeHessians(InfoGDOP& info) :
    A(generateHessianPattern(info.exc_jac->A)),
    B(generateHessianPattern(info.exc_jac->B)),
    C(generateHessianPattern(info.exc_jac->C)),
    D(generateHessianPattern(info.exc_jac->D)),

    A_exists(info.exc_jac->A_exists),
    B_exists(info.exc_jac->B_exists),
    C_exists(info.exc_jac->C_exists),
    D_exists(info.exc_jac->D_exists),

    A_extr(!A_exists ? nullptr : initExtrapolationData(A->lnnz, 5)),
    B_extr(!B_exists ? nullptr : initExtrapolationData(B->lnnz, 5)),
    C_extr(!C_exists ? nullptr : initExtrapolationData(C->lnnz, 5)),
    D_extr(!D_exists ? nullptr : initExtrapolationData(D->lnnz, 5)),

    A_buffer(FixedVector<modelica_real>(!A_exists ? 0 : A->lnnz)),
    B_buffer(FixedVector<modelica_real>(!B_exists ? 0 : B->lnnz)),
    C_buffer(FixedVector<modelica_real>(!C_exists ? 0 : C->lnnz)),
    D_buffer(FixedVector<modelica_real>(!D_exists ? 0 : D->lnnz)),

    A_lambda(FixedVector<modelica_real>(!A_exists ? 0 : A->numFuncs)),
    B_lambda(FixedVector<modelica_real>(!B_exists ? 0 : B->numFuncs)),
    C_lambda(FixedVector<modelica_real>(!C_exists ? 0 : C->numFuncs)),
    D_lambda(FixedVector<modelica_real>(!D_exists ? 0 : D->numFuncs)),

    A_args{info.data, info.threadData, A, info.u_indices_real_vars.raw(), A_lambda.raw(), nullptr},
    B_args{info.data, info.threadData, B, info.u_indices_real_vars.raw(), B_lambda.raw(), nullptr},
    C_args{info.data, info.threadData, C, info.u_indices_real_vars.raw(), C_lambda.raw(), nullptr},
    D_args{info.data, info.threadData, D, info.u_indices_real_vars.raw(), D_lambda.raw(), nullptr} {
    /* attach global, heap allocated C structs to auto free */
    info.auto_free.attach({A, B, C, D}, freeHessianPattern);
    info.auto_free.attach({A_extr, B_extr, C_extr, D_extr}, freeExtrapolationData);
}

} // namespace OpenModelica
