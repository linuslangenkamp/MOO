#include "gdop_problem_impl.h"

// Dummy implementation of FullSweep_OM
FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds, 
                           bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size)
    : FullSweep(std::move(lfg), mesh, std::move(g_bounds), has_lagrange, f_size, g_size, x_size, u_size, p_size) {
    // Dummy constructor (does nothing)
}
void FullSweep_OM::callback_eval(const F64* xu_nlp, const F64* p) {
    // Dummy evaluation (does nothing)
}

void FullSweep_OM::callback_jac(const F64* xu_nlp, const F64* p) {
    // Dummy Jacobian (does nothing)
}

void FullSweep_OM::callback_hes(const F64* xu_nlp, const F64* p) {
    // Dummy Hessian (does nothing)
}

// Dummy implementation of BoundarySweep_OM
BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                                   bool has_mayer, int r_size, int x_size, int p_size)
    : BoundarySweep(std::move(mr), mesh, std::move(r_bounds), has_mayer, r_size, x_size, p_size) {
    // Dummy constructor (does nothing)
}

void BoundarySweep_OM::callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    // Dummy evaluation (does nothing)
}

void BoundarySweep_OM::callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    // Dummy Jacobian (does nothing)
}

void BoundarySweep_OM::callback_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    // Dummy Hessian (does nothing)
}

/* initialize OpenModelica Jacobian with sparsity, seeds, etc. */
Jacobians_OM::Jacobians_OM(DATA* data, threadData_t* threadData, bool mayer_exists, bool lagrange_exists, int x_size) :
    // set ptrs
    A(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_A])),
    B(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_B])),
    C(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_C])),
    D(&(data->simulationInfo->analyticJacobians[data->callback->INDEX_JAC_D])),

    // allocate memory and initilization of Jacobians
    A_exists((bool)(data->callback->initialAnalyticJacobianA(data, threadData, A) == 0)),
    B_exists((bool)(data->callback->initialAnalyticJacobianB(data, threadData, B) == 0)),
    C_exists((bool)(data->callback->initialAnalyticJacobianC(data, threadData, C) == 0)),
    D_exists((bool)(data->callback->initialAnalyticJacobianD(data, threadData, D) == 0)),

    // remember csc -> coo format + the conversion arrays
    B_coo(Exchange_COO_CSC::from_csc((int*)B->sparsePattern->leadindex, (int*)B->sparsePattern->index,
                                     (int)B->sizeCols, (int)B->sparsePattern->numberOfNonZeros)),
    C_coo(Exchange_COO_CSC::from_csc((int*)C->sparsePattern->leadindex, (int*)C->sparsePattern->index,
                                     (int)C->sizeCols, (int)C->sparsePattern->numberOfNonZeros)),
    D_coo(Exchange_COO_CSC::from_csc((int*)D->sparsePattern->leadindex, (int*)D->sparsePattern->index,
                                     (int)D->sizeCols, (int)D->sparsePattern->numberOfNonZeros)) {
    if (mayer_exists) {
        C_mayer_coo = std::make_unique<RowExchange_COO_CSC>(RowExchange_COO_CSC::extract_row(C_coo, (int)lagrange_exists + x_size));
    }
}

Problem* create_gdop_om(DATA* data, threadData_t* threadData, Mesh& mesh) {
    /* variable sizes */
    int x_size = data->modelData->nStates;
    int u_size = data->modelData->nInputVars;
    int p_size = 0; // TODO: add this feature

    int g_size = data->modelData->nOptimizeConstraints;
    int r_size = data->modelData->nOptimizeFinalConstraints; // TODO: add *generic boundary* constraints later also at t=t0

    /* variable bounds */
    FixedVector<Bounds> x_bounds(x_size);
    FixedVector<Bounds> u_bounds(u_size);
    FixedVector<Bounds> p_bounds(p_size);

    for (int x = 0; x < x_size; x++) {
        x_bounds[x].lb = data->modelData->realVarsData[x].attribute.min;
        x_bounds[x].ub = data->modelData->realVarsData[x].attribute.max;
    }

    /* new generated function getInputVarIndices, just fills a index list */
    auto u_indices_real_vars = std::make_unique<int>(u_size);
    data->callback->getInputVarIndices(data, u_indices_real_vars.get());
    for (int u = 0; u < u_size; u++) {
        int u_index = *(u_indices_real_vars.get() + u);
        u_bounds[u].lb = data->modelData->realVarsData[u_index].attribute.min;
        u_bounds[u].ub = data->modelData->realVarsData[u_index].attribute.max;
    }

    /* constraint bounds */
    FixedVector<Bounds> g_bounds(g_size);
    FixedVector<Bounds> r_bounds(r_size);

    int first_index_g = data->modelData->nVariablesReal - (g_size + r_size);
    int first_index_r = data->modelData->nVariablesReal - r_size;
    for (int g = 0; g < g_size; g++) {
        g_bounds[g].lb = data->modelData->realVarsData[first_index_g + g].attribute.min;
        g_bounds[g].ub = data->modelData->realVarsData[first_index_g + g].attribute.max;
    }

    for (int r = 0; r < r_size; r++) {
        r_bounds[r].lb = data->modelData->realVarsData[first_index_r + r].attribute.min;
        r_bounds[r].ub = data->modelData->realVarsData[first_index_r + r].attribute.max;
    }

    /* for now we ignore xf fixed (need some steps in Backend to detect)
     * and also ignore x0 non fixed, since too complicated
     * => assume x(t_0) = x0 fixed, x(t_f) free to r constraint / maybe the old BE can do that already?!
     * option: generate fixed final states individually
     */
    FixedVector<std::optional<F64>> x0_fixed(x_size);
    FixedVector<std::optional<F64>> xf_fixed(x_size);

    /* set *fixed* initial, final states */
    for (int x = 0; x < x_size; x++) {
        x0_fixed[x] = data->modelData->realVarsData[x].attribute.start;
    }

    // TODO: figure this out; we get some derivative ptrs i guess?
    short der_index_mayer_realVars = -1;
    short der_indices_lagrange_realVars[2] = {-1, -1};
    F64* address_mayer_realVars;
    F64* address_lagrange_realVars;

    // this is really ugly IMO, fix this when ready for master!
    bool mayer_exists = (data->callback->mayer(data, &address_mayer_realVars, &der_index_mayer_realVars) >= 0);
    bool lagrange_exists = (data->callback->lagrange(data, &address_lagrange_realVars, &der_indices_lagrange_realVars[0], &der_indices_lagrange_realVars[1]) >= 0);

    /* create functions and bounds */
    FixedVector<FunctionMR> mr((int)mayer_exists + r_size);
    FixedVector<FunctionLFG> lfg((int)lagrange_exists + x_size + g_size);

    /* create sparsity patterns */
    Jacobians_OM jac_om(data, threadData, mayer_exists, lagrange_exists, x_size);

    print_jacobian_sparsity(jac_om.B, true, "B");
    print_jacobian_sparsity(jac_om.C, true, "C");

    jac_om.C_mayer_coo->col.print();
    jac_om.C_mayer_coo->coo_to_csc.print();
    /*
    print_jacobian_sparsity(jac_om.A, true, "A");
    print_jacobian_sparsity(jac_om.D, true, "D");
    */

    return NULL;
}
