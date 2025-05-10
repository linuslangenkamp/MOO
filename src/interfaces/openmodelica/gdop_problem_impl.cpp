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

/**
 * @brief Constructs and initializes OpenModelica Jacobians and their COO mappings for optimization.
 *
 * This constructor sets up Jacobian matrices (A, B, C, D) using OpenModelica's internal structures
 * and converts them from CSC (Compressed Sparse Column) format to COO (Coordinate) format.
 * The COO format is required by the optimization, which expects different symbolic block
 * orderings and row arrangements than OpenModelica provides.
 *
 * The OM matrices represent the following derivatives in CSC:
 *   - A = (f_{xu})^T
 *   - B = (f_{xu}, L_{xu}, g_{xu})^T
 *   - C = (f_{xu}, L_{xu}, M_{xu}, g_{xu})^T
 *   - D = (r_{xu})^T
 *
 * However, OPT expects the COO format with:
 *   1. [L, f, g] — i.e. Lagrange terms first in the row ordering of B
 *   2. [M, r]    — i.e. Mayer terms followed by boundary constraints in D
 *
 * Because of this mismatch between OpenModelica's CSC block structure and OPT's COO expectations,
 * the following transformation steps are performed:
 *
 * 1. **Initialize OpenModelica Jacobian objects** using their dedicated function pointers:
 *    These include A, B, C, and D — each initialized via `initialAnalyticJacobianX(...)`.
 *
 * 2. **Convert each Jacobian from CSC to COO** using `Exchange_COO_CSC::from_csc(...)`:
 *    - **A_coo**: a direct conversion from CSC to COO.
 *    - **B_coo**: optionally reorders the Lagrange row (`L_{xu}`) to appear first if `lagrange_exists`.
 *    - **C_coo**: optionally reorders the Mayer row (`M_{xu}`) to appear first if `mayer_exists`.
 *
 * 3. **Extract the Mayer sub-block** from `C_coo` using `RowExchange_COO_CSC::extract_row(...)`.
 *    The reordered Mayer row is now at row index 0 and is extracted into `C_mayer_coo`, thus
 *    forming the first part of OPT’s [M, r] structure.
 *
 * 4. **Construct `D_coo`** after `C_mayer_coo`, taking into account the already-extracted
 *    Mayer term. The D matrix contains the `r_{xu}` rows, but for consistency with the optimizer’s
 *    expected [M, r] structure, we **shift the value indices** using the `nnz_offset`
 *    passed into `from_csc(...)`. This ensures correct alignment and value mapping for `[*, r]`.
 *
 * The permutation arrays `csc_to_coo` and `coo_to_csc` constructed during each step allow consistent
 * mapping of values between the original CSC (OpenModelica) and the reordered COO (optimizer) forms.
 * 
 * Hopefully, we soon get block D = (r_{xu}, M_{xu})^T or (M_{xu}, r_{xu})^T, making it easier!
 */
Jacobians_OM::Jacobians_OM(DATA* data, threadData_t* threadData, bool mayer_exists, bool lagrange_exists, int x_size) :
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
                                     lagrange_exists ? x_size : -1)),
    C_coo(Exchange_COO_CSC::from_csc((int*)C->sparsePattern->leadindex, (int*)C->sparsePattern->index,
                                     (int)C->sizeCols, (int)C->sparsePattern->numberOfNonZeros,
                                     mayer_exists ? x_size + (int)(lagrange_exists) : -1)),
    C_mayer_coo(mayer_exists ? std::make_unique<RowExchange_COO_CSC>(RowExchange_COO_CSC::extract_row(C_coo, 0)) : NULL),
    D_coo(Exchange_COO_CSC::from_csc((int*)D->sparsePattern->leadindex, (int*)D->sparsePattern->index,
                                     (int)D->sizeCols, (int)D->sparsePattern->numberOfNonZeros,
                                     -1, C_mayer_coo->nnz)) {
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

    print_jacobian_sparsity(jac_om.D, true, "D");

    jac_om.D_coo.row.print();
    jac_om.D_coo.col.print();

    //     print_jacobian_sparsity(jac_om.C, true, "C");

    //jac_om.C_mayer_coo->col.print();
    //jac_om.C_mayer_coo->csc_to_coo.print();
    /*
    print_jacobian_sparsity(jac_om.A, true, "A");
    print_jacobian_sparsity(jac_om.D, true, "D");
    */

    return NULL;
}
