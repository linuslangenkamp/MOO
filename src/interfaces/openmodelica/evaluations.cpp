#include "evaluations.h"

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
 * 3. `D_coo` contains only `r_{xu}`, but OPT expects Mayer terms (`M_{xu}`) first.
 *    If `mayer_exists`, we set nnz_offset of D_coo to `C_coo.row_nnz(0)` to ensure `[M, r]` order.
 *    See `nnz_offset` in `Exchange_COO_CSC` for further info.
 *
 * The permutation arrays `csc_to_coo` and `coo_to_csc` constructed during each step allow consistent
 * mapping of values between the original CSC (OpenModelica) and the reordered COO (optimizer) forms.
 * 
 * Hopefully, we soon get block D = (r_{xu}, M_{xu})^T or (M_{xu}, r_{xu})^T, making less involved!
 */
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
    D_coo(Exchange_COO_CSC::from_csc((int*)D->sparsePattern->leadindex, (int*)D->sparsePattern->index,
                                     (int)D->sizeCols, (int)D->sparsePattern->numberOfNonZeros,
                                     -1, info.mayer_exists ? C_coo.row_nnz(0) : 0)) {
}

void ExchangeJacobians::init_jac(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr) {
    init_jac_lfg(data, threadData, info, lfg);
    init_jac_mr(data, threadData, info, mr);
}

void ExchangeJacobians::init_jac_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg) {
    for (int nz = 0; nz < B_coo.nnz; nz++) {
        int row = B_coo.row[nz];
        int col = B_coo.col[nz];
        int csc_buffer_entry = B_coo.coo_to_csc(nz); // jac_buffer == OpenModelica CSC buffer!
        if (col < info.x_size) {
            lfg[row].jac.dx.push_back(JacobianSparsity{col, csc_buffer_entry});
        }
        else if (col < info.xu_size) {
            lfg[row].jac.du.push_back(JacobianSparsity{col - info.x_size, csc_buffer_entry});
        }
        else {
            lfg[row].jac.dp.push_back(JacobianSparsity{col - info.xu_size, csc_buffer_entry});
        }
    }
}

void ExchangeJacobians::init_jac_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr) {

}