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
    D_coo(info.r_size != 0 ? Exchange_COO_CSC::from_csc((int*)D->sparsePattern->leadindex, (int*)D->sparsePattern->index,
                                     (int)D->sizeCols, (int)D->sparsePattern->numberOfNonZeros,
                                     -1, info.mayer_exists ? C_coo.row_nnz(0) : 0) : Exchange_COO_CSC()) {
}

/* just enumerate them from 0 ... #lfg - 1 abd 0 ... #mr - 1, the correct placement will be handled in eval */
void init_eval(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr) {
    init_eval_lfg(data, threadData, info, lfg);
    init_eval_mr(data, threadData, info, mr);
}

void init_eval_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg) {
    for (int i = 0; i < lfg.int_size(); i++) {
        lfg[i].buf_index = i;
    }
}

void init_eval_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr) {
    for (int i = 0; i < mr.int_size(); i++) {
        mr[i].buf_index = i;
    }
}

/* since its no problem and conversion from COO <-> CSC has been carried out, we just use the CSC ordering in the Jacobians
 * thus hopefully only using a memcpy to the OPT buffer for Lfg and 2 memcpy's for C(Mayer CSC) and D(r CSC) */
void init_jac(DATA* data, threadData_t* threadData, InfoGDOP& info, ExchangeJacobians& exc_jac, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr) {
    init_jac_lfg(data, threadData, info, exc_jac, lfg);
    init_jac_mr(data, threadData, info, exc_jac, mr);
}

void init_jac_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, ExchangeJacobians& exc_jac, FixedVector<FunctionLFG>& lfg) {
    /* full B Jacobian */
    for (int nz = 0; nz < exc_jac.B_coo.nnz; nz++) {
        int row = exc_jac.B_coo.row[nz];
        int col = exc_jac.B_coo.col[nz];
        int csc_buffer_entry_B = exc_jac.B_coo.coo_to_csc(nz); // jac_buffer == OpenModelica B CSC buffer!
        if (col < info.x_size) {
            lfg[row].jac.dx.push_back(JacobianSparsity{col, csc_buffer_entry_B});
        }
        else if (col < info.xu_size) {
            lfg[row].jac.du.push_back(JacobianSparsity{col - info.x_size, csc_buffer_entry_B});
        }
        else {
            lfg[row].jac.dp.push_back(JacobianSparsity{col - info.xu_size, csc_buffer_entry_B});
        }
    }
}

void init_jac_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, ExchangeJacobians& exc_jac, FixedVector<FunctionMR>& mr) {
    /* M (first row) in C(COO) Jacobian */
    int nz_C = 0;
    if (info.mayer_exists) {
        while (exc_jac.C_coo.row[nz_C] == 0) {
            int col = exc_jac.C_coo.col[nz_C];
            int csc_buffer_entry_C = exc_jac.C_coo.coo_to_csc(nz_C); // jac_buffer == OpenModelica C(row=0) CSC buffer!

            /* for now only final constraints! no parameters, no dx0, esp. no du! */
            if (col < info.x_size) {
                mr[0].jac.dxf.push_back(JacobianSparsity{col, csc_buffer_entry_C});
            }

            nz_C++;
        }
    }

    /* r in D Jacobian */
    int r_start = (int)info.mayer_exists;
    for (int nz_D = 0; nz_D < exc_jac.D_coo.nnz; nz_D++) {
        int row = exc_jac.D_coo.row[nz_D];
        int col = exc_jac.D_coo.col[nz_D];
        int csc_buffer_entry_D = exc_jac.D_coo.coo_to_csc(nz_D); // jac_buffer == OpenModelica D CSC buffer!
        if (col < info.x_size) {
            /* add the Mayer offset, since the values F64* is [M, r] */
            mr[r_start + row].jac.dxf.push_back(JacobianSparsity{col, exc_jac.D_coo.nnz_offset + csc_buffer_entry_D});
        }
    }
}

/* TODO: add me */
void set_parameters(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* p) {
    return;
}

void set_states(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* x_ij) {
    for (int x = 0; x < info.x_size; x++) {
        data->localData[0]->realVars[info.index_x_real_vars + x] = (modelica_real) x_ij[x];
    }
}

void set_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* u_ij) {
    for (int u = 0; u < info.u_size; u++) {
        data->localData[0]->realVars[info.u_indices_real_vars[u]] = (modelica_real) u_ij[info.x_size + u];
    }
}

void set_states_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* xu_ij) {
    for (int x = 0; x < info.x_size; x++) {
        data->localData[0]->realVars[info.index_x_real_vars + x] = (modelica_real) xu_ij[x];
    }

    for (int u = 0; u < info.u_size; u++) {
        data->localData[0]->realVars[info.u_indices_real_vars[u]] = (modelica_real) xu_ij[info.x_size + u];
    }
}

void set_time(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64 t_ij) {
    data->localData[0]->timeValue = (modelica_real) t_ij;
}

void eval_lfg_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_lfg_buffer) {
    int nz = 0;
    /* L */
    if (info.lagrange_exists) {
        eval_lfg_buffer[nz++] = data->localData[0]->realVars[info.index_lagrange_real_vars];
    }
    /* f */
    for (int der_x = 0; der_x < info.f_size; der_x++) {
        eval_lfg_buffer[nz++] = data->localData[0]->realVars[info.index_der_x_real_vars + der_x];
    }
    /* g */
    for (int g = 0; g < info.g_size; g++) {
        eval_lfg_buffer[nz++] = data->localData[0]->realVars[info.index_g_real_vars + g];
    }
}

void eval_mr_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_mr_buffer) {
    int nz = 0;
    /* M */
    if (info.mayer_exists) {
        eval_mr_buffer[nz++] = data->localData[0]->realVars[info.index_mayer_real_vars];
    }
    /* r */
    for (int r = 0; r < info.r_size; r++) {
        eval_mr_buffer[nz++] = data->localData[0]->realVars[info.index_r_real_vars + r];
    }
}

void jac_CD_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<F64>& eval_jac_mr_buffer) {

}
