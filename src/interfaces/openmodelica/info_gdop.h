#ifndef OPT_OM_INFO_GDOP_H
#define OPT_OM_INFO_GDOP_H

#include "simulation_data.h"

#include <base/nlp_structs.h>
#include <base/fixed_vector.h>
#include <memory>

/* foward decl */
struct ExchangeJacobians;

struct InfoGDOP {
    int x_size;
    int u_size;
    int p_size;
    int xu_size;
    int f_size;
    int g_size;
    int r_size;

    bool mayer_exists;
    bool lagrange_exists;

    modelica_real* __address_mayer_real_vars;
    modelica_real* __address_lagrange_real_vars;

    const int index_x_real_vars = 0;
    int index_der_x_real_vars;
    int index_mayer_real_vars = -1;
    int index_lagrange_real_vars = -1;
    FixedVector<int> u_indices_real_vars;
    int index_g_real_vars;
    int index_r_real_vars;

    std::unique_ptr<ExchangeJacobians> exc_jac;
};

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

struct ExchangeJacobians {
    JACOBIAN* A;
    JACOBIAN* B;
    JACOBIAN* C;
    JACOBIAN* D;

    bool A_exists;
    bool B_exists;
    bool C_exists;
    bool D_exists;

    Exchange_COO_CSC A_coo;
    Exchange_COO_CSC B_coo;
    Exchange_COO_CSC C_coo;
    Exchange_COO_CSC D_coo;

    ExchangeJacobians(DATA* data, threadData_t* threadData, InfoGDOP& info);
};

#endif // OPT_OM_INFO_GDOP_H
