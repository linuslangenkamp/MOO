#ifndef OPT_OM_INFO_GDOP_H
#define OPT_OM_INFO_GDOP_H

#include "simulation_data.h"
#include "simulation/solver/gbode_main.h"
#include "simulation/solver/external_input.h"

#include "sim_runtime_ext.h"

#include <base/nlp_structs.h>
#include <base/fixed_vector.h>
#include <nlp/solvers/nlp_solver_settings.h>

#include <memory>

namespace OpenModelica {

/* foward decl */
struct ExchangeJacobians;
struct ExchangeHessians;

struct InfoGDOP {
    /* custom attaching and auto freeing of C-style mallocs / callocs */
    AutoFree auto_free;

    DATA* data;               // pointer to OM data object
    threadData_t* threadData; // pointer to OM threadData object
    int argc;                 // command-line arg count OM
    char** argv;              // command-line args OM

    /* problem sizes */
    int x_size;
    int u_size;
    int p_size;
    int xu_size;
    int f_size;
    int g_size;
    int r_size;

    /* objective structure */
    bool mayer_exists;
    bool lagrange_exists;

    /* addresses in realVars */
    modelica_real* address_mayer_real_vars;
    modelica_real* address_lagrange_real_vars;

    /* realVars variables indices */
    const int index_x_real_vars = 0;
    int index_der_x_real_vars;
    int index_mayer_real_vars = -1;
    int index_lagrange_real_vars = -1;
    FixedVector<int> u_indices_real_vars;
    int index_g_real_vars;
    int index_r_real_vars;

    /* exchange format for Jacobians OM <-> OPT */
    std::unique_ptr<ExchangeJacobians> exc_jac;

    /* numerical augmented Hessians */
    std::unique_ptr<ExchangeHessians> exc_hes;

    /* time horizon */
    f64 model_start_time; // model start_time
    f64 model_stop_time;  // model stop_time
    f64 tf;               // tf = start - stop, since t0 = 0 for OPT
    int intervals;        // model interval count
    int stages;           // stage count

    InfoGDOP(DATA* data, threadData_t* threadData, int argc, char** argv);

    void set_omc_flags(NLP::NLPSolverSettings& nlp_solver_settings);
    void set_time_horizon(int steps);
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
 * 4. Optional Jacobian buffer with size nnz(Matrix) are also included and can be used if no in-place
 *    buffer is applicable
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

    FixedVector<modelica_real> A_buffer;
    FixedVector<modelica_real> B_buffer;
    FixedVector<modelica_real> C_buffer;
    FixedVector<modelica_real> D_buffer;

    ExchangeJacobians(InfoGDOP& info);
};

struct ExchangeHessians {
    HESSIAN_PATTERN* A;
    HESSIAN_PATTERN* B;
    HESSIAN_PATTERN* C;
    HESSIAN_PATTERN* D;

    bool A_exists;
    bool B_exists;
    bool C_exists;
    bool D_exists;

    ExtrapolationData* A_extr;
    ExtrapolationData* B_extr;
    ExtrapolationData* C_extr;
    ExtrapolationData* D_extr;

    /* some workspace memory for output of Hessian, use if needed */
    FixedVector<modelica_real> A_buffer;
    FixedVector<modelica_real> B_buffer;
    FixedVector<modelica_real> C_buffer;
    FixedVector<modelica_real> D_buffer;

    /* some workspace memory for dual multiplicators, use if needed */
    FixedVector<modelica_real> A_lambda;
    FixedVector<modelica_real> B_lambda;
    FixedVector<modelica_real> C_lambda;
    FixedVector<modelica_real> D_lambda;

    /* wrapper args for Richardson extrapolation */
    HessianFiniteDiffArgs A_args;
    HessianFiniteDiffArgs B_args;
    HessianFiniteDiffArgs C_args;
    HessianFiniteDiffArgs D_args;

    /* mapping of OM C and D HESSIAN_PATTERN indices -> Mr buffer indices, will be set in 'init_hes_mr()' */
    FixedVector<std::pair<int, int>> C_to_Mr_buffer;
    FixedVector<std::pair<int, int>> D_to_Mr_buffer;

    ExchangeHessians(InfoGDOP& info);
};

} // namespace OpenModelica

#endif // OPT_OM_INFO_GDOP_H
