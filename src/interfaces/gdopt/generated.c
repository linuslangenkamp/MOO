#include <src/interfaces/c/structures.h>

// === problem sizes (compile const) ===

#define X_SIZE 1
#define U_SIZE 1
#define P_SIZE 0
#define RP_SIZE 0

#define R_SIZE 0
#define G_SIZE 0

#define HAS_MAYER false
#define HAS_LAGRANGE true

// === declare global variables (values can be influenced by runtime parameters - _rp) ===

bounds_t globl_x_bounds[X_SIZE] = { { -100.0, 100.0 } };
bounds_t globl_u_bounds[U_SIZE] = { { 0.5, 1.0 } };
bounds_t globl_p_bounds[P_SIZE];

bounds_t globl_g_bounds[G_SIZE];
bounds_t globl_r_bounds[R_SIZE];

optional_value_t globl_x0_fixed[X_SIZE] = { {1.0, true} };
optional_value_t globl_xf_fixed[X_SIZE];

f64 globl_x_nominal[X_SIZE];
f64 globl_u_nominal[U_SIZE];
f64 globl_p_nominal[P_SIZE];

f64 globl_obj_nominal;
f64 globl_f_nominal[X_SIZE];
f64 globl_g_nominal[G_SIZE];
f64 globl_r_nominal[R_SIZE];

f64 _rp[RP_SIZE];

// === sparsity and evaluation structures (compile const) ===

eval_structure_t globl_lfg_eval = {
    .buf_index = (int[]){0, 1}
};

coo_t globl_lfg_jac = {
    .row = (int[]){0, 0, 1, 1},
    .col = (int[]){0, 1, 0, 1},
    .buf_index = (int[]){0, 1, 2, 3},
    .nnz = 4
};

coo_t globl_lfg_lt_aug_hes = {
    .row = (int[]){0, 1},
    .col = (int[]){0, 1},
    .buf_index = (int[]){0, 1},
    .nnz = 2
};

eval_structure_t globl_mr_eval;
coo_t globl_mr_jac;
coo_t globl_mr_lt_aug_hes;

c_problem_t globl_c_problem = {
    .x_size = X_SIZE,
    .u_size = U_SIZE,
    .xu_size = X_SIZE + U_SIZE,
    .p_size = P_SIZE,
    .rp_size = RP_SIZE,
    .r_size = R_SIZE,
    .g_size = G_SIZE,
    .has_mayer = HAS_MAYER,
    .has_lagrange = HAS_LAGRANGE,
    .x_bounds = globl_x_bounds,
    .u_bounds = globl_u_bounds,
    .p_bounds = globl_p_bounds,
    .r_bounds = globl_r_bounds,
    .g_bounds = globl_g_bounds,
    .x0_fixed = globl_x0_fixed,
    .xf_fixed = globl_xf_fixed,
    .x_nominal = globl_x_nominal,
    .u_nominal = globl_u_nominal,
    .p_nominal = globl_p_nominal,
    .obj_nominal = &globl_obj_nominal,
    .f_nominal = globl_f_nominal,
    .g_nominal = globl_g_nominal,
    .r_nominal = globl_r_nominal,
    .lfg_eval  = &globl_lfg_eval,
    .lfg_jac  = &globl_lfg_jac,
    .lfg_lt_aug_hes  = &globl_lfg_lt_aug_hes,
    .mr_eval = &globl_mr_eval,
    .mr_jac = &globl_mr_jac,
    .mr_lt_aug_hes = &globl_mr_lt_aug_hes,
};

c_problem_t* get_update_c_problem() {
    return &globl_c_problem;
}

// [L, f, g]
void eval_lfg(const f64* xu, const f64* p, f64* out) {
    out[0] /* L */ = 0.5 * (xu[0] * xu[0] /* x^2 */ + xu[1] * xu[1] /* u */);
    out[1] /* L */ = xu[0] /* x */ + xu[1] /* u */;
}

// ∇ [L, f, g]
void jac_lfg(const f64* xu, const f64* p, f64* out) {
    out[0] = xu[0]; /* L_x = x */
    out[1] = xu[1]; /* L_u = u */
    out[2] = 1; /* f_x = 1 */
    out[3] = 1; /* f_u = 1 */
}

// σ ∇² L + λ^T ∇² [f, g] (lower triangle)
void hes_lfg(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64* out) {
    out[0] = obj_factor; /* ℒ_xx = 1 * sigma */
    out[1] = obj_factor; /* ℒ_uu = 1 * sigma */
}

// [M, r]
void eval_mr(const f64* x0, const f64* xuf, const f64* p, f64* out) {

}

// ∇ [M, r]
void jac_mr(const f64* x0, const f64* xuf, const f64* p, f64* out) {

}

// σ ∇² M + λ^T ∇² r (lower triangle)
void hes_mr(const f64* x0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64* out) {

}
