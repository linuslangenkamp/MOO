#ifndef MOO_C_STRUCTS_H
#define MOO_C_STRUCTS_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdbool.h>

typedef double f64;

typedef struct bounds_t {
    f64 lb;
    f64 ub;
} bounds_t;

typedef struct optional_value_t {
    f64 value;
    bool is_set;
} optional_value_t;

// eval_structure_t[fn] -> buf_index of out in eval_lfg() - fn is sorted as L -> f -> g
typedef struct eval_structure_t {
    int* buf_index;  // buf_index
} eval_structure_t;

typedef struct coo_t {
    int* row;        // row indices
    int* col;        // col indices
    int* buf_index;  // non-zero index == buf_index
    int nnz;         // total nnz
} coo_t;

typedef struct c_problem_t {
    const int x_size;
    const int u_size;
    const int xu_size;
    const int p_size;
    const int rp_size;

    const int r_size;
    const int g_size;

    const bool has_mayer;
    const bool has_lagrange;

    bounds_t* x_bounds;
    bounds_t* u_bounds;
    bounds_t* p_bounds;

    bounds_t* r_bounds;
    bounds_t* g_bounds;

    optional_value_t* x0_fixed;
    optional_value_t* xf_fixed;

    f64* x_nominal;
    f64* u_nominal;
    f64* p_nominal;

    f64* obj_nominal;
    f64* f_nominal;
    f64* g_nominal;
    f64* r_nominal;

    eval_structure_t* lfg_eval;
    coo_t* lfg_jac;
    coo_t* lfg_lt_aug_hes;

    eval_structure_t* mr_eval;
    coo_t* mr_jac;
    coo_t* mr_lt_aug_hes;
} c_problem_t;

// pass these in some callback struct
c_problem_t* get_update_c_problem();

void eval_lfg(const f64* xu, const f64* p, f64 t, f64* out);
void jac_lfg(const f64* xu, const f64* p, f64 t, f64* out);
void hes_lfg(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64 t, f64* out);
void eval_mr(const f64* x0, const f64* xuf, const f64* p, f64 t0, f64 tf, f64* out);
void jac_mr(const f64* x0, const f64* xuf, const f64* p, f64 t0, f64 tf, f64* out);
void hes_mr(const f64* x0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64 t0, f64 tf, f64* out);

#ifdef __cplusplus
}
#endif

#endif // MOO_C_STRUCTS_H
