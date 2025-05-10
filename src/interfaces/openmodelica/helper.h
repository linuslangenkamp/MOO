#ifndef OPT_OM_HELPER
#define OPT_OM_HELPER

#include <base/util.h>

/* just simple helper structs */

struct InfoGDOP {
    int x_size;
    int u_size;
    int p_size;

    int xu_size;

    int f_size;
    int g_size;
    int r_size;

    int real_vars_g_start_index;
    int real_vars_r_start_index;

    bool mayer_exists;
    bool lagrange_exists;

    /* TODO: better save the index in real_vars => parallel computations */
    F64* address_mayer_real_vars;
    F64* address_lagrange_real_vars;
};

#endif // OPT_OM_HELPER
