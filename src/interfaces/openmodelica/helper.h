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

    bool mayer_exists;
    bool lagrange_exists;

    /* TODO: better save the index in real_vars => parallel computations */
    modelica_real* __address_mayer_real_vars;
    modelica_real* __address_lagrange_real_vars;

    const int index_x_real_vars = 0;
    int index_der_x_real_vars;
    int index_mayer_real_vars = -1;
    int index_lagrange_real_vars = -1;
    FixedVector<int> u_indices_real_vars;
    int index_g_real_vars;
    int index_r_real_vars;
};

#endif // OPT_OM_HELPER
