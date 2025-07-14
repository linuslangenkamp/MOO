#include "scaling.h"

namespace OpenModelica {

NominalScaling create_gdop_nominal_scaling(GDOP::GDOP& gdop, InfoGDOP& info) {
    // x, g, f of the NLP { min f(x) s.t. g_l <= g(x) <= g_l }
    auto x_nominal = FixedVector<f64>(gdop.number_vars);
    auto g_nominal = FixedVector<f64>(gdop.number_constraints);
    f64  f_nominal = 1;

    // get problem sizes
    auto x_size  = info.x_size;
    auto u_size  = info.u_size;
    auto xu_size = info.xu_size;
    auto f_size = info.f_size;
    auto g_size = info.g_size;
    auto r_size = info.r_size;
    auto fg_size = f_size + g_size;

    auto real_vars_data = info.data->modelData->realVarsData;

    auto has_mayer = gdop.problem.boundary->has_mayer;
    auto has_lagrange = gdop.problem.full->has_lagrange;

    if (has_mayer && has_lagrange) {
        f_nominal = (real_vars_data[info.index_mayer_real_vars].attribute.nominal +
                     real_vars_data[info.index_lagrange_real_vars].attribute.nominal) / 2;
    }
    else if (has_lagrange) {
        f_nominal = real_vars_data[info.index_lagrange_real_vars].attribute.nominal;
    }
    else if (has_mayer) {
        f_nominal = real_vars_data[info.index_mayer_real_vars].attribute.nominal;
    }
    else {
        // use default 1 - no objective set!
    }

    // x(t_0)
    for (int x = 0; x < info.x_size; x++) {
        x_nominal[x] = real_vars_data[x].attribute.nominal;
    }

    // (x, u)_(t_node)
    for (int node = 0; node < gdop.mesh.node_count; node++) {
        for (int x = 0; x < x_size; x++) {
            x_nominal[x_size + node * xu_size + x] = real_vars_data[x].attribute.nominal;
        }

        for (int u = 0; u < u_size; u++) {
            int u_real_vars = info.u_indices_real_vars[u];
            x_nominal[2 * x_size + node * xu_size + u] = real_vars_data[u_real_vars].attribute.nominal;
        }
    }

    for (int node = 0; node < gdop.mesh.node_count; node++) {
        for (int f = 0; f < f_size; f++) {
            g_nominal[node * fg_size + f] = x_nominal[f]; // reuse x nominal for dynamic for now!
        }

        for (int g = 0; g < g_size; g++) {
            g_nominal[f_size + node * fg_size + g] = real_vars_data[info.index_g_real_vars + g].attribute.nominal;
        }
    }

    for (int r = 0; r < r_size; r++) {
        g_nominal[gdop.off_fg_total + r] = real_vars_data[info.index_r_real_vars + r].attribute.nominal;
    }

    return NominalScaling(std::move(x_nominal), std::move(g_nominal), f_nominal);
}

} // namespace OpenModelica
