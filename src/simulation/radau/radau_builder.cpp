#include <simulation/radau/radau_builder.h>

namespace Simulation {

RadauBuilder& RadauBuilder::radau_scheme(RadauScheme radau_scheme_) {
    scheme = radau_scheme_;
    return *this;
}

RadauBuilder& RadauBuilder::radau_h0(f64 h_init_) {
    h_init = h_init_;
    return *this;
}

RadauBuilder& RadauBuilder::radau_tol(f64 atol_, f64 rtol_) {
    atol = atol_;
    rtol = rtol_;
    return *this;
}

RadauIntegrator RadauBuilder::build() const {
    f64 h0 = h_init != 0.0 ? h_init : (dense_output_grid.back() - dense_output_grid[0]) / (2 * dense_output_grid.size());

    return RadauIntegrator(
            /* every Integrator */
            ode_func, dense_output_grid,
            x_start_values, x_size,
            user_data, parameters, p_size,
            controls, data,
            jac_func, jac_fmt,
            i_row, j_col, nnz,
            /* Radau specific */
            scheme,
            h0,
            atol,
            rtol
        );
}

} // namespace Simulation
