#include "collocation.h"

/*
 * input: values - f(c_1), ..., f(c_m)
 * output: int_{0}^{1} f(t) dt \approx sum_{k=1}^{m} b_k * f(c_k), m stages, b_k weights, c_k nodes
 */
f64 Collocation::integrate(const f64* values, const int scheme) {
    return Linalg::dot(scheme, values, b[scheme].data());
};

void Collocation::diff_matrix_multiply(const int scheme, const int x_size, const int xu_size, const int fg_size,
                                       const f64* x_prev, const f64* x_new, f64* out) {
    for (int row = 1; row < scheme + 1; row++) {
        for (int x_index = 0; x_index < x_size; x_index++) {
            int out_row = (row - 1) * fg_size + x_index;

            // col = 0 -> x_pref
            out[out_row] += D[scheme][row][0] * x_prev[x_index];

            // col > 0 -> x_new
            for (int col = 1; col < scheme + 1; col++) {
                out[out_row] += D[scheme][row][col] * x_new[(col - 1) * xu_size + x_index];
            }
        }
    }
};

f64 Collocation::interpolate(int scheme, bool contains_zero, const f64* values, int increment,
                             f64 interval_start, f64 interval_end, f64 T) const {
    const auto& nodes   = contains_zero ? c0[scheme] : c[scheme];
    const auto& weights = contains_zero ? w0[scheme] : w[scheme];

    // rescale T to the [0, 1] nominal interval domain
    f64 h = interval_end - interval_start;
    f64 node_start = nodes[0];
    f64 node_end = nodes[scheme - 1];
    f64 T_hat = (T - interval_start) / h * (node_end - node_start) + node_start;

    // check for exact match with any node to avoid division by zero
    for (int j = 0; j < scheme; j++) {
        if (std::abs(T_hat - nodes[j]) < 1e-14) return values[j];
    }

    // compute the barycentric interpolant
    f64 numerator = 0.0;
    f64 denominator = 0.0;
    for (int j = 0; j < scheme; j++) {
        f64 temp = weights[j] / (T_hat - nodes[j]);
        numerator += temp * values[increment * j];
        denominator += temp;
    }

    return numerator / denominator;
};
