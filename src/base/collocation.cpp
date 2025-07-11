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

/**
 * @brief Interpolates a value at a given point T using barycentric interpolation. @runtime: O(scheme).
 *
 * This method performs barycentric interpolation based on a set of collocation nodes
 * and corresponding barycentric weights. It first rescales the target point T
 * to the nominal interval domain of the nodes and then computes the interpolated
 * value. It includes a check for exact matches with existing nodes to prevent
 * division by zero.
 *
 * @param scheme         The number of collocation points (nodes and weights).
 * @param contains_zero  A boolean indicating whether the `c0`/`w0` (containing zero)
 *                       or `c`/`w` (not containing zero) sets of nodes/weights should be used.
 * @param values         A pointer to an array of function values corresponding to the collocation nodes.
 * @param increment      The stride to use when accessing elements in the `values` array.
 * @param interval_start The start of the physical interval.
 * @param interval_end   The end of the physical interval.
 * @param point          The point at which to interpolate.
 * @return The interpolated value at point.
 */
f64 Collocation::interpolate(int scheme, bool contains_zero, const f64* values, int increment,
                             f64 interval_start, f64 interval_end, f64 point) const {
    const auto& nodes   = contains_zero ? c0[scheme] : c[scheme];
    const auto& weights = contains_zero ? w0[scheme] : w[scheme];

    // rescale T to the [0, 1] nominal interval domain
    f64 h = interval_end - interval_start;
    f64 node_start = nodes[0];
    f64 node_end = nodes[scheme - 1];
    f64 point_hat = (point - interval_start) / h * (node_end - node_start) + node_start;

    // check for exact match with any node to avoid division by zero
    for (int j = 0; j < scheme; j++) {
        if (std::abs(point_hat - nodes[j]) < 1e-14) return values[j];
    }

    // compute the barycentric interpolant
    f64 numerator = 0.0;
    f64 denominator = 0.0;
    for (int j = 0; j < scheme; j++) {
        f64 temp = weights[j] / (point_hat - nodes[j]);
        numerator += temp * values[increment * j];
        denominator += temp;
    }

    return numerator / denominator;
};
