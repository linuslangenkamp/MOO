#include "collocation.h"


double Collocation::integrate(const double* values, const int scheme) {
    /*
    input: values - f(c_1), ..., f(c_m)
    output: int_{0}^{1} f(t) dt \approx sum_{k=1}^{m} b_k * f(c_k), m steps, b_k weights, c_k nodes
    */

    return Linalg::dot(scheme, values, b[scheme].data());
};

void Collocation::diff_matrix_multiply(const int scheme, const int x_size, const int xu_size, const int fg_size,
                                       const double* x_prev, const double* x_new, double* out) {
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
