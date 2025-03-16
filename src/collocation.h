#ifndef OPT_COLLOCATION_H
#define OPT_COLLOCATION_H

#include <vector>

#include "util.h"

struct Collocation {

    // nodes, i.e. [c1, c2, ..., cm]
    const std::vector<std::vector<double>> c =
#include "../data/radauConstantsC.data"

    // nodes including -1, i.e. [-1, c1, c2, ..., cm]
    const std::vector<std::vector<double>> c0 =
#include "../data/radauConstantsC0.data"

    // quadrature weights {}, {2}, ...
    const std::vector<std::vector<double>> b =
#include "../data/radauConstantsB.data"

    // differentiation matrices
    const std::vector<std::vector<std::vector<double>>> D =
#include "../data/radauConstantsD.data"

    // integral 0 to 1 of some values
    double integrate(const double* values, const int scheme);

    // multiply given differentiation matrix scheme with x_prev, (x_i0, ui0, xi1, ui1, ..., xnm, unm) u only for offset
    void diff_matrix_multiply(const int scheme, const int x_size, const int xu_size, const int fg_size, const double* x_prev, const double* x_new, double* out);
};

#endif  // OPT_COLLOCATION_H
