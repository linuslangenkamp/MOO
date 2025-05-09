#ifndef OPT_COLLOCATION_H
#define OPT_COLLOCATION_H

#include <vector>

#include "util.h"
#include "linalg.h"

struct Collocation {

    // nodes, i.e. [c1, c2, ..., cm]
    const std::vector<std::vector<F64>> c =
#include "../data/radauConstantsC.data"

    // nodes including 0, i.e. [0, c1, c2, ..., cm]
    const std::vector<std::vector<F64>> c0 =
#include "../data/radauConstantsC0.data"

    // quadrature weights {}, {1}, ...
    const std::vector<std::vector<F64>> b =
#include "../data/radauConstantsB.data"

    // differentiation matrices
    const std::vector<std::vector<std::vector<F64>>> D =
#include "../data/radauConstantsD.data"

    // integral 0 to 1 of some values
    F64 integrate(const F64* values, const int scheme);

    // multiply given differentiation matrix scheme with x_prev, (x_i0, ui0, xi1, ui1, ..., xnm, unm) u only for offset
    void diff_matrix_multiply(const int scheme, const int x_size, const int xu_size, const int fg_size, const F64* x_prev, const F64* x_new, F64* out);
};

#endif  // OPT_COLLOCATION_H
