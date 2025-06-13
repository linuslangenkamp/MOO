#ifndef OPT_COLLOCATION_H
#define OPT_COLLOCATION_H

#include <vector>

#include "util.h"
#include "linalg.h"

struct Collocation {

    // nodes, i.e. [c1, c2, ..., cm]
    const std::vector<std::vector<f64>> c =
#include "../data/radauConstantsC.data"

    // nodes including 0, i.e. [0 = c0, c1, c2, ..., cm = 1]
    const std::vector<std::vector<f64>> c0 =
#include "../data/radauConstantsC0.data"

    // quadrature weights {}, {1}, ...
    const std::vector<std::vector<f64>> b =
#include "../data/radauConstantsB.data"

    // differentiation matrices
    const std::vector<std::vector<std::vector<f64>>> D =
#include "../data/radauConstantsD.data"

    // barycentric weights (only inner nodes [c1, ..., cm])
    const std::vector<std::vector<f64>> w0 =
#include "../data/radauConstantsW0.data"

    // barycentric weights (all nodes [0 = c0, c1, ..., cm = 1])
    const std::vector<std::vector<f64>> w =
#include "../data/radauConstantsW.data"

    // integral 0 to 1 of some values
    f64 integrate(const f64* values, const int scheme);

    // multiply given differentiation matrix scheme with x_prev, (x_i0, ui0, xi1, ui1, ..., xnm, unm) u only for offset
    void diff_matrix_multiply(const int scheme, const int x_size, const int xu_size, const int fg_size, const f64* x_prev, const f64* x_new, f64* out);

    // evaluate the interpolating polynomial at some point T
    f64 interpolate(int scheme, bool contains_zero, const f64* values, int stride, f64 interval_start, f64 interval_end, f64 T) const;
};

#endif  // OPT_COLLOCATION_H
