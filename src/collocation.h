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
    double integrate(const double* values, const size_t scheme);
};

#endif  // OPT_COLLOCATION_H
