#ifndef OPT_INTEGRATOR_H
#define OPT_INTEGRATOR_H

#include <vector>

#include "util.h"

struct Integrator {

    // nodes, i.e. [c1, c2, ..., cm]
    const std::vector<gVector> c =
#include "constantsC.data"

    // nodes including -1, i.e. [-1, c1, c2, ..., cm]
    const std::vector<gVector> c0 =
#include "constantsC0.data"

    // quadrature weights {}, {2}, ...
    const std::vector<gVector> b =
#include "constantsB.data"

    // differentiation matrices
    const std::vector<std::vector<gVector>> D =
#include "constantsD.data"

    // integral -1 to 1 of some values, scheme is chosen automatically
    gNumber integrate(const gVector& values);
};

#endif  // OPT_INTEGRATOR_H
