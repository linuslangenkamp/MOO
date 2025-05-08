#ifndef OPT_LINALG_H
#define OPT_LINALG_H

#include "util.h"

// these are possibly blas functions that can be replaced with proper BLAS later
// for now its sufficient to memorize where these subroutines are used

namespace Linalg {
    f64 dot(const int size, const f64* x, const f64* y);
    void dsaxpy(const int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out);
    void dgmv(const int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out);
}

#endif  // OPT_LINALG_H
