#ifndef OPT_LINALG_H
#define OPT_LINALG_H

#include "util.h"

// these are possibly blas functions that can be replaced with proper BLAS later
// for now its sufficient to memorize where these subroutines are used

namespace Linalg {
    F64 dot(const int size, const F64* x, const F64* y);
    void dsaxpy(const int size, const F64* x, const F64* y, const F64* D, F64 beta, bool invD, F64* out);
    void dgmv(const int size, const F64* x, const F64* y, const F64* D, F64 beta, bool invD, F64* out);
}

#endif  // OPT_LINALG_H
