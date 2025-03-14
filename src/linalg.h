#ifndef OPT_LINALG_H
#define OPT_LINALG_H

#include <cstddef> 

// these are possibly blas functions that can be replaced with proper BLAS later
// for now its sufficient to memorize where these subroutines are used

namespace Linalg {
    double dot(const size_t size, const double* x, const double* y);
    void dsaxpy(const size_t size, const double* x, const double* y, const double* D, double beta, bool invD, double* out);
    void dgmv(const size_t size, const double* x, const double* y, const double* D, double beta, bool invD, double* out);
}

#endif  // OPT_LINALG_H
