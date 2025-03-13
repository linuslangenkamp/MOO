#ifndef OPT_LINALG_H
#define OPT_LINALG_H

#include <cstddef> 

// these are possibly blas functions that can be replaced with proper BLAS later
// for now its sufficient to memorize where these subroutines are used

namespace Linalg {

    /** 
     * @brief Perform dot product: ret = transpose(x) * y
     *
     * @param size      Vector length
     * @param x         First Vector
     * @param y         Second Vector
     * @return double   transpose(vec1) * vec2
     */
    double dot(const size_t size, const double* x, const double* y) {
        size_t i;
        double ret = 0;

        for (i = 0; i < size; i++) {
            ret += x[i] * y[i];
        }
        return ret;
    }

    // these are not needed rn, scaling is delayed
    /** 
     * @brief Perform diagonal scaled axpy: out = D * (x + beta * y) or out = D^(-1) * (x + beta * y)
     *
     * @param size      Vector length
     * @param x         Vector x
     * @param y         Vector y
     * @param D         Diagonal Matrix D
     * @param beta      Scaling factor of beta
     * @param invD      Invert Matrix D
     * @param out       Vector to fill
     */
    void dsaxpy(const size_t size, const double* x, const double* y, const double* D, double beta, bool invD, double* out) {
        size_t i;

        if (invD) {
            for (i = 0; i < size; i++) {
                out[i] = (x[i] + beta * y[i]) / D[i];
            }
        }
        else {
            for (i = 0; i < size; i++) {
                out[i] = D[i] * (x[i] + beta * y[i]);
            }
        }
    }
    /** 
     * @brief Perform diagonal general matrix vector: out = D * x + beta * y or D^(-1) * x + beta * y
     *
     * @param size      Vector length
     * @param x         Vector x
     * @param y         Vector y
     * @param D         Diagonal Matrix D
     * @param beta      Scaling factor of beta
     * @param invD      Invert Matrix D
     * @param out       Vector to fill
     */
    void dgmv(const size_t size, const double* x, const double* y, const double* D, double beta, bool invD, double* out) {
        size_t i;

        if (invD) {
            for (i = 0; i < size; i++) {
                out[i] = x[i] / D[i] + beta * y[i];
            }
        }
        else {
            for (i = 0; i < size; i++) {
                out[i] =  D[i] * x[i] + beta * y[i];
            }
        }
    }

}

#endif  // OPT_LINALG_H
