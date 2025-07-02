#include "linalg.h"

namespace Linalg {

    /** 
     * @brief Perform dot product: ret = transpose(x) * y
     *
     * @param size      Vector length
     * @param x         First Vector
     * @param y         Second Vector
     * @return f64   transpose(vec1) * vec2
     */
    f64 dot(const int size, const f64* x, const f64* y) {
        f64 ret = 0;
        for (int i = 0; i < size; i++) {
            ret += x[i] * y[i];
        }
        return ret;
    }

    void Dv(const f64* D, bool invD, const f64* x, const int size, f64* out) {
        int i;

        if (invD) {
            for (i = 0; i < size; i++) {
                out[i] = x[i] / D[i];
            }
        }
        else {
            for (i = 0; i < size; i++) {
                out[i] = D[i] * x[i];
            }
        }
    }

    void Dv_inplace(const f64* D, bool invD, f64* x, const int size) {
        int i;

        if (invD) {
            for (i = 0; i < size; i++) {
                x[i] /= D[i];
            }
        }
        else {
            for (i = 0; i < size; i++) {
                x[i] *= D[i];
            }
        }
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
    void dsaxpy(const int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out) {
        int i;

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
    void dgmv(const int size, const f64* x, const f64* y, const f64* D, f64 beta, bool invD, f64* out) {
        int i;

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