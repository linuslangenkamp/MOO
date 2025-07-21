#include "linalg.h"

namespace Linalg {

/** 
 * @brief Perform dot product: ret = transpose(x) * y
 *
 * @param size      Vector length
 * @param x         First Vector ptr
 * @param y         Second Vector ptr
 * @return f64   transpose(vec1) * vec2
 */
f64 dot(const int size, const f64* x, const f64* y) {
    f64 ret = 0.0;
    for (int i = 0; i < size; i++) {
        ret += x[i] * y[i];
    }
    return ret;
}

/** 
 * @brief Perform inplace squaring (x_1^2, x_2^2, ..., x_size^2)
 *
 * @param size      Vector length
 * @param x         Vector ptr
 */
void square(const int size, f64* x) {
    for (int i = 0; i < size; i++) {
        x[i] *= x[i];
    }
}


void diagmat_vec(const f64* D, bool invD, const f64* x, const int size, f64* out) {
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

void diagmat_vec_inplace(const f64* D, bool invD, f64* x, const int size) {
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

} // namespace Linalg
