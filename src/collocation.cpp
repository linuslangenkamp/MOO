#include "collocation.h"
#include "linalg.h"
#include "util.h"

double Collocation::integrate(const double* values, const size_t scheme) {
    /*
    input: values - f(c_1), ..., f(c_m)
    output: int_{0}^{1} f(t) dt \approx sum_{k=1}^{m} b_k * f(c_k), m steps, b_k weights, c_k nodes
    */

    return Linalg::dot(scheme, values, b[scheme].data());
};