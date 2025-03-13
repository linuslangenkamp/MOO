#include "integrator.h"
#include "util.h"

gNumber Integrator::integrate(const gVector& values) {
    /*
    input: values - f(c_1), ..., f(c_m)
    output: int_{-1}^{1} f(t) dt \approx sum_{k=1}^{m} b_k * f(c_k), m steps, b_k weights, c_k nodes
    */

    gNumber integral = 0;
    for (int k = 0; k < values.size(); k++) {
        integral += b[values.size()][k] * values[k];
    }

    return integral;
};