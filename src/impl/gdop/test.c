#include "test.h"

double x;

void cosine(const double* val) {
    x = 1 - val[0]*val[0] / 2 + val[1];
}