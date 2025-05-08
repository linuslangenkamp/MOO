#ifndef OPT_UTIL_H
#define OPT_UTIL_H

#include <vector>
#include <cmath>

typedef double f64;

const f64 PLUS_INFINITY = std::numeric_limits<f64>::infinity();
const f64 MINUS_INFINITY = -std::numeric_limits<f64>::infinity();

template <typename T>
inline int int_size(const std::vector<T>& vec) {
    return static_cast<int>(vec.size());
}

#endif // OPT_UTIL_H
