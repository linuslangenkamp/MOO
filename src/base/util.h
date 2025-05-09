#ifndef OPT_UTIL_H
#define OPT_UTIL_H

#include <vector>
#include <cmath>

typedef double F64;

const F64 PLUS_INFINITY = std::numeric_limits<F64>::infinity();
const F64 MINUS_INFINITY = -std::numeric_limits<F64>::infinity();

template <typename T>
inline int int_size(const std::vector<T>& vec) {
    return static_cast<int>(vec.size());
}

#endif // OPT_UTIL_H
