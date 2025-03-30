#ifndef OPT_UTIL_H
#define OPT_UTIL_H

#include <vector>


template <typename T>
inline int int_size(const std::vector<T>& vec) {
    return static_cast<int>(vec.size());
}

#endif // OPT_UTIL_H
