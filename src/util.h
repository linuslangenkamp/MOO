#ifndef OPT_UTIL_H
#define OPT_UTIL_H

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <sstream>
#include <span>
#include <tuple>
#include <optional>
#include <vector>

const double PLUS_INFINITY = std::numeric_limits<double>::infinity();
const double MINUS_INFINITY = -std::numeric_limits<double>::infinity();

namespace Util {

    template <typename T>
    std::string vectorToString(const std::vector<T>& vec) {
        std::stringstream out;
        out << std::fixed << std::setprecision(15);
        if (!vec.empty()) {
            out << vec[0];
            for (int i = 1; i < vec.size(); i++) {
                out << ", " << vec[i];
            }
        }
        return out.str();
    }
}

#endif  // OPT_UTIL_H
