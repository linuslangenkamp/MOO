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
#include <set>
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

    struct OrderedIndexSet {
        struct Compare {
            bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
                if (a.first != b.first) {
                    return a.first < b.first;
                } else {
                    return a.second < b.second;
                }
            }
        };

        std::set<std::pair<int, int>, Compare> set;

        void insertSparsity(std::vector<HessianSparsity>& hes, int row_off, int col_off) {
            for (auto& coo : hes) {
                set.insert({coo.index1 + row_off, coo.index2 + col_off});
            }
        }

        inline size_t size() const {
            return set.size();
        }
    };
}

#endif  // OPT_UTIL_H
