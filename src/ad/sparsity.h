// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_SPARSITY_H
#define MOO_AD_SPARSITY_H

#include "symbol.h"

#include <utility>
#include <vector>

namespace ad {

class Expr;
class Function;
class Vec;

class SparsityPattern {
public:
    SparsityPattern() = default;
    SparsityPattern(int rows, int cols, std::vector<std::pair<int, int>> entries);

    int rows() const;
    int cols() const;
    const std::vector<std::pair<int, int>> &entries() const;
    bool contains(int row, int col) const;
    bool empty() const;
    int nnz() const;

private:
    int rows_ = 0;
    int cols_ = 0;
    std::vector<std::pair<int, int>> entries_;
};

SparsityPattern sparsity(const Expr &expr, const Vars &wrt);
SparsityPattern sparsity(const Vec &vec, const Vars &wrt);

} // namespace ad

#endif // MOO_AD_SPARSITY_H
