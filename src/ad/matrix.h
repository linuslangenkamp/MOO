// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_MATRIX_H
#define MOO_AD_MATRIX_H

#include <vector>

namespace ad {

struct DenseMatrix {
    int rows = 0;
    int cols = 0;
    std::vector<double> values;

    DenseMatrix() = default;
    DenseMatrix(int rows, int cols, std::vector<double> values);

    double operator()(int row, int col) const;
};

struct SparseMatrix {
    int rows = 0;
    int cols = 0;
    std::vector<int> row;
    std::vector<int> col;
    std::vector<double> values;

    SparseMatrix() = default;
    SparseMatrix(int rows, int cols, std::vector<int> row, std::vector<int> col, std::vector<double> values);

    int nnz() const;
};

} // namespace ad

#endif // MOO_AD_MATRIX_H
