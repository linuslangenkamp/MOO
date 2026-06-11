// SPDX-License-Identifier: LGPL-3.0-or-later
#include "matrix.h"

#include <stdexcept>
#include <utility>

namespace ad {

DenseMatrix::DenseMatrix(int rows_in, int cols_in, std::vector<double> values_in)
    : rows(rows_in),
      cols(cols_in),
      values(std::move(values_in)) {
    if (rows < 0 || cols < 0) {
        throw std::runtime_error("dense matrix dimensions must be non-negative");
    }
    if (static_cast<int>(values.size()) != rows * cols) {
        throw std::runtime_error("dense matrix value count does not match dimensions");
    }
}

double DenseMatrix::operator()(int row, int col) const {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::out_of_range("dense matrix index out of range");
    }
    return values[static_cast<std::size_t>(row * cols + col)];
}

SparseMatrix::SparseMatrix(int rows_in, int cols_in, std::vector<int> row_in, std::vector<int> col_in, std::vector<double> values_in)
    : rows(rows_in),
      cols(cols_in),
      row(std::move(row_in)),
      col(std::move(col_in)),
      values(std::move(values_in)) {
    if (rows < 0 || cols < 0) {
        throw std::runtime_error("sparse matrix dimensions must be non-negative");
    }
    if (row.size() != col.size() || row.size() != values.size()) {
        throw std::runtime_error("sparse matrix triplet arrays must have equal length");
    }
    for (int k = 0; k < nnz(); ++k) {
        const int r = row[static_cast<std::size_t>(k)];
        const int c = col[static_cast<std::size_t>(k)];
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::runtime_error("sparse matrix entry out of bounds");
        }
    }
}

int SparseMatrix::nnz() const {
    return static_cast<int>(values.size());
}

} // namespace ad
