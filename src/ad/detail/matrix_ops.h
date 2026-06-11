// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_DETAIL_MATRIX_OPS_H
#define MOO_AD_DETAIL_MATRIX_OPS_H

#include "../matrix.h"

#include <limits>
#include <stdexcept>
#include <utility>

namespace ad::detail {

inline int checked_matrix_size(int rows, int cols) {
    if (rows < 0 || cols < 0) {
        throw std::runtime_error("matrix dimensions must be non-negative");
    }
    const long long product = static_cast<long long>(rows) * static_cast<long long>(cols);
    if (product > std::numeric_limits<int>::max()) {
        throw std::runtime_error("matrix size overflows int");
    }
    return static_cast<int>(product);
}

inline int matrix_flat_index(int row, int col, int rows, int cols, MatrixLayout layout) {
    if (row < 0 || row >= rows || col < 0 || col >= cols) {
        throw std::runtime_error("matrix index out of bounds");
    }
    if (layout == MatrixLayout::ColumnMajor) {
        return row + col * rows;
    }
    return row * cols + col;
}

inline std::pair<int, int> matrix_coordinates_from_flat(int index, int rows, int cols, MatrixLayout layout) {
    if (index < 0 || index >= checked_matrix_size(rows, cols)) {
        throw std::runtime_error("matrix flat index out of bounds");
    }
    if (layout == MatrixLayout::ColumnMajor) {
        return {index % rows, index / rows};
    }
    return {index / cols, index % cols};
}

} // namespace ad::detail

#endif // MOO_AD_DETAIL_MATRIX_OPS_H
