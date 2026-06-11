// SPDX-License-Identifier: LGPL-3.0-or-later
#include "sparse_mat.h"

#include "detail/matrix_ops.h"
#include "detail/simplify.h"

#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace ad {
namespace {

struct CanonicalPattern {
    std::vector<int> row;
    std::vector<int> col;
    std::vector<int> order;
    bool identity_order = true;
};

CanonicalPattern canonicalize_pattern(int rows, int cols, const std::vector<int> &row, const std::vector<int> &col) {
    if (rows < 0 || cols < 0) {
        throw std::runtime_error("sparse matrix dimensions must be non-negative");
    }
    if (row.size() != col.size()) {
        throw std::runtime_error("sparse matrix row and column arrays must have equal length");
    }

    std::vector<int> order(row.size());
    for (std::size_t i = 0; i < order.size(); ++i) {
        const int r = row[i];
        const int c = col[i];
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            throw std::runtime_error("sparse matrix entry out of bounds");
        }
        order[i] = static_cast<int>(i);
    }

    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return std::tie(row[static_cast<std::size_t>(lhs)], col[static_cast<std::size_t>(lhs)], lhs) <
               std::tie(row[static_cast<std::size_t>(rhs)], col[static_cast<std::size_t>(rhs)], rhs);
    });

    CanonicalPattern result;
    result.row.reserve(row.size());
    result.col.reserve(col.size());
    result.order = order;
    for (std::size_t out = 0; out < order.size(); ++out) {
        const int source = order[out];
        if (source != static_cast<int>(out)) {
            result.identity_order = false;
        }
        const int r = row[static_cast<std::size_t>(source)];
        const int c = col[static_cast<std::size_t>(source)];
        if (!result.row.empty() && result.row.back() == r && result.col.back() == c) {
            throw std::runtime_error("sparse matrix duplicate structural entry");
        }
        result.row.push_back(r);
        result.col.push_back(c);
    }
    return result;
}

void require_valid(const SparseMat &matrix) {
    if (!matrix.valid()) {
        throw std::runtime_error("invalid sparse matrix expression");
    }
}

void require_same_shape(const SparseMat &lhs, const SparseMat &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
        throw std::runtime_error("sparse matrix shape mismatch");
    }
}

bool less_coord(int lhs_row, int lhs_col, int rhs_row, int rhs_col) {
    return lhs_row < rhs_row || (lhs_row == rhs_row && lhs_col < rhs_col);
}

SparseMat combine_union(const SparseMat &lhs, const SparseMat &rhs, bool subtract_rhs) {
    require_same_shape(lhs, rhs);

    std::vector<int> row;
    std::vector<int> col;
    std::vector<Expr> elements;
    int i = 0;
    int j = 0;
    while (i < lhs.nnz() || j < rhs.nnz()) {
        if (j == rhs.nnz() ||
            (i < lhs.nnz() && less_coord(lhs.row_indices()[static_cast<std::size_t>(i)],
                                         lhs.col_indices()[static_cast<std::size_t>(i)],
                                         rhs.row_indices()[static_cast<std::size_t>(j)],
                                         rhs.col_indices()[static_cast<std::size_t>(j)]))) {
            row.push_back(lhs.row_indices()[static_cast<std::size_t>(i)]);
            col.push_back(lhs.col_indices()[static_cast<std::size_t>(i)]);
            elements.push_back(lhs.values()[i]);
            ++i;
        } else if (i == lhs.nnz() ||
                   less_coord(rhs.row_indices()[static_cast<std::size_t>(j)],
                              rhs.col_indices()[static_cast<std::size_t>(j)],
                              lhs.row_indices()[static_cast<std::size_t>(i)],
                              lhs.col_indices()[static_cast<std::size_t>(i)])) {
            row.push_back(rhs.row_indices()[static_cast<std::size_t>(j)]);
            col.push_back(rhs.col_indices()[static_cast<std::size_t>(j)]);
            elements.push_back(subtract_rhs ? -rhs.values()[j] : rhs.values()[j]);
            ++j;
        } else {
            row.push_back(lhs.row_indices()[static_cast<std::size_t>(i)]);
            col.push_back(lhs.col_indices()[static_cast<std::size_t>(i)]);
            elements.push_back(subtract_rhs ? lhs.values()[i] - rhs.values()[j]
                                            : lhs.values()[i] + rhs.values()[j]);
            ++i;
            ++j;
        }
    }
    return SparseMat(vec(std::move(elements)), lhs.rows(), lhs.cols(), std::move(row), std::move(col));
}

} // namespace

SparseMat::SparseMat(Vec values, int rows, int cols, std::vector<int> row, std::vector<int> col)
    : values_(std::move(values)),
      rows_(rows),
      cols_(cols) {
    if (!values_.valid()) {
        throw std::runtime_error("sparse matrix values must be a valid Vec");
    }
    CanonicalPattern pattern = canonicalize_pattern(rows_, cols_, row, col);
    if (values_.size() != static_cast<int>(pattern.row.size())) {
        throw std::runtime_error("sparse matrix value count must match structural nonzero count");
    }
    if (!pattern.identity_order) {
        values_ = gather(values_, pattern.order);
    }
    row_ = std::move(pattern.row);
    col_ = std::move(pattern.col);
}

bool SparseMat::valid() const {
    return values_.valid();
}

int SparseMat::rows() const {
    return rows_;
}

int SparseMat::cols() const {
    return cols_;
}

int SparseMat::nnz() const {
    return static_cast<int>(row_.size());
}

const Vec &SparseMat::values() const {
    return values_;
}

const std::vector<int> &SparseMat::row_indices() const {
    return row_;
}

const std::vector<int> &SparseMat::col_indices() const {
    return col_;
}

Expr SparseMat::operator()(int row, int col) const {
    require_valid(*this);
    if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
        throw std::runtime_error("sparse matrix index out of bounds");
    }
    const auto found = std::lower_bound(row_.begin(), row_.end(), row);
    for (auto it = found; it != row_.end() && *it == row; ++it) {
        const std::size_t index = static_cast<std::size_t>(std::distance(row_.begin(), it));
        if (col_[index] == col) {
            return values_[static_cast<int>(index)];
        }
    }
    return constant(0.0);
}

SparseMat SparseMat::transpose() const {
    require_valid(*this);
    return SparseMat(values_, cols_, rows_, col_, row_);
}

Mat SparseMat::to_dense(MatrixLayout layout) const {
    require_valid(*this);
    std::vector<int> indices;
    indices.reserve(row_.size());
    for (std::size_t k = 0; k < row_.size(); ++k) {
        indices.push_back(detail::matrix_flat_index(row_[k], col_[k], rows_, cols_, layout));
    }
    return Mat(scatter_add(values_, std::move(indices), detail::checked_matrix_size(rows_, cols_)),
               rows_,
               cols_,
               layout);
}

void SparseMat::eval(const Values &values, EvalWorkspace &workspace, double *out, int out_size) const {
    require_valid(*this);
    values_.eval(values, workspace, out, out_size);
}

SparseMat sparse_mat(Vec values, int rows, int cols, std::vector<int> row, std::vector<int> col) {
    return SparseMat(std::move(values), rows, cols, std::move(row), std::move(col));
}

SparseMat sparse_mat_variable(const std::string &label, int rows, int cols, std::vector<int> row, std::vector<int> col) {
    CanonicalPattern pattern = canonicalize_pattern(rows, cols, row, col);
    const int nnz = static_cast<int>(pattern.row.size());
    return SparseMat(vec_variable(label, nnz), rows, cols, std::move(pattern.row), std::move(pattern.col));
}

SparseMat sparse_mat_parameter(const std::string &label, int rows, int cols, std::vector<int> row, std::vector<int> col) {
    CanonicalPattern pattern = canonicalize_pattern(rows, cols, row, col);
    const int nnz = static_cast<int>(pattern.row.size());
    return SparseMat(vec_parameter(label, nnz), rows, cols, std::move(pattern.row), std::move(pattern.col));
}

SparseMat sparse_mat_constant(int rows, int cols, std::vector<int> row, std::vector<int> col, std::vector<double> values) {
    return SparseMat(vec_constant(std::move(values)), rows, cols, std::move(row), std::move(col));
}

SparseMat operator+(const SparseMat &lhs, const SparseMat &rhs) {
    return combine_union(lhs, rhs, false);
}

SparseMat operator-(const SparseMat &lhs, const SparseMat &rhs) {
    return combine_union(lhs, rhs, true);
}

SparseMat operator*(const Expr &lhs, const SparseMat &rhs) {
    require_valid(rhs);
    return SparseMat(lhs * rhs.values(), rhs.rows(), rhs.cols(), rhs.row_indices(), rhs.col_indices());
}

SparseMat operator*(const SparseMat &lhs, const Expr &rhs) {
    require_valid(lhs);
    return SparseMat(lhs.values() * rhs, lhs.rows(), lhs.cols(), lhs.row_indices(), lhs.col_indices());
}

SparseMat operator*(double lhs, const SparseMat &rhs) {
    require_valid(rhs);
    return SparseMat(lhs * rhs.values(), rhs.rows(), rhs.cols(), rhs.row_indices(), rhs.col_indices());
}

SparseMat operator*(const SparseMat &lhs, double rhs) {
    require_valid(lhs);
    return SparseMat(lhs.values() * rhs, lhs.rows(), lhs.cols(), lhs.row_indices(), lhs.col_indices());
}

SparseMat operator-(const SparseMat &expr) {
    require_valid(expr);
    return SparseMat(-expr.values(), expr.rows(), expr.cols(), expr.row_indices(), expr.col_indices());
}

SparseMat hadamard(const SparseMat &lhs, const SparseMat &rhs) {
    require_same_shape(lhs, rhs);

    std::vector<int> row;
    std::vector<int> col;
    std::vector<Expr> elements;
    int i = 0;
    int j = 0;
    while (i < lhs.nnz() && j < rhs.nnz()) {
        if (less_coord(lhs.row_indices()[static_cast<std::size_t>(i)],
                       lhs.col_indices()[static_cast<std::size_t>(i)],
                       rhs.row_indices()[static_cast<std::size_t>(j)],
                       rhs.col_indices()[static_cast<std::size_t>(j)])) {
            ++i;
        } else if (less_coord(rhs.row_indices()[static_cast<std::size_t>(j)],
                              rhs.col_indices()[static_cast<std::size_t>(j)],
                              lhs.row_indices()[static_cast<std::size_t>(i)],
                              lhs.col_indices()[static_cast<std::size_t>(i)])) {
            ++j;
        } else {
            row.push_back(lhs.row_indices()[static_cast<std::size_t>(i)]);
            col.push_back(lhs.col_indices()[static_cast<std::size_t>(i)]);
            elements.push_back(lhs.values()[i] * rhs.values()[j]);
            ++i;
            ++j;
        }
    }
    return SparseMat(vec(std::move(elements)), lhs.rows(), lhs.cols(), std::move(row), std::move(col));
}

Vec operator*(const SparseMat &lhs, const Vec &rhs) {
    require_valid(lhs);
    return detail::make_symbolic_sparse_matvec(lhs.values(), lhs.rows(), lhs.cols(), lhs.row_indices(), lhs.col_indices(), rhs);
}

Mat operator*(const SparseMat &lhs, const Mat &rhs) {
    require_valid(lhs);
    if (!rhs.valid()) {
        throw std::runtime_error("invalid dense matrix expression");
    }
    return Mat(detail::make_symbolic_sparse_matmul(lhs.values(),
                                                   lhs.rows(),
                                                   lhs.cols(),
                                                   lhs.row_indices(),
                                                   lhs.col_indices(),
                                                   rhs.vec(),
                                                   rhs.rows(),
                                                   rhs.cols(),
                                                   rhs.layout(),
                                                   MatrixLayout::ColumnMajor,
                                                   true),
               lhs.rows(),
               rhs.cols(),
               MatrixLayout::ColumnMajor);
}

Mat operator*(const Mat &lhs, const SparseMat &rhs) {
    if (!lhs.valid()) {
        throw std::runtime_error("invalid dense matrix expression");
    }
    require_valid(rhs);
    return Mat(detail::make_symbolic_sparse_matmul(rhs.values(),
                                                   rhs.rows(),
                                                   rhs.cols(),
                                                   rhs.row_indices(),
                                                   rhs.col_indices(),
                                                   lhs.vec(),
                                                   lhs.rows(),
                                                   lhs.cols(),
                                                   lhs.layout(),
                                                   lhs.layout(),
                                                   false),
               lhs.rows(),
               rhs.cols(),
               lhs.layout());
}

} // namespace ad
