// SPDX-License-Identifier: LGPL-3.0-or-later
#include "mat.h"

#include "detail/matrix_ops.h"
#include "detail/simplify.h"

#include <stdexcept>
#include <utility>

namespace ad {
namespace {

void require_valid(const Mat &matrix) {
    if (!matrix.valid()) {
        throw std::runtime_error("invalid matrix expression");
    }
}

void require_same_shape(const Mat &lhs, const Mat &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    if (lhs.rows() != rhs.rows() || lhs.cols() != rhs.cols()) {
        throw std::runtime_error("matrix shape mismatch");
    }
}

void require_matmul_shape(const Mat &lhs, const Mat &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    if (lhs.cols() != rhs.rows()) {
        throw std::runtime_error("matrix multiplication shape mismatch");
    }
}

Mat wrap_like(Vec values, const Mat &model) {
    return Mat(std::move(values), model.rows(), model.cols(), model.layout());
}

Vec rhs_values_in_lhs_layout(const Mat &lhs, const Mat &rhs) {
    return rhs.as_layout(lhs.layout()).vec();
}

std::vector<int> indices_for_layout_repack(const Mat &source, MatrixLayout target_layout) {
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(source.size()));
    for (int out = 0; out < source.size(); ++out) {
        const auto coords = detail::matrix_coordinates_from_flat(out, source.rows(), source.cols(), target_layout);
        indices.push_back(detail::matrix_flat_index(coords.first, coords.second, source.rows(), source.cols(), source.layout()));
    }
    return indices;
}

std::vector<int> indices_for_transpose(const Mat &source) {
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(source.size()));
    for (int out = 0; out < source.size(); ++out) {
        const auto coords = detail::matrix_coordinates_from_flat(out, source.cols(), source.rows(), source.layout());
        indices.push_back(detail::matrix_flat_index(coords.second, coords.first, source.rows(), source.cols(), source.layout()));
    }
    return indices;
}

} // namespace

Mat::Mat(Vec values, int rows, int cols, MatrixLayout layout)
    : values_(std::move(values)), rows_(rows), cols_(cols), layout_(layout) {
    if (!values_.valid()) {
        throw std::runtime_error("matrix values must be a valid Vec");
    }
    if (detail::checked_matrix_size(rows_, cols_) != values_.size()) {
        throw std::runtime_error("matrix dimensions must match Vec size");
    }
}

bool Mat::valid() const {
    return values_.valid();
}

int Mat::rows() const {
    return rows_;
}

int Mat::cols() const {
    return cols_;
}

int Mat::size() const {
    return detail::checked_matrix_size(rows_, cols_);
}

MatrixLayout Mat::layout() const {
    return layout_;
}

const Vec &Mat::vec() const {
    return values_;
}

Expr Mat::operator()(int row, int col) const {
    require_valid(*this);
    return values_[detail::matrix_flat_index(row, col, rows_, cols_, layout_)];
}

Vec Mat::row(int index) const {
    require_valid(*this);
    if (index < 0 || index >= rows_) {
        throw std::runtime_error("matrix row index out of bounds");
    }
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(cols_));
    for (int col_index = 0; col_index < cols_; ++col_index) {
        indices.push_back(detail::matrix_flat_index(index, col_index, rows_, cols_, layout_));
    }
    return gather(values_, std::move(indices));
}

Vec Mat::col(int index) const {
    require_valid(*this);
    if (index < 0 || index >= cols_) {
        throw std::runtime_error("matrix column index out of bounds");
    }
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(rows_));
    for (int row_index = 0; row_index < rows_; ++row_index) {
        indices.push_back(detail::matrix_flat_index(row_index, index, rows_, cols_, layout_));
    }
    return gather(values_, std::move(indices));
}

Mat Mat::block(int row, int col, int rows, int cols) const {
    require_valid(*this);
    if (rows < 0 || cols < 0 || row < 0 || col < 0 || row + rows > rows_ || col + cols > cols_) {
        throw std::runtime_error("matrix block out of bounds");
    }
    const int out_size = detail::checked_matrix_size(rows, cols);
    std::vector<int> indices;
    indices.reserve(static_cast<std::size_t>(out_size));
    for (int out = 0; out < out_size; ++out) {
        const auto coords = detail::matrix_coordinates_from_flat(out, rows, cols, layout_);
        indices.push_back(detail::matrix_flat_index(row + coords.first, col + coords.second, rows_, cols_, layout_));
    }
    return Mat(gather(values_, std::move(indices)), rows, cols, layout_);
}

Mat Mat::transpose() const {
    require_valid(*this);
    return Mat(gather(values_, indices_for_transpose(*this)), cols_, rows_, layout_);
}

Mat Mat::reshape(int rows, int cols) const {
    require_valid(*this);
    if (detail::checked_matrix_size(rows, cols) != size()) {
        throw std::runtime_error("matrix reshape must preserve size");
    }
    return Mat(values_, rows, cols, layout_);
}

Mat Mat::reshape(int rows, int cols, MatrixLayout layout) const {
    return reshape(rows, cols).as_layout(layout);
}

Mat Mat::as_layout(MatrixLayout layout) const {
    require_valid(*this);
    if (layout == layout_) {
        return *this;
    }
    return Mat(gather(values_, indices_for_layout_repack(*this, layout)), rows_, cols_, layout);
}

void Mat::eval(const Values &values, EvalWorkspace &workspace, double *out, int out_size) const {
    require_valid(*this);
    values_.eval(values, workspace, out, out_size);
}

Mat mat(Vec values, int rows, int cols, MatrixLayout layout) {
    return Mat(std::move(values), rows, cols, layout);
}

Mat mat_variable(const std::string &label, int rows, int cols, MatrixLayout layout) {
    return Mat(vec_variable(label, detail::checked_matrix_size(rows, cols)), rows, cols, layout);
}

Mat mat_parameter(const std::string &label, int rows, int cols, MatrixLayout layout) {
    return Mat(vec_parameter(label, detail::checked_matrix_size(rows, cols)), rows, cols, layout);
}

Mat mat_constant(int rows, int cols, const std::vector<double> &values, MatrixLayout layout) {
    return Mat(vec_constant(values), rows, cols, layout);
}

Mat operator+(const Mat &lhs, const Mat &rhs) {
    require_same_shape(lhs, rhs);
    return wrap_like(lhs.vec() + rhs_values_in_lhs_layout(lhs, rhs), lhs);
}

Mat operator-(const Mat &lhs, const Mat &rhs) {
    require_same_shape(lhs, rhs);
    return wrap_like(lhs.vec() - rhs_values_in_lhs_layout(lhs, rhs), lhs);
}

Mat operator*(const Mat &lhs, const Mat &rhs) {
    require_matmul_shape(lhs, rhs);
    return Mat(detail::make_symbolic_matmul(lhs.vec(),
                                            lhs.rows(),
                                            lhs.cols(),
                                            lhs.layout(),
                                            rhs.vec(),
                                            rhs.rows(),
                                            rhs.cols(),
                                            rhs.layout(),
                                            lhs.layout()),
               lhs.rows(),
               rhs.cols(),
               lhs.layout());
}

Vec operator*(const Mat &lhs, const Vec &rhs) {
    require_valid(lhs);
    return detail::make_symbolic_matvec(lhs.vec(), lhs.rows(), lhs.cols(), lhs.layout(), rhs);
}

Mat operator+(const Mat &lhs, double rhs) {
    require_valid(lhs);
    return wrap_like(lhs.vec() + rhs, lhs);
}

Mat operator-(const Mat &lhs, double rhs) {
    require_valid(lhs);
    return wrap_like(lhs.vec() - rhs, lhs);
}

Mat operator*(const Mat &lhs, double rhs) {
    require_valid(lhs);
    return wrap_like(lhs.vec() * rhs, lhs);
}

Mat operator/(const Mat &lhs, double rhs) {
    require_valid(lhs);
    return wrap_like(lhs.vec() / rhs, lhs);
}

Mat operator+(double lhs, const Mat &rhs) {
    require_valid(rhs);
    return wrap_like(lhs + rhs.vec(), rhs);
}

Mat operator-(double lhs, const Mat &rhs) {
    require_valid(rhs);
    return wrap_like(lhs - rhs.vec(), rhs);
}

Mat operator*(double lhs, const Mat &rhs) {
    require_valid(rhs);
    return wrap_like(lhs * rhs.vec(), rhs);
}

Mat operator/(double lhs, const Mat &rhs) {
    require_valid(rhs);
    return wrap_like(lhs / rhs.vec(), rhs);
}

Mat operator*(const Expr &lhs, const Mat &rhs) {
    require_valid(rhs);
    return wrap_like(lhs * rhs.vec(), rhs);
}

Mat operator*(const Mat &lhs, const Expr &rhs) {
    require_valid(lhs);
    return wrap_like(lhs.vec() * rhs, lhs);
}

Mat operator-(const Mat &expr) {
    require_valid(expr);
    return wrap_like(-expr.vec(), expr);
}

Mat sin(const Mat &expr) {
    require_valid(expr);
    return wrap_like(sin(expr.vec()), expr);
}

Mat cos(const Mat &expr) {
    require_valid(expr);
    return wrap_like(cos(expr.vec()), expr);
}

Mat tan(const Mat &expr) {
    require_valid(expr);
    return wrap_like(tan(expr.vec()), expr);
}

Mat exp(const Mat &expr) {
    require_valid(expr);
    return wrap_like(exp(expr.vec()), expr);
}

Mat log(const Mat &expr) {
    require_valid(expr);
    return wrap_like(log(expr.vec()), expr);
}

Mat sigmoid(const Mat &expr) {
    require_valid(expr);
    return wrap_like(sigmoid(expr.vec()), expr);
}

Mat abs(const Mat &expr) {
    require_valid(expr);
    return wrap_like(abs(expr.vec()), expr);
}

Mat sqrt(const Mat &expr) {
    require_valid(expr);
    return wrap_like(sqrt(expr.vec()), expr);
}

Mat asin(const Mat &expr) {
    require_valid(expr);
    return wrap_like(asin(expr.vec()), expr);
}

Mat acos(const Mat &expr) {
    require_valid(expr);
    return wrap_like(acos(expr.vec()), expr);
}

Mat atan(const Mat &expr) {
    require_valid(expr);
    return wrap_like(atan(expr.vec()), expr);
}

Mat sinh(const Mat &expr) {
    require_valid(expr);
    return wrap_like(sinh(expr.vec()), expr);
}

Mat cosh(const Mat &expr) {
    require_valid(expr);
    return wrap_like(cosh(expr.vec()), expr);
}

Mat tanh(const Mat &expr) {
    require_valid(expr);
    return wrap_like(tanh(expr.vec()), expr);
}

Mat log10(const Mat &expr) {
    require_valid(expr);
    return wrap_like(log10(expr.vec()), expr);
}

Mat pow(const Mat &base, double exponent) {
    require_valid(base);
    return wrap_like(pow(base.vec(), exponent), base);
}

Mat hadamard(const Mat &lhs, const Mat &rhs) {
    return elem_mul(lhs, rhs);
}

Mat elem_mul(const Mat &lhs, const Mat &rhs) {
    require_same_shape(lhs, rhs);
    return wrap_like(lhs.vec() * rhs_values_in_lhs_layout(lhs, rhs), lhs);
}

Mat elem_div(const Mat &lhs, const Mat &rhs) {
    require_same_shape(lhs, rhs);
    return wrap_like(lhs.vec() / rhs_values_in_lhs_layout(lhs, rhs), lhs);
}

Mat elem_pow(const Mat &base, const Mat &exponent) {
    require_same_shape(base, exponent);
    return wrap_like(pow(base.vec(), rhs_values_in_lhs_layout(base, exponent)), base);
}

Mat elem_min(const Mat &lhs, const Mat &rhs) {
    require_same_shape(lhs, rhs);
    return wrap_like(min(lhs.vec(), rhs_values_in_lhs_layout(lhs, rhs)), lhs);
}

Mat elem_max(const Mat &lhs, const Mat &rhs) {
    require_same_shape(lhs, rhs);
    return wrap_like(max(lhs.vec(), rhs_values_in_lhs_layout(lhs, rhs)), lhs);
}

Mat outer(const Vec &lhs, const Vec &rhs, MatrixLayout layout) {
    return Mat(detail::make_outer_product(lhs, rhs, layout), lhs.size(), rhs.size(), layout);
}

Expr bilinear_form(const Vec &lhs, const Mat &matrix, const Vec &rhs) {
    require_valid(matrix);
    if (lhs.size() != matrix.rows()) {
        throw std::runtime_error("bilinear form lhs size must match matrix row count");
    }
    return dot(lhs, matrix * rhs);
}

Expr quadratic_form(const Vec &values, const Mat &matrix) {
    return bilinear_form(values, matrix, values);
}

Vec solve(const Mat &matrix, const Vec &rhs, LinearSolveOptions options) {
    require_valid(matrix);
    if (!rhs.valid()) {
        throw std::runtime_error("linear solve rhs must be a valid Vec");
    }
    return detail::make_linear_solve(matrix.vec(),
                                     matrix.rows(),
                                     matrix.cols(),
                                     matrix.layout(),
                                     rhs,
                                     options,
                                     false);
}

Vec solve_transpose(const Mat &matrix, const Vec &rhs, LinearSolveOptions options) {
    require_valid(matrix);
    if (!rhs.valid()) {
        throw std::runtime_error("linear transpose solve rhs must be a valid Vec");
    }
    return detail::make_linear_solve(matrix.vec(),
                                     matrix.rows(),
                                     matrix.cols(),
                                     matrix.layout(),
                                     rhs,
                                     options,
                                     true);
}

} // namespace ad
