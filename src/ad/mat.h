// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_MAT_H
#define MOO_AD_MAT_H

#include "vec.h"

#include <string>
#include <vector>

namespace ad {

class Mat {
public:
    Mat() = default;
    Mat(Vec values, int rows, int cols, MatrixLayout layout = MatrixLayout::ColumnMajor);

    bool valid() const;
    int rows() const;
    int cols() const;
    int size() const;
    MatrixLayout layout() const;
    const Vec &vec() const;

    Expr operator()(int row, int col) const;
    Vec row(int index) const;
    Vec col(int index) const;
    Mat block(int row, int col, int rows, int cols) const;
    Mat transpose() const;
    Mat reshape(int rows, int cols) const;
    Mat reshape(int rows, int cols, MatrixLayout layout) const;
    Mat as_layout(MatrixLayout layout) const;

    void eval(const Values &values, EvalWorkspace &workspace, double *out, int out_size) const;

private:
    Vec values_;
    int rows_ = 0;
    int cols_ = 0;
    MatrixLayout layout_ = MatrixLayout::ColumnMajor;
};

Mat mat(Vec values, int rows, int cols, MatrixLayout layout = MatrixLayout::ColumnMajor);
Mat mat_variable(const std::string &label, int rows, int cols, MatrixLayout layout = MatrixLayout::ColumnMajor);
Mat mat_parameter(const std::string &label, int rows, int cols, MatrixLayout layout = MatrixLayout::ColumnMajor);
Mat mat_constant(int rows, int cols, const std::vector<double> &values, MatrixLayout layout = MatrixLayout::ColumnMajor);

Mat operator+(const Mat &lhs, const Mat &rhs);
Mat operator-(const Mat &lhs, const Mat &rhs);
Mat operator*(const Mat &lhs, const Mat &rhs);
Vec operator*(const Mat &lhs, const Vec &rhs);

Mat operator+(const Mat &lhs, double rhs);
Mat operator-(const Mat &lhs, double rhs);
Mat operator*(const Mat &lhs, double rhs);
Mat operator/(const Mat &lhs, double rhs);
Mat operator+(double lhs, const Mat &rhs);
Mat operator-(double lhs, const Mat &rhs);
Mat operator*(double lhs, const Mat &rhs);
Mat operator/(double lhs, const Mat &rhs);

Mat operator*(const Expr &lhs, const Mat &rhs);
Mat operator*(const Mat &lhs, const Expr &rhs);
Mat operator-(const Mat &expr);

Mat sin(const Mat &expr);
Mat cos(const Mat &expr);
Mat tan(const Mat &expr);
Mat exp(const Mat &expr);
Mat log(const Mat &expr);
Mat sigmoid(const Mat &expr);
Mat abs(const Mat &expr);
Mat sqrt(const Mat &expr);
Mat asin(const Mat &expr);
Mat acos(const Mat &expr);
Mat atan(const Mat &expr);
Mat sinh(const Mat &expr);
Mat cosh(const Mat &expr);
Mat tanh(const Mat &expr);
Mat log10(const Mat &expr);
Mat pow(const Mat &base, double exponent);
Mat hadamard(const Mat &lhs, const Mat &rhs);
Mat elem_mul(const Mat &lhs, const Mat &rhs);
Mat elem_div(const Mat &lhs, const Mat &rhs);
Mat elem_pow(const Mat &base, const Mat &exponent);
Mat elem_min(const Mat &lhs, const Mat &rhs);
Mat elem_max(const Mat &lhs, const Mat &rhs);
Mat outer(const Vec &lhs, const Vec &rhs, MatrixLayout layout = MatrixLayout::ColumnMajor);
Expr bilinear_form(const Vec &lhs, const Mat &matrix, const Vec &rhs);
Expr quadratic_form(const Vec &values, const Mat &matrix);
Vec solve(const Mat &matrix, const Vec &rhs, LinearSolveOptions options = {});
Vec solve_transpose(const Mat &matrix, const Vec &rhs, LinearSolveOptions options = {});

} // namespace ad

#endif // MOO_AD_MAT_H
