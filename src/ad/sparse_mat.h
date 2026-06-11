// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_SPARSE_MAT_H
#define MOO_AD_SPARSE_MAT_H

#include "mat.h"

#include <string>
#include <vector>

namespace ad {

class SparseMat {
public:
    SparseMat() = default;
    SparseMat(Vec values, int rows, int cols, std::vector<int> row, std::vector<int> col);

    bool valid() const;
    int rows() const;
    int cols() const;
    int nnz() const;

    const Vec &values() const;
    const std::vector<int> &row_indices() const;
    const std::vector<int> &col_indices() const;

    Expr operator()(int row, int col) const;
    SparseMat transpose() const;
    Mat to_dense(MatrixLayout layout = MatrixLayout::ColumnMajor) const;

    void eval(const Values &values, EvalWorkspace &workspace, double *out, int out_size) const;

private:
    Vec values_;
    int rows_ = 0;
    int cols_ = 0;
    std::vector<int> row_;
    std::vector<int> col_;
};

SparseMat sparse_mat(Vec values, int rows, int cols, std::vector<int> row, std::vector<int> col);
SparseMat sparse_mat_variable(const std::string &label, int rows, int cols, std::vector<int> row, std::vector<int> col);
SparseMat sparse_mat_parameter(const std::string &label, int rows, int cols, std::vector<int> row, std::vector<int> col);
SparseMat sparse_mat_constant(int rows, int cols, std::vector<int> row, std::vector<int> col, std::vector<double> values);

SparseMat operator+(const SparseMat &lhs, const SparseMat &rhs);
SparseMat operator-(const SparseMat &lhs, const SparseMat &rhs);
SparseMat operator*(const Expr &lhs, const SparseMat &rhs);
SparseMat operator*(const SparseMat &lhs, const Expr &rhs);
SparseMat operator*(double lhs, const SparseMat &rhs);
SparseMat operator*(const SparseMat &lhs, double rhs);
SparseMat operator-(const SparseMat &expr);
SparseMat hadamard(const SparseMat &lhs, const SparseMat &rhs);

Vec operator*(const SparseMat &lhs, const Vec &rhs);
Mat operator*(const SparseMat &lhs, const Mat &rhs);
Mat operator*(const Mat &lhs, const SparseMat &rhs);

} // namespace ad

#endif // MOO_AD_SPARSE_MAT_H
