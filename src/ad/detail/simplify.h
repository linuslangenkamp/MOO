// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_DETAIL_SIMPLIFY_H
#define MOO_AD_DETAIL_SIMPLIFY_H

#include "../expr.h"
#include "../matrix.h"
#include "../vec.h"
#include "node.h"

#include <memory>
#include <string>
#include <vector>

namespace ad::detail {

Expr make_scalar_constant(double value);
Expr make_scalar_variable(const std::string &label);
Expr make_scalar_parameter(const std::string &label);
Expr make_scalar_unary(ScalarUnaryOp op, const Expr &arg);
Expr make_scalar_pow_const(const Expr &arg, double exponent);
Expr make_scalar_binary(ScalarBinaryOp op, const Expr &lhs, const Expr &rhs);
Expr make_vector_element(const Vec &vec, int index);
Expr make_sum(const Vec &values);
Expr make_dot(const Vec &lhs, const Vec &rhs);

Vec make_vec_from_elements(std::vector<Expr> elements);
Vec make_vec_variable(const std::string &label, int size);
Vec make_vec_parameter(const std::string &label, int size);
Vec make_vec_constant(std::vector<double> values);
Vec make_vec_unary(VecUnaryOp op, const Vec &arg);
Vec make_vec_binary(VecBinaryOp op, const Vec &lhs, const Vec &rhs);
Vec make_vec_scale(const Expr &scale, const Vec &vec);
Vec make_dense_matvec(const DenseMatrix &matrix, const Vec &rhs);
Vec make_sparse_matvec(const SparseMatrix &matrix, const Vec &rhs);
Vec make_symbolic_matvec(const Vec &matrix, int rows, int cols, MatrixLayout layout, const Vec &rhs);
Vec make_symbolic_matmul(const Vec &lhs, int lhs_rows, int lhs_cols, MatrixLayout lhs_layout,
                         const Vec &rhs, int rhs_rows, int rhs_cols, MatrixLayout rhs_layout,
                         MatrixLayout result_layout);
Vec make_symbolic_sparse_matvec(const Vec &values, int rows, int cols, std::vector<int> row, std::vector<int> col, const Vec &rhs);
Vec make_symbolic_sparse_matmul(const Vec &sparse_values, int sparse_rows, int sparse_cols,
                                std::vector<int> sparse_row, std::vector<int> sparse_col,
                                const Vec &dense, int dense_rows, int dense_cols, MatrixLayout dense_layout,
                                MatrixLayout result_layout, bool sparse_lhs);
Vec make_outer_product(const Vec &lhs, const Vec &rhs, MatrixLayout result_layout);
Vec make_linear_solve(const Vec &matrix, int rows, int cols, MatrixLayout layout, const Vec &rhs,
                      LinearSolveOptions options, bool transpose);
Vec make_slice(const Vec &source, int start, int length);
Vec make_scatter_slice(const Vec &values, int start, int output_size);
Vec make_concat(const Vec &lhs, const Vec &rhs);
Vec make_gather(const Vec &source, std::vector<int> indices);
Vec make_scatter_add(const Vec &values, std::vector<int> indices, int output_size);

Expr simplify_expr(const Expr &expr);
Vec simplify_vec(const Vec &vec);

} // namespace ad::detail

#endif // MOO_AD_DETAIL_SIMPLIFY_H
