// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2026 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef MOO_AD_CORE_H
#define MOO_AD_CORE_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ad {

struct Graph;

using NodeId = int;

inline constexpr NodeId invalid_node = -1;

bool nearly_zero(double x);
bool exactly(double a, double b);

enum class Op {
    Constant,
    Input,
    Param,
    Add,
    Sub,
    Mul,
    Div,
    Neg,
    Sin,
    Cos,
    Tan,
    Exp,
    Log,
    PowConst
};

const char *op_name(Op op);

struct Node {
    Op op = Op::Constant;
    NodeId a = invalid_node;
    NodeId b = invalid_node;
    double value = 0.0;
    std::string name;
    int index = -1;
};

struct Expr {
    Graph *g = nullptr;
    NodeId id = invalid_node;

    Expr() = default;
    Expr(Graph *graph, NodeId node);

    explicit operator bool() const;
};

struct NamedVector {
    std::string name;
    std::vector<Expr> values;

    std::size_t size() const;
    Expr &operator[](std::size_t i);
    const Expr &operator[](std::size_t i) const;
};

struct DenseMatrix {
    int rows = 0;
    int cols = 0;
    std::vector<double> values;

    DenseMatrix() = default;
    DenseMatrix(int rows_in, int cols_in, std::vector<double> values_in);

    double operator()(int row, int col) const;
};

struct SparseMatrix {
    int rows = 0;
    int cols = 0;
    std::vector<int> row_indices;
    std::vector<int> col_indices;
    std::vector<double> values;

    SparseMatrix() = default;
    SparseMatrix(int rows_in, int cols_in, std::vector<int> rows_in_values, std::vector<int> cols_in_values, std::vector<double> values_in);

    int nnz() const;
};

struct Graph {
    std::vector<Node> nodes;

    Expr constant(double v);
    Expr input(const std::string &name, int index);
    Expr param(const std::string &name, int index);
    Expr unary(Op op, Expr x);
    Expr binary(Op op, Expr x, Expr y);
    Expr pow_const(Expr x, double p);
    void check_same(Expr e) const;
    NamedVector inputs(const std::string &name, int n);
    NamedVector params(const std::string &name, int n);
};

void require_same_graph(Expr a, Expr b);
Expr operator+(Expr a, Expr b);
Expr operator-(Expr a, Expr b);
Expr operator*(Expr a, Expr b);
Expr operator/(Expr a, Expr b);
Expr operator-(Expr a);
Expr operator+(Expr a, double b);
Expr operator-(Expr a, double b);
Expr operator*(Expr a, double b);
Expr operator/(Expr a, double b);
Expr operator+(double a, Expr b);
Expr operator-(double a, Expr b);
Expr operator*(double a, Expr b);
Expr operator/(double a, Expr b);
Expr sin(Expr x);
Expr cos(Expr x);
Expr tan(Expr x);
Expr exp(Expr x);
Expr log(Expr x);
Expr pow_const(Expr x, double p);
Expr sqr(Expr x);

void require_vector_size(const std::vector<Expr> &v, int expected, const std::string &what);
void require_same_graph(const std::vector<Expr> &v);
Expr balanced_sum(Graph &g, const std::vector<Expr> &terms);
std::vector<Expr> dense_matvec(const DenseMatrix &matrix, const std::vector<Expr> &x);
std::vector<Expr> sparse_matvec(const SparseMatrix &matrix, const std::vector<Expr> &x);
std::vector<Expr> kron_eye_matvec(const DenseMatrix &base, int eye_size, const std::vector<Expr> &x);
std::vector<Expr> vector_add(const std::vector<Expr> &a, const std::vector<Expr> &b);
std::vector<Expr> vector_sub(const std::vector<Expr> &a, const std::vector<Expr> &b);
std::vector<Expr> vector_scale(double factor, const std::vector<Expr> &x);

enum class VectorOp {
    Values,
    Add,
    Sub,
    Scale,
    DenseMatVec,
    SparseMatVec,
    KronEyeMatVec,
    Slice,
    Concat,
};

struct FunctionVectorNode {
    VectorOp op = VectorOp::Values;
    int size = 0;
    std::vector<NodeId> values;
    int a = -1;
    int b = -1;
    double scale = 1.0;
    DenseMatrix matrix;
    SparseMatrix sparse_matrix;
    int eye_size = 0;
    int start = 0;
    int stride = 1;
    std::string name;
    bool is_input_group = false;
    bool is_param_group = false;
};

struct GraphFunction {
    Graph graph;
    std::vector<NodeId> inputs;
    std::vector<NodeId> params;
    std::vector<NodeId> outputs;

    std::vector<std::pair<std::string, int>> input_groups;
    std::vector<std::pair<std::string, int>> param_groups;
    std::vector<FunctionVectorNode> vector_nodes;
    int input_vector = -1;
    int output_vector = -1;
    std::vector<int> param_vectors;
    bool vector_structure_valid = false;

    int input_size() const;
    int param_size() const;
    int output_size() const;
    int vector_node_count() const;
    bool has_vector_structure() const;
    int input_group_size(const std::string &name) const;
    int param_group_size(const std::string &name) const;
};

GraphFunction optimize(const GraphFunction &f);

void add_unique_group(std::vector<std::pair<std::string, int>> &groups, const std::string &name, int size);
int input_group_size(const GraphFunction &f, const std::string &name);
int param_group_size(const GraphFunction &f, const std::string &name);
GraphFunction function_from(Graph &&graph, const NamedVector &inputs, const std::vector<Expr> &outputs, const std::vector<NamedVector> &params = {});
GraphFunction function_from(Graph &&graph, const NamedVector &inputs, const NamedVector &outputs, const std::vector<NamedVector> &params = {});

struct GraphFunctionBuilder;

struct VectorExpr {
    GraphFunctionBuilder *builder = nullptr;
    int id = -1;
    int size = 0;

    VectorExpr() = default;
    VectorExpr(GraphFunctionBuilder *builder_in, int id_in, int size_in);

    explicit operator bool() const;
};

struct VectorNode {
    VectorOp op = VectorOp::Values;
    int size = 0;
    std::vector<Expr> values;
    int a = -1;
    int b = -1;
    double scale = 1.0;
    DenseMatrix matrix;
    SparseMatrix sparse_matrix;
    int eye_size = 0;
    int start = 0;
    int stride = 1;
    std::string name;
    bool is_input_group = false;
    bool is_param_group = false;
};

struct GraphFunctionBuilder {
    Graph graph;
    std::vector<VectorNode> vectors;

    Expr constant(double value);
    VectorExpr inputs(const std::string &name, int size);
    VectorExpr params(const std::string &name, int size);
    VectorExpr vector(std::vector<Expr> values);
    VectorExpr add(VectorExpr lhs, VectorExpr rhs);
    VectorExpr sub(VectorExpr lhs, VectorExpr rhs);
    VectorExpr scale(double factor, VectorExpr rhs);
    VectorExpr dense_matvec(DenseMatrix matrix, VectorExpr rhs);
    VectorExpr sparse_matvec(SparseMatrix matrix, VectorExpr rhs);
    VectorExpr kron_eye_matvec(DenseMatrix base, int eye_size, VectorExpr rhs);
    VectorExpr slice(VectorExpr rhs, int start, int length, int stride = 1);
    VectorExpr concat(VectorExpr lhs, VectorExpr rhs);
    Expr at(VectorExpr rhs, int index);
    GraphFunction function(VectorExpr inputs, VectorExpr outputs, const std::vector<VectorExpr> &params = {});

private:
    VectorExpr add_vector(VectorNode node);
    void require_vector(VectorExpr expr) const;
    std::vector<FunctionVectorNode> freeze_vector_nodes() const;
    std::vector<Expr> lower(VectorExpr expr, std::vector<std::vector<Expr>> &cache);
};

} // namespace ad

#endif // MOO_AD_CORE_H
