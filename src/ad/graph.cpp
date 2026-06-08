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
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "core.h"

namespace ad {

bool nearly_zero(double x) {
    return std::abs(x) <= 0.0;
}

bool exactly(double a, double b) {
    return a == b;
}

const char *op_name(Op op) {
    switch (op) {
        case Op::Constant:
            return "Constant";
        case Op::Input:
            return "Input";
        case Op::Param:
            return "Param";
        case Op::Add:
            return "Add";
        case Op::Sub:
            return "Sub";
        case Op::Mul:
            return "Mul";
        case Op::Div:
            return "Div";
        case Op::Neg:
            return "Neg";
        case Op::Sin:
            return "Sin";
        case Op::Cos:
            return "Cos";
        case Op::Tan:
            return "Tan";
        case Op::Exp:
            return "Exp";
        case Op::Log:
            return "Log";
        case Op::PowConst:
            return "PowConst";
    }
    return "Unknown";
}

Expr::Expr(Graph *graph, NodeId node)
    : g(graph),
      id(node) {}

Expr::operator bool() const {
    return g != nullptr && id != invalid_node;
}

std::size_t NamedVector::size() const {
    return values.size();
}

Expr &NamedVector::operator[](std::size_t i) {
    return values[i];
}

const Expr &NamedVector::operator[](std::size_t i) const {
    return values[i];
}

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
    return values[static_cast<std::size_t>(row * cols + col)];
}

SparseMatrix::SparseMatrix(int rows_in, int cols_in, std::vector<int> rows_in_values, std::vector<int> cols_in_values, std::vector<double> values_in)
    : rows(rows_in),
      cols(cols_in),
      row_indices(std::move(rows_in_values)),
      col_indices(std::move(cols_in_values)),
      values(std::move(values_in)) {
    if (rows < 0 || cols < 0) {
        throw std::runtime_error("sparse matrix dimensions must be nonnegative");
    }
    if (row_indices.size() != col_indices.size() || row_indices.size() != values.size()) {
        throw std::runtime_error("sparse matrix row, column, and value arrays must have equal length");
    }
    for (int k = 0; k < static_cast<int>(values.size()); ++k) {
        int row = row_indices[static_cast<std::size_t>(k)];
        int col = col_indices[static_cast<std::size_t>(k)];
        if (row < 0 || row >= rows || col < 0 || col >= cols) {
            throw std::runtime_error("sparse matrix entry out of bounds");
        }
    }
}

int SparseMatrix::nnz() const {
    return static_cast<int>(values.size());
}

Expr Graph::constant(double v) {
    Node n;
    n.op = Op::Constant;
    n.value = v;
    nodes.push_back(std::move(n));
    return Expr(this, static_cast<NodeId>(nodes.size() - 1));
}

Expr Graph::input(const std::string &name, int index) {
    Node n;
    n.op = Op::Input;
    n.name = name;
    n.index = index;
    nodes.push_back(std::move(n));
    return Expr(this, static_cast<NodeId>(nodes.size() - 1));
}

Expr Graph::param(const std::string &name, int index) {
    Node n;
    n.op = Op::Param;
    n.name = name;
    n.index = index;
    nodes.push_back(std::move(n));
    return Expr(this, static_cast<NodeId>(nodes.size() - 1));
}

Expr Graph::unary(Op op, Expr x) {
    check_same(x);
    Node n;
    n.op = op;
    n.a = x.id;
    nodes.push_back(std::move(n));
    return Expr(this, static_cast<NodeId>(nodes.size() - 1));
}

Expr Graph::binary(Op op, Expr x, Expr y) {
    check_same(x);
    check_same(y);
    Node n;
    n.op = op;
    n.a = x.id;
    n.b = y.id;
    nodes.push_back(std::move(n));
    return Expr(this, static_cast<NodeId>(nodes.size() - 1));
}

Expr Graph::pow_const(Expr x, double p) {
    check_same(x);
    Node n;
    n.op = Op::PowConst;
    n.a = x.id;
    n.value = p;
    nodes.push_back(std::move(n));
    return Expr(this, static_cast<NodeId>(nodes.size() - 1));
}

void Graph::check_same(Expr e) const {
    if (e.g != this || e.id < 0 || e.id >= static_cast<NodeId>(nodes.size())) {
        throw std::runtime_error("expression belongs to a different graph or is invalid");
    }
}

NamedVector Graph::inputs(const std::string &name, int n) {
    NamedVector v;
    v.name = name;
    for (int i = 0; i < n; ++i) {
        v.values.push_back(input(name, i));
    }
    return v;
}

NamedVector Graph::params(const std::string &name, int n) {
    NamedVector v;
    v.name = name;
    for (int i = 0; i < n; ++i) {
        v.values.push_back(param(name, i));
    }
    return v;
}

void require_same_graph(Expr a, Expr b) {
    if (a.g != b.g) {
        throw std::runtime_error("expressions belong to different graphs");
    }
}

Expr operator+(Expr a, Expr b) {
    require_same_graph(a, b);
    return a.g->binary(Op::Add, a, b);
}

Expr operator-(Expr a, Expr b) {
    require_same_graph(a, b);
    return a.g->binary(Op::Sub, a, b);
}

Expr operator*(Expr a, Expr b) {
    require_same_graph(a, b);
    return a.g->binary(Op::Mul, a, b);
}

Expr operator/(Expr a, Expr b) {
    require_same_graph(a, b);
    return a.g->binary(Op::Div, a, b);
}

Expr operator-(Expr a) {
    return a.g->unary(Op::Neg, a);
}

Expr operator+(Expr a, double b) {
    return a + a.g->constant(b);
}

Expr operator-(Expr a, double b) {
    return a - a.g->constant(b);
}

Expr operator*(Expr a, double b) {
    return a * a.g->constant(b);
}

Expr operator/(Expr a, double b) {
    return a / a.g->constant(b);
}

Expr operator+(double a, Expr b) {
    return b.g->constant(a) + b;
}

Expr operator-(double a, Expr b) {
    return b.g->constant(a) - b;
}

Expr operator*(double a, Expr b) {
    return b.g->constant(a) * b;
}

Expr operator/(double a, Expr b) {
    return b.g->constant(a) / b;
}

Expr sin(Expr x) {
    return x.g->unary(Op::Sin, x);
}

Expr cos(Expr x) {
    return x.g->unary(Op::Cos, x);
}

Expr tan(Expr x) {
    return x.g->unary(Op::Tan, x);
}

Expr exp(Expr x) {
    return x.g->unary(Op::Exp, x);
}

Expr log(Expr x) {
    return x.g->unary(Op::Log, x);
}

Expr pow_const(Expr x, double p) {
    return x.g->pow_const(x, p);
}

Expr sqr(Expr x) {
    return x * x;
}

void require_vector_size(const std::vector<Expr> &v, int expected, const std::string &what) {
    if (static_cast<int>(v.size()) != expected) {
        throw std::runtime_error(what + " has wrong size");
    }
}

void require_same_graph(const std::vector<Expr> &v) {
    if (v.empty()) {
        return;
    }
    Graph *g = v.front().g;
    for (const auto &expr : v) {
        if (expr.g != g) {
            throw std::runtime_error("vector expressions belong to different graphs");
        }
    }
}

Expr balanced_sum(Graph &g, const std::vector<Expr> &terms) {
    if (terms.empty()) {
        return g.constant(0.0);
    }
    std::vector<Expr> layer = terms;
    while (layer.size() > 1) {
        std::vector<Expr> next;
        next.reserve((layer.size() + 1) / 2);
        for (std::size_t i = 0; i < layer.size(); i += 2) {
            if (i + 1 < layer.size()) {
                next.push_back(layer[i] + layer[i + 1]);
            } else {
                next.push_back(layer[i]);
            }
        }
        layer = std::move(next);
    }
    return layer.front();
}

std::vector<Expr> dense_matvec(const DenseMatrix &matrix, const std::vector<Expr> &x) {
    require_vector_size(x, matrix.cols, "dense matvec input");
    require_same_graph(x);
    if (x.empty() && matrix.cols != 0) {
        throw std::runtime_error("dense matvec input is empty");
    }
    Graph *g = x.empty() ? nullptr : x.front().g;
    if (g == nullptr) {
        throw std::runtime_error("dense matvec requires expressions from a graph");
    }

    std::vector<Expr> out;
    out.reserve(static_cast<std::size_t>(matrix.rows));
    for (int row = 0; row < matrix.rows; ++row) {
        std::vector<Expr> terms;
        terms.reserve(static_cast<std::size_t>(matrix.cols));
        for (int col = 0; col < matrix.cols; ++col) {
            double coeff = matrix(row, col);
            if (exactly(coeff, 0.0)) {
                continue;
            }
            terms.push_back(g->constant(coeff) * x[static_cast<std::size_t>(col)]);
        }
        out.push_back(balanced_sum(*g, terms));
    }
    return out;
}

std::vector<Expr> sparse_matvec(const SparseMatrix &matrix, const std::vector<Expr> &x) {
    require_vector_size(x, matrix.cols, "sparse matvec input");
    require_same_graph(x);
    if (x.empty() && matrix.cols != 0) {
        throw std::runtime_error("sparse matvec input is empty");
    }
    Graph *g = x.empty() ? nullptr : x.front().g;
    if (g == nullptr) {
        throw std::runtime_error("sparse matvec requires expressions from a graph");
    }

    std::vector<std::vector<Expr>> terms(static_cast<std::size_t>(matrix.rows));
    for (int k = 0; k < matrix.nnz(); ++k) {
        double coeff = matrix.values[static_cast<std::size_t>(k)];
        if (exactly(coeff, 0.0)) {
            continue;
        }
        int row = matrix.row_indices[static_cast<std::size_t>(k)];
        int col = matrix.col_indices[static_cast<std::size_t>(k)];
        terms[static_cast<std::size_t>(row)].push_back(g->constant(coeff) * x[static_cast<std::size_t>(col)]);
    }

    std::vector<Expr> out;
    out.reserve(static_cast<std::size_t>(matrix.rows));
    for (auto &row_terms : terms) {
        out.push_back(balanced_sum(*g, row_terms));
    }
    return out;
}

std::vector<Expr> kron_eye_matvec(const DenseMatrix &base, int eye_size, const std::vector<Expr> &x) {
    if (eye_size <= 0) {
        throw std::runtime_error("kron-eye matvec identity size must be positive");
    }
    require_vector_size(x, base.cols * eye_size, "kron-eye matvec input");
    require_same_graph(x);
    if (x.empty() && base.cols * eye_size != 0) {
        throw std::runtime_error("kron-eye matvec input is empty");
    }
    Graph *g = x.empty() ? nullptr : x.front().g;
    if (g == nullptr) {
        throw std::runtime_error("kron-eye matvec requires expressions from a graph");
    }

    std::vector<Expr> out;
    out.reserve(static_cast<std::size_t>(base.rows * eye_size));
    for (int row = 0; row < base.rows; ++row) {
        for (int inner = 0; inner < eye_size; ++inner) {
            std::vector<Expr> terms;
            terms.reserve(static_cast<std::size_t>(base.cols));
            for (int col = 0; col < base.cols; ++col) {
                double coeff = base(row, col);
                if (exactly(coeff, 0.0)) {
                    continue;
                }
                terms.push_back(g->constant(coeff) * x[static_cast<std::size_t>(col * eye_size + inner)]);
            }
            out.push_back(balanced_sum(*g, terms));
        }
    }
    return out;
}

std::vector<Expr> vector_add(const std::vector<Expr> &a, const std::vector<Expr> &b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("vector sizes must match");
    }
    require_same_graph(a);
    require_same_graph(b);
    std::vector<Expr> out;
    out.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        require_same_graph(a[i], b[i]);
        out.push_back(a[i] + b[i]);
    }
    return out;
}

std::vector<Expr> vector_sub(const std::vector<Expr> &a, const std::vector<Expr> &b) {
    if (a.size() != b.size()) {
        throw std::runtime_error("vector sizes must match");
    }
    require_same_graph(a);
    require_same_graph(b);
    std::vector<Expr> out;
    out.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        require_same_graph(a[i], b[i]);
        out.push_back(a[i] - b[i]);
    }
    return out;
}

std::vector<Expr> vector_scale(double factor, const std::vector<Expr> &x) {
    require_same_graph(x);
    std::vector<Expr> out;
    out.reserve(x.size());
    for (const auto &expr : x) {
        out.push_back(factor * expr);
    }
    return out;
}

int GraphFunction::input_size() const {
    return static_cast<int>(inputs.size());
}

int GraphFunction::param_size() const {
    return static_cast<int>(params.size());
}

int GraphFunction::output_size() const {
    return static_cast<int>(outputs.size());
}

int GraphFunction::vector_node_count() const {
    return static_cast<int>(vector_nodes.size());
}

bool GraphFunction::has_vector_structure() const {
    return vector_structure_valid && output_vector >= 0;
}

int GraphFunction::input_group_size(const std::string &name) const {
    for (auto [n, s] : input_groups) {
        if (n == name) {
            return s;
        }
    }
    return -1;
}

int GraphFunction::param_group_size(const std::string &name) const {
    for (auto [n, s] : param_groups) {
        if (n == name) {
            return s;
        }
    }
    return -1;
}

void add_unique_group(std::vector<std::pair<std::string, int>> &groups, const std::string &name, int size) {
    for (auto &g : groups) {
        if (g.first == name) {
            g.second = std::max(g.second, size);
            return;
        }
    }
    groups.push_back({name, size});
}

int input_group_size(const GraphFunction &f, const std::string &name) {
    return std::max(f.input_group_size(name), 0);
}

int param_group_size(const GraphFunction &f, const std::string &name) {
    return std::max(f.param_group_size(name), 0);
}

GraphFunction function_from(Graph &&graph, const NamedVector &inputs, const std::vector<Expr> &outputs, const std::vector<NamedVector> &params) {
    GraphFunction f;
    f.graph = std::move(graph);
    for (auto e : inputs.values) {
        f.inputs.push_back(e.id);
    }
    for (auto e : outputs) {
        f.outputs.push_back(e.id);
    }
    add_unique_group(f.input_groups, inputs.name, static_cast<int>(inputs.size()));
    for (const auto &p : params) {
        add_unique_group(f.param_groups, p.name, static_cast<int>(p.size()));
        for (auto e : p.values) {
            f.params.push_back(e.id);
        }
    }
    return optimize(f);
}

GraphFunction function_from(Graph &&graph, const NamedVector &inputs, const NamedVector &outputs, const std::vector<NamedVector> &params) {
    return function_from(std::move(graph), inputs, outputs.values, params);
}

VectorExpr::VectorExpr(GraphFunctionBuilder *builder_in, int id_in, int size_in)
    : builder(builder_in),
      id(id_in),
      size(size_in) {}

VectorExpr::operator bool() const {
    return builder != nullptr && id >= 0;
}

Expr GraphFunctionBuilder::constant(double value) {
    return graph.constant(value);
}

VectorExpr GraphFunctionBuilder::inputs(const std::string &name, int size) {
    auto named = graph.inputs(name, size);
    VectorNode node;
    node.op = VectorOp::Values;
    node.size = size;
    node.values = std::move(named.values);
    node.name = name;
    node.is_input_group = true;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::params(const std::string &name, int size) {
    auto named = graph.params(name, size);
    VectorNode node;
    node.op = VectorOp::Values;
    node.size = size;
    node.values = std::move(named.values);
    node.name = name;
    node.is_param_group = true;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::vector(std::vector<Expr> values) {
    require_same_graph(values);
    VectorNode node;
    node.op = VectorOp::Values;
    node.size = static_cast<int>(values.size());
    node.values = std::move(values);
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::add(VectorExpr lhs, VectorExpr rhs) {
    require_vector(lhs);
    require_vector(rhs);
    if (lhs.size != rhs.size) {
        throw std::runtime_error("vector sizes must match");
    }
    VectorNode node;
    node.op = VectorOp::Add;
    node.size = lhs.size;
    node.a = lhs.id;
    node.b = rhs.id;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::sub(VectorExpr lhs, VectorExpr rhs) {
    require_vector(lhs);
    require_vector(rhs);
    if (lhs.size != rhs.size) {
        throw std::runtime_error("vector sizes must match");
    }
    VectorNode node;
    node.op = VectorOp::Sub;
    node.size = lhs.size;
    node.a = lhs.id;
    node.b = rhs.id;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::scale(double factor, VectorExpr rhs) {
    require_vector(rhs);
    VectorNode node;
    node.op = VectorOp::Scale;
    node.size = rhs.size;
    node.a = rhs.id;
    node.scale = factor;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::dense_matvec(DenseMatrix matrix, VectorExpr rhs) {
    require_vector(rhs);
    if (matrix.cols != rhs.size) {
        throw std::runtime_error("dense matrix/vector dimensions do not match");
    }
    VectorNode node;
    node.op = VectorOp::DenseMatVec;
    node.size = matrix.rows;
    node.a = rhs.id;
    node.matrix = std::move(matrix);
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::sparse_matvec(SparseMatrix matrix, VectorExpr rhs) {
    require_vector(rhs);
    if (matrix.cols != rhs.size) {
        throw std::runtime_error("sparse matrix/vector dimensions do not match");
    }
    VectorNode node;
    node.op = VectorOp::SparseMatVec;
    node.size = matrix.rows;
    node.a = rhs.id;
    node.sparse_matrix = std::move(matrix);
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::kron_eye_matvec(DenseMatrix base, int eye_size, VectorExpr rhs) {
    require_vector(rhs);
    if (eye_size <= 0) {
        throw std::runtime_error("kron-eye matvec identity size must be positive");
    }
    if (base.cols * eye_size != rhs.size) {
        throw std::runtime_error("kron-eye matrix/vector dimensions do not match");
    }
    VectorNode node;
    node.op = VectorOp::KronEyeMatVec;
    node.size = base.rows * eye_size;
    node.a = rhs.id;
    node.matrix = std::move(base);
    node.eye_size = eye_size;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::slice(VectorExpr rhs, int start, int length, int stride) {
    require_vector(rhs);
    if (length < 0 || stride <= 0 || start < 0 || (length > 0 && start + (length - 1) * stride >= rhs.size)) {
        throw std::runtime_error("vector slice out of bounds");
    }
    VectorNode node;
    node.op = VectorOp::Slice;
    node.size = length;
    node.a = rhs.id;
    node.start = start;
    node.stride = stride;
    return add_vector(std::move(node));
}

VectorExpr GraphFunctionBuilder::concat(VectorExpr lhs, VectorExpr rhs) {
    require_vector(lhs);
    require_vector(rhs);
    VectorNode node;
    node.op = VectorOp::Concat;
    node.size = lhs.size + rhs.size;
    node.a = lhs.id;
    node.b = rhs.id;
    return add_vector(std::move(node));
}

Expr GraphFunctionBuilder::at(VectorExpr rhs, int index) {
    require_vector(rhs);
    if (index < 0 || index >= rhs.size) {
        throw std::runtime_error("vector index out of bounds");
    }
    std::vector<std::vector<Expr>> cache(vectors.size());
    return lower(rhs, cache)[static_cast<std::size_t>(index)];
}

GraphFunction GraphFunctionBuilder::function(VectorExpr inputs, VectorExpr outputs, const std::vector<VectorExpr> &params) {
    require_vector(inputs);
    require_vector(outputs);
    std::vector<std::vector<Expr>> cache(vectors.size());
    auto input_values = lower(inputs, cache);
    auto output_values = lower(outputs, cache);
    std::vector<std::vector<Expr>> param_values_all;
    param_values_all.reserve(params.size());
    for (const auto &param : params) {
        require_vector(param);
        param_values_all.push_back(lower(param, cache));
    }

    GraphFunction f;
    for (const auto &e : input_values) {
        f.inputs.push_back(e.id);
    }
    for (const auto &e : output_values) {
        f.outputs.push_back(e.id);
    }
    const auto &input_node = vectors[static_cast<std::size_t>(inputs.id)];
    add_unique_group(f.input_groups, input_node.name, inputs.size);
    for (std::size_t i = 0; i < params.size(); ++i) {
        const auto &param = params[i];
        const auto &param_values = param_values_all[i];
        const auto &param_node = vectors[static_cast<std::size_t>(param.id)];
        add_unique_group(f.param_groups, param_node.name, param.size);
        for (const auto &e : param_values) {
            f.params.push_back(e.id);
        }
    }
    f.vector_nodes = freeze_vector_nodes();
    f.input_vector = inputs.id;
    f.output_vector = outputs.id;
    f.param_vectors.reserve(params.size());
    for (const auto &param : params) {
        f.param_vectors.push_back(param.id);
    }
    f.vector_structure_valid = true;
    f.graph = std::move(graph);
    return optimize(f);
}

VectorExpr GraphFunctionBuilder::add_vector(VectorNode node) {
    int id = static_cast<int>(vectors.size());
    int size = node.size;
    vectors.push_back(std::move(node));
    return VectorExpr(this, id, size);
}

void GraphFunctionBuilder::require_vector(VectorExpr expr) const {
    if (expr.builder != this || expr.id < 0 || expr.id >= static_cast<int>(vectors.size())) {
        throw std::runtime_error("vector expression belongs to a different graph function builder or is invalid");
    }
}

std::vector<FunctionVectorNode> GraphFunctionBuilder::freeze_vector_nodes() const {
    std::vector<FunctionVectorNode> out;
    out.reserve(vectors.size());
    for (const auto &node : vectors) {
        FunctionVectorNode frozen;
        frozen.op = node.op;
        frozen.size = node.size;
        frozen.values.reserve(node.values.size());
        for (const auto &expr : node.values) {
            frozen.values.push_back(expr.id);
        }
        frozen.a = node.a;
        frozen.b = node.b;
        frozen.scale = node.scale;
        frozen.matrix = node.matrix;
        frozen.sparse_matrix = node.sparse_matrix;
        frozen.eye_size = node.eye_size;
        frozen.start = node.start;
        frozen.stride = node.stride;
        frozen.name = node.name;
        frozen.is_input_group = node.is_input_group;
        frozen.is_param_group = node.is_param_group;
        out.push_back(std::move(frozen));
    }
    return out;
}

std::vector<Expr> GraphFunctionBuilder::lower(VectorExpr expr, std::vector<std::vector<Expr>> &cache) {
    require_vector(expr);
    auto &cached = cache[static_cast<std::size_t>(expr.id)];
    if (!cached.empty() || expr.size == 0) {
        return cached;
    }
    const auto &node = vectors[static_cast<std::size_t>(expr.id)];
    switch (node.op) {
        case VectorOp::Values:
            cached = node.values;
            break;
        case VectorOp::Add:
            cached = vector_add(lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache),
                                lower(VectorExpr(this, node.b, vectors[static_cast<std::size_t>(node.b)].size), cache));
            break;
        case VectorOp::Sub:
            cached = vector_sub(lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache),
                                lower(VectorExpr(this, node.b, vectors[static_cast<std::size_t>(node.b)].size), cache));
            break;
        case VectorOp::Scale:
            cached = vector_scale(node.scale, lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache));
            break;
        case VectorOp::DenseMatVec:
            cached = ad::dense_matvec(node.matrix, lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache));
            break;
        case VectorOp::SparseMatVec:
            cached = ad::sparse_matvec(node.sparse_matrix, lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache));
            break;
        case VectorOp::KronEyeMatVec:
            cached = ad::kron_eye_matvec(node.matrix, node.eye_size, lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache));
            break;
        case VectorOp::Slice: {
            auto base = lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache);
            cached.reserve(static_cast<std::size_t>(node.size));
            for (int i = 0; i < node.size; ++i) {
                cached.push_back(base[static_cast<std::size_t>(node.start + i * node.stride)]);
            }
            break;
        }
        case VectorOp::Concat: {
            auto lhs = lower(VectorExpr(this, node.a, vectors[static_cast<std::size_t>(node.a)].size), cache);
            auto rhs = lower(VectorExpr(this, node.b, vectors[static_cast<std::size_t>(node.b)].size), cache);
            cached.reserve(lhs.size() + rhs.size());
            cached.insert(cached.end(), lhs.begin(), lhs.end());
            cached.insert(cached.end(), rhs.begin(), rhs.end());
            break;
        }
    }
    return cached;
}

} // namespace ad
