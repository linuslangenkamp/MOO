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

#include "diff.h"

namespace ad {

// -----------------------------------------------------------------------------
// Copy graph into OptimizingBuilder
// -----------------------------------------------------------------------------

std::vector<Expr> clone_nodes(const Graph &src, OptimizingBuilder &dst, const std::optional<std::string> &input_as_param) {
    std::vector<Expr> map(src.nodes.size());
    for (NodeId i = 0; i < static_cast<NodeId>(src.nodes.size()); ++i) {
        const auto &n = src.nodes[i];
        switch (n.op) {
            case Op::Constant:
                map[i] = dst.constant(n.value);
                break;
            case Op::Input:
                if (input_as_param && n.name == *input_as_param) {
                    map[i] = dst.param(n.name, n.index);
                } else {
                    map[i] = dst.input(n.name, n.index);
                }
                break;
            case Op::Param:
                map[i] = dst.param(n.name, n.index);
                break;
            case Op::Add:
                map[i] = dst.add(map[n.a], map[n.b]);
                break;
            case Op::Sub:
                map[i] = dst.sub(map[n.a], map[n.b]);
                break;
            case Op::Mul:
                map[i] = dst.mul(map[n.a], map[n.b]);
                break;
            case Op::Div:
                map[i] = dst.div(map[n.a], map[n.b]);
                break;
            case Op::Neg:
                map[i] = dst.neg(map[n.a]);
                break;
            case Op::Sin:
                map[i] = dst.unary(Op::Sin, map[n.a]);
                break;
            case Op::Cos:
                map[i] = dst.unary(Op::Cos, map[n.a]);
                break;
            case Op::Tan:
                map[i] = dst.unary(Op::Tan, map[n.a]);
                break;
            case Op::Exp:
                map[i] = dst.unary(Op::Exp, map[n.a]);
                break;
            case Op::Log:
                map[i] = dst.unary(Op::Log, map[n.a]);
                break;
            case Op::PowConst:
                map[i] = dst.pow_const(map[n.a], n.value);
                break;
        }
    }
    return map;
}

// -----------------------------------------------------------------------------
// Forward-mode graph transform: JVP.
// -----------------------------------------------------------------------------

namespace {

bool infer_direct_values_group(FunctionVectorNode &node, const Graph &graph) {
    if (node.op != VectorOp::Values || static_cast<int>(node.values.size()) != node.size || node.size <= 0) {
        return false;
    }

    const auto &first = graph.nodes[static_cast<std::size_t>(node.values.front())];
    if (first.op != Op::Input && first.op != Op::Param) {
        return false;
    }

    for (int i = 0; i < node.size; ++i) {
        const auto &scalar = graph.nodes[static_cast<std::size_t>(node.values[static_cast<std::size_t>(i)])];
        if (scalar.op != first.op || scalar.name != first.name || scalar.index != i) {
            return false;
        }
    }

    node.name = first.name;
    node.is_input_group = first.op == Op::Input;
    node.is_param_group = first.op == Op::Param;
    return true;
}

std::vector<FunctionVectorNode> tangent_vector_nodes(const GraphFunction &f, const Graph &graph, const std::vector<Expr> &tangent) {
    std::vector<FunctionVectorNode> out;
    out.reserve(f.vector_nodes.size());
    for (const auto &node : f.vector_nodes) {
        FunctionVectorNode transformed;
        transformed.op = node.op;
        transformed.size = node.size;
        transformed.a = node.a;
        transformed.b = node.b;
        transformed.scale = node.scale;
        transformed.matrix = node.matrix;
        transformed.sparse_matrix = node.sparse_matrix;
        transformed.eye_size = node.eye_size;
        transformed.start = node.start;
        transformed.stride = node.stride;

        if (node.op == VectorOp::Values) {
            transformed.values.reserve(node.values.size());
            for (NodeId value : node.values) {
                if (value < 0 || value >= static_cast<NodeId>(tangent.size())) {
                    throw std::runtime_error("vector metadata references scalar node " + std::to_string(value) + " outside graph of size " + std::to_string(tangent.size()) +
                                             " during forward differentiation");
                }
                transformed.values.push_back(tangent[static_cast<std::size_t>(value)].id);
            }
            infer_direct_values_group(transformed, graph);
        }

        out.push_back(std::move(transformed));
    }
    return out;
}

DenseMatrix transpose(const DenseMatrix &matrix) {
    std::vector<double> values(static_cast<std::size_t>(matrix.rows * matrix.cols), 0.0);
    for (int row = 0; row < matrix.rows; ++row) {
        for (int col = 0; col < matrix.cols; ++col) {
            values[static_cast<std::size_t>(col * matrix.rows + row)] = matrix(row, col);
        }
    }
    return DenseMatrix(matrix.cols, matrix.rows, std::move(values));
}

SparseMatrix transpose(const SparseMatrix &matrix) {
    return SparseMatrix(matrix.cols, matrix.rows, matrix.col_indices, matrix.row_indices, matrix.values);
}

bool direct_values_group(const GraphFunction &f, const FunctionVectorNode &node, Op op, const std::string &name) {
    if (node.op != VectorOp::Values || static_cast<int>(node.values.size()) != node.size || node.size <= 0) {
        return false;
    }
    for (int i = 0; i < node.size; ++i) {
        const auto &scalar = f.graph.nodes[static_cast<std::size_t>(node.values[static_cast<std::size_t>(i)])];
        if (scalar.op != op || scalar.name != name || scalar.index != i) {
            return false;
        }
    }
    return true;
}

std::vector<Expr> scalar_vjp_for_values(const GraphFunction &f,
                                        OptimizingBuilder &b,
                                        const std::vector<Expr> &primal,
                                        const std::vector<NodeId> &values,
                                        const std::vector<Expr> &seeds,
                                        const std::string &wrt_input_name) {
    if (values.size() != seeds.size()) {
        throw std::runtime_error("scalar_vjp_for_values: value and seed sizes must match");
    }

    auto zero = b.constant(0.0);
    std::vector<std::vector<Expr>> adj_terms(f.graph.nodes.size());
    auto add_adj = [&](NodeId id, Expr term) {
        if (id >= 0) {
            adj_terms[static_cast<std::size_t>(id)].push_back(term);
        }
    };
    auto sum_terms = [&](const std::vector<Expr> &terms) -> Expr {
        if (terms.empty()) {
            return zero;
        }
        Expr s = terms[0];
        for (std::size_t i = 1; i < terms.size(); ++i) {
            s = b.add(s, terms[i]);
        }
        return s;
    };

    for (std::size_t i = 0; i < values.size(); ++i) {
        add_adj(values[i], seeds[i]);
    }

    for (NodeId i = static_cast<NodeId>(f.graph.nodes.size()) - 1; i >= 0; --i) {
        if (adj_terms[static_cast<std::size_t>(i)].empty()) {
            if (i == 0) {
                break;
            }
            continue;
        }
        const auto &n = f.graph.nodes[static_cast<std::size_t>(i)];
        Expr adj = sum_terms(adj_terms[static_cast<std::size_t>(i)]);
        switch (n.op) {
            case Op::Constant:
            case Op::Input:
            case Op::Param:
                break;
            case Op::Add:
                add_adj(n.a, adj);
                add_adj(n.b, adj);
                break;
            case Op::Sub:
                add_adj(n.a, adj);
                add_adj(n.b, b.neg(adj));
                break;
            case Op::Mul:
                add_adj(n.a, b.mul(adj, primal[static_cast<std::size_t>(n.b)]));
                add_adj(n.b, b.mul(adj, primal[static_cast<std::size_t>(n.a)]));
                break;
            case Op::Div:
                add_adj(n.a, b.div(adj, primal[static_cast<std::size_t>(n.b)]));
                add_adj(n.b, b.neg(b.div(b.mul(adj, primal[static_cast<std::size_t>(n.a)]),
                                         b.mul(primal[static_cast<std::size_t>(n.b)], primal[static_cast<std::size_t>(n.b)]))));
                break;
            case Op::Neg:
                add_adj(n.a, b.neg(adj));
                break;
            case Op::Sin:
                add_adj(n.a, b.mul(adj, b.unary(Op::Cos, primal[static_cast<std::size_t>(n.a)])));
                break;
            case Op::Cos:
                add_adj(n.a, b.neg(b.mul(adj, b.unary(Op::Sin, primal[static_cast<std::size_t>(n.a)]))));
                break;
            case Op::Tan:
                add_adj(n.a, b.mul(adj, b.add(b.constant(1.0), b.mul(primal[static_cast<std::size_t>(i)], primal[static_cast<std::size_t>(i)]))));
                break;
            case Op::Exp:
                add_adj(n.a, b.mul(adj, primal[static_cast<std::size_t>(i)]));
                break;
            case Op::Log:
                add_adj(n.a, b.div(adj, primal[static_cast<std::size_t>(n.a)]));
                break;
            case Op::PowConst:
                add_adj(n.a, b.mul(adj, b.mul(b.constant(n.value), b.pow_const(primal[static_cast<std::size_t>(n.a)], n.value - 1.0))));
                break;
        }
        if (i == 0) {
            break;
        }
    }

    std::vector<Expr> out;
    out.reserve(static_cast<std::size_t>(input_group_size(f, wrt_input_name)));
    for (auto id : f.inputs) {
        const auto &n = f.graph.nodes[static_cast<std::size_t>(id)];
        if (n.op == Op::Input && n.name == wrt_input_name) {
            out.push_back(sum_terms(adj_terms[static_cast<std::size_t>(id)]));
        }
    }
    return out;
}

struct ReverseVectorBuilder {
    OptimizingBuilder &builder;
    std::vector<FunctionVectorNode> nodes;
    std::vector<std::vector<Expr>> lowered;

    int add_values(const std::vector<Expr> &values) {
        FunctionVectorNode node;
        node.op = VectorOp::Values;
        node.size = static_cast<int>(values.size());
        node.values.reserve(values.size());
        for (const auto &value : values) {
            node.values.push_back(value.id);
        }
        infer_direct_values_group(node, builder.g);
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int add_add(int lhs, int rhs) {
        return add_binary(VectorOp::Add, lhs, rhs);
    }

    int add_sub(int lhs, int rhs) {
        return add_binary(VectorOp::Sub, lhs, rhs);
    }

    int add_concat(int lhs, int rhs) {
        FunctionVectorNode node;
        node.op = VectorOp::Concat;
        node.size = nodes[static_cast<std::size_t>(lhs)].size + nodes[static_cast<std::size_t>(rhs)].size;
        node.a = lhs;
        node.b = rhs;
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int add_scale(double scale, int rhs) {
        FunctionVectorNode node;
        node.op = VectorOp::Scale;
        node.size = nodes[static_cast<std::size_t>(rhs)].size;
        node.a = rhs;
        node.scale = scale;
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int add_dense_matvec(DenseMatrix matrix, int rhs) {
        if (matrix.cols != nodes[static_cast<std::size_t>(rhs)].size) {
            throw std::runtime_error("reverse vector dense matvec dimension mismatch");
        }
        FunctionVectorNode node;
        node.op = VectorOp::DenseMatVec;
        node.size = matrix.rows;
        node.a = rhs;
        node.matrix = std::move(matrix);
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int add_sparse_matvec(SparseMatrix matrix, int rhs) {
        if (matrix.cols != nodes[static_cast<std::size_t>(rhs)].size) {
            throw std::runtime_error("reverse vector sparse matvec dimension mismatch");
        }
        FunctionVectorNode node;
        node.op = VectorOp::SparseMatVec;
        node.size = matrix.rows;
        node.a = rhs;
        node.sparse_matrix = std::move(matrix);
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int add_kron_eye_matvec(DenseMatrix matrix, int eye_size, int rhs) {
        if (matrix.cols * eye_size != nodes[static_cast<std::size_t>(rhs)].size) {
            throw std::runtime_error("reverse vector kron-eye matvec dimension mismatch");
        }
        FunctionVectorNode node;
        node.op = VectorOp::KronEyeMatVec;
        node.size = matrix.rows * eye_size;
        node.a = rhs;
        node.matrix = std::move(matrix);
        node.eye_size = eye_size;
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int add_scatter_slice(int source_size, int start, int stride, int adj) {
        auto adj_values = lower(adj);
        std::vector<Expr> values(static_cast<std::size_t>(source_size), builder.constant(0.0));
        for (int i = 0; i < static_cast<int>(adj_values.size()); ++i) {
            values[static_cast<std::size_t>(start + i * stride)] = adj_values[static_cast<std::size_t>(i)];
        }
        return add_values(values);
    }

    int add_slice(int source, int start, int length, int stride = 1) {
        FunctionVectorNode node;
        node.op = VectorOp::Slice;
        node.size = length;
        node.a = source;
        node.start = start;
        node.stride = stride;
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }

    int combine(const std::vector<int> &terms) {
        if (terms.empty()) {
            return -1;
        }
        int out = terms.front();
        for (std::size_t i = 1; i < terms.size(); ++i) {
            out = add_add(out, terms[i]);
        }
        return out;
    }

    std::vector<Expr> lower(int id) {
        auto &cached = lowered[static_cast<std::size_t>(id)];
        if (!cached.empty() || nodes[static_cast<std::size_t>(id)].size == 0) {
            return cached;
        }
        const auto &node = nodes[static_cast<std::size_t>(id)];
        switch (node.op) {
            case VectorOp::Values:
                cached.reserve(node.values.size());
                for (NodeId value : node.values) {
                    cached.emplace_back(&builder.g, value);
                }
                break;
            case VectorOp::Add:
                cached = vector_add(lower(node.a), lower(node.b));
                break;
            case VectorOp::Sub:
                cached = vector_sub(lower(node.a), lower(node.b));
                break;
            case VectorOp::Scale:
                cached = vector_scale(node.scale, lower(node.a));
                break;
            case VectorOp::DenseMatVec:
                cached = dense_matvec(node.matrix, lower(node.a));
                break;
            case VectorOp::SparseMatVec:
                cached = sparse_matvec(node.sparse_matrix, lower(node.a));
                break;
            case VectorOp::KronEyeMatVec:
                cached = kron_eye_matvec(node.matrix, node.eye_size, lower(node.a));
                break;
            case VectorOp::Slice: {
                auto base = lower(node.a);
                cached.reserve(static_cast<std::size_t>(node.size));
                for (int i = 0; i < node.size; ++i) {
                    cached.push_back(base[static_cast<std::size_t>(node.start + i * node.stride)]);
                }
                break;
            }
            case VectorOp::Concat: {
                auto lhs = lower(node.a);
                auto rhs = lower(node.b);
                cached.reserve(lhs.size() + rhs.size());
                cached.insert(cached.end(), lhs.begin(), lhs.end());
                cached.insert(cached.end(), rhs.begin(), rhs.end());
                break;
            }
        }
        return cached;
    }

private:
    int add_binary(VectorOp op, int lhs, int rhs) {
        if (nodes[static_cast<std::size_t>(lhs)].size != nodes[static_cast<std::size_t>(rhs)].size) {
            throw std::runtime_error("reverse vector binary dimension mismatch");
        }
        FunctionVectorNode node;
        node.op = op;
        node.size = nodes[static_cast<std::size_t>(lhs)].size;
        node.a = lhs;
        node.b = rhs;
        int id = static_cast<int>(nodes.size());
        nodes.push_back(std::move(node));
        lowered.emplace_back();
        return id;
    }
};

struct ReverseVectorOutput {
    int output = -1;
    std::vector<FunctionVectorNode> nodes;
};

ReverseVectorOutput reverse_vector_output(const GraphFunction &f,
                                          OptimizingBuilder &builder,
                                          const std::vector<Expr> &primal,
                                          const std::string &lambda_name,
                                          const std::string &wrt_input_name) {
    ReverseVectorOutput result;
    if (!f.has_vector_structure()) {
        return result;
    }
    if (f.output_vector < 0 || f.output_vector >= static_cast<int>(f.vector_nodes.size())) {
        return result;
    }
    if (f.vector_nodes[static_cast<std::size_t>(f.output_vector)].op == VectorOp::Values) {
        return result;
    }

    ReverseVectorBuilder vb{builder};
    std::vector<Expr> lambda_values;
    lambda_values.reserve(static_cast<std::size_t>(f.output_size()));
    for (int i = 0; i < f.output_size(); ++i) {
        lambda_values.push_back(builder.param(lambda_name, i));
    }

    std::vector<std::vector<int>> adj_terms(f.vector_nodes.size());
    adj_terms[static_cast<std::size_t>(f.output_vector)].push_back(vb.add_values(lambda_values));
    int wrt_output = -1;
    bool ok = true;

    for (int vector_id = static_cast<int>(f.vector_nodes.size()) - 1; vector_id >= 0; --vector_id) {
        if (!ok || adj_terms[static_cast<std::size_t>(vector_id)].empty()) {
            continue;
        }
        const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
        for (int term : adj_terms[static_cast<std::size_t>(vector_id)]) {
            if (vb.nodes[static_cast<std::size_t>(term)].size != node.size) {
                throw std::runtime_error("reverse vector adjoint size mismatch at vector node " + std::to_string(vector_id) + ": node size " + std::to_string(node.size) +
                                         ", adjoint size " + std::to_string(vb.nodes[static_cast<std::size_t>(term)].size));
            }
        }
        int adj = vb.combine(adj_terms[static_cast<std::size_t>(vector_id)]);
        switch (node.op) {
            case VectorOp::Values:
                if (direct_values_group(f, node, Op::Input, wrt_input_name) && node.size == input_group_size(f, wrt_input_name)) {
                    wrt_output = wrt_output < 0 ? adj : vb.add_add(wrt_output, adj);
                } else if (!node.is_input_group && !node.is_param_group) {
                    auto scalar_grad = scalar_vjp_for_values(f, builder, primal, node.values, vb.lower(adj), wrt_input_name);
                    if (static_cast<int>(scalar_grad.size()) != input_group_size(f, wrt_input_name)) {
                        ok = false;
                    } else {
                        int scalar_adj = vb.add_values(scalar_grad);
                        wrt_output = wrt_output < 0 ? scalar_adj : vb.add_add(wrt_output, scalar_adj);
                    }
                }
                break;
            case VectorOp::Add:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(adj);
                adj_terms[static_cast<std::size_t>(node.b)].push_back(adj);
                break;
            case VectorOp::Sub:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(adj);
                adj_terms[static_cast<std::size_t>(node.b)].push_back(vb.add_scale(-1.0, adj));
                break;
            case VectorOp::Scale:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(vb.add_scale(node.scale, adj));
                break;
            case VectorOp::DenseMatVec:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(vb.add_dense_matvec(transpose(node.matrix), adj));
                break;
            case VectorOp::SparseMatVec:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(vb.add_sparse_matvec(transpose(node.sparse_matrix), adj));
                break;
            case VectorOp::KronEyeMatVec:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(vb.add_kron_eye_matvec(transpose(node.matrix), node.eye_size, adj));
                break;
            case VectorOp::Slice:
                adj_terms[static_cast<std::size_t>(node.a)].push_back(vb.add_scatter_slice(f.vector_nodes[static_cast<std::size_t>(node.a)].size, node.start, node.stride, adj));
                break;
            case VectorOp::Concat: {
                const int lhs_size = f.vector_nodes[static_cast<std::size_t>(node.a)].size;
                const int rhs_size = f.vector_nodes[static_cast<std::size_t>(node.b)].size;
                adj_terms[static_cast<std::size_t>(node.a)].push_back(vb.add_slice(adj, 0, lhs_size));
                adj_terms[static_cast<std::size_t>(node.b)].push_back(vb.add_slice(adj, lhs_size, rhs_size));
                break;
            }
        }
    }

    if (!ok || wrt_output < 0 || vb.nodes[static_cast<std::size_t>(wrt_output)].size != input_group_size(f, wrt_input_name)) {
        return result;
    }

    result.output = wrt_output;
    result.nodes = std::move(vb.nodes);
    return result;
}

} // namespace

GraphFunction forward_diff(const GraphFunction &f, const std::string &wrt_input_name, const std::string &direction_name) {
    OptimizingBuilder b;
    auto primal = clone_nodes(f.graph, b);
    std::vector<Expr> tangent(f.graph.nodes.size());
    auto zero = b.constant(0.0);

    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        switch (n.op) {
            case Op::Constant:
            case Op::Param:
                tangent[i] = zero;
                break;
            case Op::Input:
                if (n.name == wrt_input_name) {
                    tangent[i] = b.param(direction_name, n.index);
                } else {
                    tangent[i] = zero;
                }
                break;
            case Op::Add:
                tangent[i] = b.add(tangent[n.a], tangent[n.b]);
                break;
            case Op::Sub:
                tangent[i] = b.sub(tangent[n.a], tangent[n.b]);
                break;
            case Op::Mul:
                tangent[i] = b.add(b.mul(tangent[n.a], primal[n.b]), b.mul(primal[n.a], tangent[n.b]));
                break;
            case Op::Div: {
                auto num = b.sub(b.mul(tangent[n.a], primal[n.b]), b.mul(primal[n.a], tangent[n.b]));
                auto den = b.mul(primal[n.b], primal[n.b]);
                tangent[i] = b.div(num, den);
                break;
            }
            case Op::Neg:
                tangent[i] = b.neg(tangent[n.a]);
                break;
            case Op::Sin:
                tangent[i] = b.mul(b.unary(Op::Cos, primal[n.a]), tangent[n.a]);
                break;
            case Op::Cos:
                tangent[i] = b.neg(b.mul(b.unary(Op::Sin, primal[n.a]), tangent[n.a]));
                break;
            case Op::Tan:
                tangent[i] = b.mul(tangent[n.a], b.add(b.constant(1.0), b.mul(primal[i], primal[i])));
                break;
            case Op::Exp:
                tangent[i] = b.mul(primal[i], tangent[n.a]);
                break;
            case Op::Log:
                tangent[i] = b.div(tangent[n.a], primal[n.a]);
                break;
            case Op::PowConst: {
                if (exactly(n.value, 0.0)) {
                    tangent[i] = zero;
                } else {
                    auto c = b.constant(n.value);
                    auto p = b.pow_const(primal[n.a], n.value - 1.0);
                    tangent[i] = b.mul(b.mul(c, p), tangent[n.a]);
                }
                break;
            }
        }
    }

    GraphFunction out;
    for (auto id : f.inputs) {
        out.inputs.push_back(primal[id].id);
    }

    for (auto id : f.params) {
        out.params.push_back(primal[id].id);
    }

    for (auto [name, size] : f.param_groups) {
        add_unique_group(out.param_groups, name, size);
    }

    for (auto [name, size] : f.input_groups) {
        add_unique_group(out.input_groups, name, size);
    }

    int wrt_size = 0;
    for (auto [name, size] : f.input_groups) {
        if (name == wrt_input_name) {
            wrt_size = size;
        }
    }
    add_unique_group(out.param_groups, direction_name, wrt_size);

    for (auto id : f.outputs) {
        out.outputs.push_back(tangent[id].id);
    }
    if (f.has_vector_structure() && f.output_vector >= 0 && f.output_vector < static_cast<int>(f.vector_nodes.size()) &&
        f.vector_nodes[static_cast<std::size_t>(f.output_vector)].op != VectorOp::Values) {
        out.vector_nodes = tangent_vector_nodes(f, b.g, tangent);
        out.input_vector = f.input_vector;
        out.output_vector = f.output_vector;
        out.param_vectors = f.param_vectors;
        out.vector_structure_valid = true;
    }
    out.graph = std::move(b.g);
    return optimize(out);
}

// -----------------------------------------------------------------------------
// Reverse-mode graph transform: VJP for vector-valued functions.
// Returns grad_x(lambda^T f(x)).
// -----------------------------------------------------------------------------

GraphFunction reverse_diff(const GraphFunction &f, const std::string &lambda_name, const std::string &wrt_input_name) {
    OptimizingBuilder b;
    auto primal = clone_nodes(f.graph, b);
    auto zero = b.constant(0.0);

    std::vector<std::vector<Expr>> adj_terms(f.graph.nodes.size());
    auto add_adj = [&](NodeId id, Expr term) {
        if (id >= 0) {
            adj_terms[id].push_back(term);
        }
    };

    for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i) {
        add_adj(f.outputs[i], b.param(lambda_name, i));
    }

    auto sum_terms = [&](const std::vector<Expr> &terms) -> Expr {
        if (terms.empty()) {
            return zero;
        }
        Expr s = terms[0];
        for (std::size_t i = 1; i < terms.size(); ++i) {
            s = b.add(s, terms[i]);
        }
        return s;
    };

    for (NodeId i = static_cast<NodeId>(f.graph.nodes.size()) - 1; i >= 0; --i) {
        if (adj_terms[i].empty()) {
            continue;
        }
        const auto &n = f.graph.nodes[i];
        Expr adj = sum_terms(adj_terms[i]);
        switch (n.op) {
            case Op::Constant:
            case Op::Input:
            case Op::Param:
                break;
            case Op::Add:
                add_adj(n.a, adj);
                add_adj(n.b, adj);
                break;
            case Op::Sub:
                add_adj(n.a, adj);
                add_adj(n.b, b.neg(adj));
                break;
            case Op::Mul:
                add_adj(n.a, b.mul(adj, primal[n.b]));
                add_adj(n.b, b.mul(adj, primal[n.a]));
                break;
            case Op::Div:
                add_adj(n.a, b.div(adj, primal[n.b]));
                add_adj(n.b, b.neg(b.div(b.mul(adj, primal[n.a]), b.mul(primal[n.b], primal[n.b]))));
                break;
            case Op::Neg:
                add_adj(n.a, b.neg(adj));
                break;
            case Op::Sin:
                add_adj(n.a, b.mul(adj, b.unary(Op::Cos, primal[n.a])));
                break;
            case Op::Cos:
                add_adj(n.a, b.neg(b.mul(adj, b.unary(Op::Sin, primal[n.a]))));
                break;
            case Op::Tan:
                add_adj(n.a, b.mul(adj, b.add(b.constant(1.0), b.mul(primal[i], primal[i]))));
                break;
            case Op::Exp:
                add_adj(n.a, b.mul(adj, primal[i]));
                break;
            case Op::Log:
                add_adj(n.a, b.div(adj, primal[n.a]));
                break;
            case Op::PowConst:
                add_adj(n.a, b.mul(adj, b.mul(b.constant(n.value), b.pow_const(primal[n.a], n.value - 1.0))));
                break;
        }
        if (i == 0) {
            break;
        }
    }

    GraphFunction out;

    int out_size = static_cast<int>(f.outputs.size());
    add_unique_group(out.param_groups, lambda_name, out_size);
    for (int i = 0; i < out_size; ++i) {
        // lambda seed nodes were created while seeding output adjoints.
        // They do not need to be duplicated here.
    }
    for (auto [name, size] : f.param_groups) {
        add_unique_group(out.param_groups, name, size);
    }
    for (auto [name, size] : f.input_groups) {
        add_unique_group(out.input_groups, name, size);
    }

    for (auto id : f.inputs) {
        out.inputs.push_back(primal[id].id);
    }
    for (auto id : f.params) {
        out.params.push_back(primal[id].id);
    }

    for (auto id : f.inputs) {
        const auto &n = f.graph.nodes[id];
        if (n.op == Op::Input && n.name == wrt_input_name) {
            out.outputs.push_back(sum_terms(adj_terms[id]).id);
        }
    }
    auto vector_output = reverse_vector_output(f, b, primal, lambda_name, wrt_input_name);
    if (vector_output.output >= 0) {
        out.vector_nodes = std::move(vector_output.nodes);
        out.output_vector = vector_output.output;
        out.vector_structure_valid = true;
    }
    out.graph = std::move(b.g);
    return optimize(out);
}


} // namespace ad
