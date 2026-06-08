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

#include "optimize.h"

namespace ad {

bool NodeKey::operator==(const NodeKey &o) const {
    return op == o.op && a == o.a && b == o.b && value == o.value && name == o.name && index == o.index;
}

std::size_t NodeKeyHash::operator()(const NodeKey &k) const {
    std::size_t h = std::hash<int>{}(static_cast<int>(k.op));
    auto mix = [&](std::size_t x) {
        h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    };
    mix(std::hash<int>{}(k.a));
    mix(std::hash<int>{}(k.b));
    mix(std::hash<double>{}(k.value));
    mix(std::hash<std::string>{}(k.name));
    mix(std::hash<int>{}(k.index));
    return h;
}

Expr OptimizingBuilder::make_raw(const Node &n) {
    NodeKey k{n.op, n.a, n.b, n.value, n.name, n.index};
    if (n.op == Op::Add || n.op == Op::Mul) {
        if (k.b < k.a) {
            std::swap(k.a, k.b);
        }
    }
    auto it = cse.find(k);
    if (it != cse.end()) {
        return Expr(&g, it->second);
    }
    Node nn = n;
    nn.a = k.a;
    nn.b = k.b;
    g.nodes.push_back(std::move(nn));
    NodeId id = static_cast<NodeId>(g.nodes.size() - 1);
    cse.emplace(std::move(k), id);
    return Expr(&g, id);
}

Expr OptimizingBuilder::constant(double v) {
    Node n;
    n.op = Op::Constant;
    n.value = v;
    return make_raw(n);
}

Expr OptimizingBuilder::input(const std::string &name, int index) {
    Node n;
    n.op = Op::Input;
    n.name = name;
    n.index = index;
    return make_raw(n);
}

Expr OptimizingBuilder::param(const std::string &name, int index) {
    Node n;
    n.op = Op::Param;
    n.name = name;
    n.index = index;
    return make_raw(n);
}

bool OptimizingBuilder::is_const(Expr e, double *value) const {
    const auto &n = g.nodes[e.id];
    if (n.op != Op::Constant) {
        return false;
    }
    if (value) {
        *value = n.value;
    }
    return true;
}

Expr OptimizingBuilder::add(Expr x, Expr y) {
    double cx, cy;
    if (is_const(x, &cx) && is_const(y, &cy)) {
        return constant(cx + cy);
    }
    if (is_const(x, &cx) && exactly(cx, 0.0)) {
        return y;
    }
    if (is_const(y, &cy) && exactly(cy, 0.0)) {
        return x;
    }
    if (x.id == y.id) {
        return mul(constant(2.0), x);
    }
    Node n;
    n.op = Op::Add;
    n.a = x.id;
    n.b = y.id;
    return make_raw(n);
}

Expr OptimizingBuilder::sub(Expr x, Expr y) {
    double cx, cy;
    if (is_const(x, &cx) && is_const(y, &cy)) {
        return constant(cx - cy);
    }
    if (is_const(y, &cy) && exactly(cy, 0.0)) {
        return x;
    }
    if (is_const(x, &cx) && exactly(cx, 0.0)) {
        return neg(y);
    }
    if (x.id == y.id) {
        return constant(0.0);
    }
    Node n;
    n.op = Op::Sub;
    n.a = x.id;
    n.b = y.id;
    return make_raw(n);
}

Expr OptimizingBuilder::mul(Expr x, Expr y) {
    double cx, cy;
    if (is_const(x, &cx) && is_const(y, &cy)) {
        return constant(cx * cy);
    }
    if (is_const(x, &cx)) {
        if (exactly(cx, 0.0)) {
            return constant(0.0);
        }
        if (exactly(cx, 1.0)) {
            return y;
        }
        if (exactly(cx, -1.0)) {
            return neg(y);
        }
    }
    if (is_const(y, &cy)) {
        if (exactly(cy, 0.0)) {
            return constant(0.0);
        }
        if (exactly(cy, 1.0)) {
            return x;
        }
        if (exactly(cy, -1.0)) {
            return neg(x);
        }
    }
    Node n;
    n.op = Op::Mul;
    n.a = x.id;
    n.b = y.id;
    return make_raw(n);
}

Expr OptimizingBuilder::div(Expr x, Expr y) {
    double cx, cy;
    if (is_const(x, &cx) && is_const(y, &cy)) {
        return constant(cx / cy);
    }
    if (is_const(x, &cx) && exactly(cx, 0.0)) {
        return constant(0.0);
    }
    if (is_const(y, &cy) && exactly(cy, 1.0)) {
        return x;
    }
    Node n;
    n.op = Op::Div;
    n.a = x.id;
    n.b = y.id;
    return make_raw(n);
}

Expr OptimizingBuilder::neg(Expr x) {
    double cx;
    if (is_const(x, &cx)) {
        return constant(-cx);
    }
    const auto &nx = g.nodes[x.id];
    if (nx.op == Op::Neg) {
        return Expr(&g, nx.a);
    }
    Node n;
    n.op = Op::Neg;
    n.a = x.id;
    return make_raw(n);
}

Expr OptimizingBuilder::unary(Op op, Expr x) {
    double cx;
    if (is_const(x, &cx)) {
        switch (op) {
            case Op::Sin:
                return constant(std::sin(cx));
            case Op::Cos:
                return constant(std::cos(cx));
            case Op::Tan:
                return constant(std::tan(cx));
            case Op::Exp:
                return constant(std::exp(cx));
            case Op::Log:
                return constant(std::log(cx));
            default:
                break;
        }
    }
    Node n;
    n.op = op;
    n.a = x.id;
    return make_raw(n);
}

Expr OptimizingBuilder::pow_const(Expr x, double p) {
    double cx;
    if (exactly(p, 0.0)) {
        return constant(1.0);
    }
    if (exactly(p, 1.0)) {
        return x;
    }
    if (exactly(p, 2.0)) {
        return mul(x, x);
    }
    if (is_const(x, &cx)) {
        return constant(std::pow(cx, p));
    }
    Node n;
    n.op = Op::PowConst;
    n.a = x.id;
    n.value = p;
    return make_raw(n);
}

Expr ob_binary(OptimizingBuilder &b, Op op, Expr x, Expr y) {
    switch (op) {
        case Op::Add:
            return b.add(x, y);
        case Op::Sub:
            return b.sub(x, y);
        case Op::Mul:
            return b.mul(x, y);
        case Op::Div:
            return b.div(x, y);
        default:
            throw std::runtime_error("not a binary op");
    }
}

std::vector<NodeId> topo_used(const Graph &g, const std::vector<NodeId> &outputs) {
    std::vector<char> seen(g.nodes.size(), false);
    std::vector<NodeId> order;
    std::function<void(NodeId)> dfs = [&](NodeId id) {
        if (id < 0 || seen[id]) {
            return;
        }
        seen[id] = true;
        const auto &n = g.nodes[id];
        dfs(n.a);
        dfs(n.b);
        order.push_back(id);
    };
    for (auto o : outputs) {
        dfs(o);
    }
    return order;
}

GraphFunction optimize(const GraphFunction &f) {
    std::vector<NodeId> roots = f.outputs;
    roots.insert(roots.end(), f.inputs.begin(), f.inputs.end());
    roots.insert(roots.end(), f.params.begin(), f.params.end());
    if (f.vector_structure_valid) {
        std::vector<char> reachable_vectors(f.vector_nodes.size(), false);
        std::function<void(int)> mark_vector = [&](int vector_id) {
            if (vector_id < 0 || vector_id >= static_cast<int>(f.vector_nodes.size())) {
                return;
            }
            if (reachable_vectors[static_cast<std::size_t>(vector_id)]) {
                return;
            }
            reachable_vectors[static_cast<std::size_t>(vector_id)] = true;
            const auto &vector_node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
            roots.insert(roots.end(), vector_node.values.begin(), vector_node.values.end());
            if (vector_node.op == VectorOp::Add || vector_node.op == VectorOp::Sub || vector_node.op == VectorOp::Concat) {
                mark_vector(vector_node.a);
                mark_vector(vector_node.b);
            } else if (vector_node.op == VectorOp::Scale || vector_node.op == VectorOp::DenseMatVec || vector_node.op == VectorOp::SparseMatVec ||
                       vector_node.op == VectorOp::KronEyeMatVec || vector_node.op == VectorOp::Slice) {
                mark_vector(vector_node.a);
            }
        };
        mark_vector(f.output_vector);
        for (int vector_id : f.param_vectors) {
            mark_vector(vector_id);
        }
        mark_vector(f.input_vector);
    }
    auto order = topo_used(f.graph, roots);
    OptimizingBuilder b;
    std::unordered_map<NodeId, Expr> map;

    for (NodeId old : order) {
        const auto &n = f.graph.nodes[old];
        Expr e;
        switch (n.op) {
            case Op::Constant:
                e = b.constant(n.value);
                break;
            case Op::Input:
                e = b.input(n.name, n.index);
                break;
            case Op::Param:
                e = b.param(n.name, n.index);
                break;
            case Op::Add:
                e = b.add(map.at(n.a), map.at(n.b));
                break;
            case Op::Sub:
                e = b.sub(map.at(n.a), map.at(n.b));
                break;
            case Op::Mul:
                e = b.mul(map.at(n.a), map.at(n.b));
                break;
            case Op::Div:
                e = b.div(map.at(n.a), map.at(n.b));
                break;
            case Op::Neg:
                e = b.neg(map.at(n.a));
                break;
            case Op::Sin:
                e = b.unary(Op::Sin, map.at(n.a));
                break;
            case Op::Cos:
                e = b.unary(Op::Cos, map.at(n.a));
                break;
            case Op::Tan:
                e = b.unary(Op::Tan, map.at(n.a));
                break;
            case Op::Exp:
                e = b.unary(Op::Exp, map.at(n.a));
                break;
            case Op::Log:
                e = b.unary(Op::Log, map.at(n.a));
                break;
            case Op::PowConst:
                e = b.pow_const(map.at(n.a), n.value);
                break;
        }
        map.emplace(old, e);
    }

    GraphFunction out;
    out.graph = std::move(b.g);
    for (auto id : f.inputs) {
        out.inputs.push_back(map.at(id).id);
    }
    for (auto id : f.params) {
        out.params.push_back(map.at(id).id);
    }
    for (auto id : f.outputs) {
        out.outputs.push_back(map.at(id).id);
    }
    out.input_groups = f.input_groups;
    out.param_groups = f.param_groups;
    out.vector_nodes = f.vector_nodes;
    out.input_vector = f.input_vector;
    out.output_vector = f.output_vector;
    out.param_vectors = f.param_vectors;
    out.vector_structure_valid = f.vector_structure_valid;
    if (out.vector_structure_valid) {
        std::vector<char> reachable_vectors(out.vector_nodes.size(), false);
        std::function<void(int)> mark_vector = [&](int vector_id) {
            if (vector_id < 0 || vector_id >= static_cast<int>(out.vector_nodes.size())) {
                out.vector_structure_valid = false;
                return;
            }
            if (reachable_vectors[static_cast<std::size_t>(vector_id)]) {
                return;
            }
            reachable_vectors[static_cast<std::size_t>(vector_id)] = true;
            const auto &vector_node = out.vector_nodes[static_cast<std::size_t>(vector_id)];
            if (vector_node.op == VectorOp::Add || vector_node.op == VectorOp::Sub || vector_node.op == VectorOp::Concat) {
                mark_vector(vector_node.a);
                mark_vector(vector_node.b);
            } else if (vector_node.op == VectorOp::Scale || vector_node.op == VectorOp::DenseMatVec || vector_node.op == VectorOp::SparseMatVec ||
                       vector_node.op == VectorOp::KronEyeMatVec || vector_node.op == VectorOp::Slice) {
                mark_vector(vector_node.a);
            }
        };
        mark_vector(out.output_vector);

        for (int vector_id = 0; vector_id < static_cast<int>(out.vector_nodes.size()); ++vector_id) {
            if (!reachable_vectors[static_cast<std::size_t>(vector_id)]) {
                auto &vector_node = out.vector_nodes[static_cast<std::size_t>(vector_id)];
                vector_node = FunctionVectorNode{};
                continue;
            }
            auto &vector_node = out.vector_nodes[static_cast<std::size_t>(vector_id)];
            for (auto &id : vector_node.values) {
                auto it = map.find(id);
                if (it == map.end()) {
                    out.vector_structure_valid = false;
                    break;
                }
                id = it->second.id;
            }
            if (!out.vector_structure_valid) {
                break;
            }
        }
    }
    return out;
}

} // namespace ad
