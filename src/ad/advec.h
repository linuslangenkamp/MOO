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

#ifndef MOO_ADVEC_H
#define MOO_ADVEC_H

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

namespace advec
{

// -----------------------------------------------------------------------------
// Header-only symbolic graph AD for vector-valued functions.
//
// Design summary:
//   - Build vector functions as DAGs.
//   - reverse_diff(F, lambda) builds graph for grad_x(lambda^T F(x)).
//   - forward_diff(G, wrt, v) builds graph for directional derivative of G.
//   - VM evaluates optimized graphs.
//   - CEmitter emits one-shot C and staged C. Staged C separates nodes that do
//     not depend on direction variables from nodes that do depend on them.
// -----------------------------------------------------------------------------

struct Graph;

using NodeId = int;

inline constexpr NodeId invalid_node = -1;

inline bool nearly_zero(double x)
{
    return std::abs(x) <= 0.0;
}

inline bool exactly(double a, double b)
{
    return a == b;
}

enum class Op
{
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

inline const char *op_name(Op op)
{
    switch (op)
    {
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

struct Node
{
    Op op = Op::Constant;
    NodeId a = invalid_node;
    NodeId b = invalid_node;
    double value = 0.0; // constants and pow exponent
    std::string name;   // input/param group name
    int index = -1;     // input/param index in group
};

struct Expr
{
    Graph *g = nullptr;
    NodeId id = invalid_node;

    Expr() = default;
    Expr(Graph *graph, NodeId node)
        : g(graph),
          id(node)
    {
    }

    explicit operator bool() const { return g != nullptr && id != invalid_node; }
};

struct NamedVector
{
    std::string name;
    std::vector<Expr> values;

    std::size_t size() const { return values.size(); }
    Expr &operator[](std::size_t i) { return values[i]; }
    const Expr &operator[](std::size_t i) const { return values[i]; }
};

struct Graph
{
    std::vector<Node> nodes;

    Expr constant(double v)
    {
        Node n;
        n.op = Op::Constant;
        n.value = v;
        nodes.push_back(std::move(n));
        return Expr(this, static_cast<NodeId>(nodes.size() - 1));
    }

    Expr input(const std::string &name, int index)
    {
        Node n;
        n.op = Op::Input;
        n.name = name;
        n.index = index;
        nodes.push_back(std::move(n));
        return Expr(this, static_cast<NodeId>(nodes.size() - 1));
    }

    Expr param(const std::string &name, int index)
    {
        Node n;
        n.op = Op::Param;
        n.name = name;
        n.index = index;
        nodes.push_back(std::move(n));
        return Expr(this, static_cast<NodeId>(nodes.size() - 1));
    }

    Expr unary(Op op, Expr x)
    {
        check_same(x);
        Node n;
        n.op = op;
        n.a = x.id;
        nodes.push_back(std::move(n));
        return Expr(this, static_cast<NodeId>(nodes.size() - 1));
    }

    Expr binary(Op op, Expr x, Expr y)
    {
        check_same(x);
        check_same(y);
        Node n;
        n.op = op;
        n.a = x.id;
        n.b = y.id;
        nodes.push_back(std::move(n));
        return Expr(this, static_cast<NodeId>(nodes.size() - 1));
    }

    Expr pow_const(Expr x, double p)
    {
        check_same(x);
        Node n;
        n.op = Op::PowConst;
        n.a = x.id;
        n.value = p;
        nodes.push_back(std::move(n));
        return Expr(this, static_cast<NodeId>(nodes.size() - 1));
    }

    void check_same(Expr e) const
    {
        if (e.g != this || e.id < 0 || e.id >= static_cast<NodeId>(nodes.size()))
        {
            throw std::runtime_error("expression belongs to a different graph or is invalid");
        }
    }

    NamedVector inputs(const std::string &name, int n)
    {
        NamedVector v;
        v.name = name;
        for (int i = 0; i < n; ++i)
        {
            v.values.push_back(input(name, i));
        }
        return v;
    }

    NamedVector params(const std::string &name, int n)
    {
        NamedVector v;
        v.name = name;
        for (int i = 0; i < n; ++i)
        {
            v.values.push_back(param(name, i));
        }
        return v;
    }
};

inline void require_same_graph(Expr a, Expr b)
{
    if (a.g != b.g)
    {
        throw std::runtime_error("expressions belong to different graphs");
    }
}

inline Expr operator+(Expr a, Expr b)
{
    require_same_graph(a, b);
    return a.g->binary(Op::Add, a, b);
}

inline Expr operator-(Expr a, Expr b)
{
    require_same_graph(a, b);
    return a.g->binary(Op::Sub, a, b);
}

inline Expr operator*(Expr a, Expr b)
{
    require_same_graph(a, b);
    return a.g->binary(Op::Mul, a, b);
}

inline Expr operator/(Expr a, Expr b)
{
    require_same_graph(a, b);
    return a.g->binary(Op::Div, a, b);
}

inline Expr operator-(Expr a)
{
    return a.g->unary(Op::Neg, a);
}

inline Expr operator+(Expr a, double b)
{
    return a + a.g->constant(b);
}

inline Expr operator-(Expr a, double b)
{
    return a - a.g->constant(b);
}

inline Expr operator*(Expr a, double b)
{
    return a * a.g->constant(b);
}

inline Expr operator/(Expr a, double b)
{
    return a / a.g->constant(b);
}

inline Expr operator+(double a, Expr b)
{
    return b.g->constant(a) + b;
}

inline Expr operator-(double a, Expr b)
{
    return b.g->constant(a) - b;
}

inline Expr operator*(double a, Expr b)
{
    return b.g->constant(a) * b;
}

inline Expr operator/(double a, Expr b)
{
    return b.g->constant(a) / b;
}

inline Expr sin(Expr x)
{
    return x.g->unary(Op::Sin, x);
}

inline Expr cos(Expr x)
{
    return x.g->unary(Op::Cos, x);
}

inline Expr tan(Expr x)
{
    return x.g->unary(Op::Tan, x);
}

inline Expr exp(Expr x)
{
    return x.g->unary(Op::Exp, x);
}

inline Expr log(Expr x)
{
    return x.g->unary(Op::Log, x);
}

inline Expr pow_const(Expr x, double p)
{
    return x.g->pow_const(x, p);
}

inline Expr sqr(Expr x)
{
    return x * x;
}

struct GraphFunction
{
    Graph graph;
    std::vector<NodeId> inputs;
    std::vector<NodeId> params;
    std::vector<NodeId> outputs;

    std::vector<std::pair<std::string, int>> input_groups;
    std::vector<std::pair<std::string, int>> param_groups;

    int input_size() const { return static_cast<int>(inputs.size()); }
    int param_size() const { return static_cast<int>(params.size()); }
    int output_size() const { return static_cast<int>(outputs.size()); }
};

inline void add_unique_group(std::vector<std::pair<std::string, int>> &groups, const std::string &name, int size)
{
    for (auto &g : groups)
    {
        if (g.first == name)
        {
            g.second = std::max(g.second, size);
            return;
        }
    }
    groups.push_back({name, size});
}

inline GraphFunction function_from(Graph &&graph, const NamedVector &inputs, const std::vector<Expr> &outputs, const std::vector<NamedVector> &params = {})
{
    GraphFunction f;
    f.graph = std::move(graph);
    for (auto e : inputs.values)
    {
        f.inputs.push_back(e.id);
    }
    for (auto e : outputs)
    {
        f.outputs.push_back(e.id);
    }
    add_unique_group(f.input_groups, inputs.name, static_cast<int>(inputs.size()));
    for (const auto &p : params)
    {
        add_unique_group(f.param_groups, p.name, static_cast<int>(p.size()));
        for (auto e : p.values)
        {
            f.params.push_back(e.id);
        }
    }
    return f;
}

// -----------------------------------------------------------------------------
// Graph cloning / optimized building helpers
// -----------------------------------------------------------------------------

struct NodeKey
{
    Op op{};
    NodeId a = invalid_node;
    NodeId b = invalid_node;
    double value = 0.0;
    std::string name;
    int index = -1;

    bool operator==(const NodeKey &o) const { return op == o.op && a == o.a && b == o.b && value == o.value && name == o.name && index == o.index; }
};

struct NodeKeyHash
{
    std::size_t operator()(const NodeKey &k) const
    {
        std::size_t h = std::hash<int>{}(static_cast<int>(k.op));
        auto mix = [&](std::size_t x)
        {
            h ^= x + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        };
        mix(std::hash<int>{}(k.a));
        mix(std::hash<int>{}(k.b));
        mix(std::hash<double>{}(k.value));
        mix(std::hash<std::string>{}(k.name));
        mix(std::hash<int>{}(k.index));
        return h;
    }
};

struct OptimizingBuilder
{
    Graph g;
    std::unordered_map<NodeKey, NodeId, NodeKeyHash> cse;

    Expr make_raw(const Node &n)
    {
        NodeKey k{n.op, n.a, n.b, n.value, n.name, n.index};
        if (n.op == Op::Add || n.op == Op::Mul)
        {
            if (k.b < k.a)
            {
                std::swap(k.a, k.b);
            }
        }
        auto it = cse.find(k);
        if (it != cse.end())
        {
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

    Expr constant(double v)
    {
        Node n;
        n.op = Op::Constant;
        n.value = v;
        return make_raw(n);
    }

    Expr input(const std::string &name, int index)
    {
        Node n;
        n.op = Op::Input;
        n.name = name;
        n.index = index;
        return make_raw(n);
    }

    Expr param(const std::string &name, int index)
    {
        Node n;
        n.op = Op::Param;
        n.name = name;
        n.index = index;
        return make_raw(n);
    }

    bool is_const(Expr e, double *value = nullptr) const
    {
        const auto &n = g.nodes[e.id];
        if (n.op != Op::Constant)
        {
            return false;
        }
        if (value)
        {
            *value = n.value;
        }
        return true;
    }

    Expr add(Expr x, Expr y)
    {
        double cx, cy;
        if (is_const(x, &cx) && is_const(y, &cy))
        {
            return constant(cx + cy);
        }
        if (is_const(x, &cx) && exactly(cx, 0.0))
        {
            return y;
        }
        if (is_const(y, &cy) && exactly(cy, 0.0))
        {
            return x;
        }
        if (x.id == y.id)
        {
            return mul(constant(2.0), x);
        }
        Node n;
        n.op = Op::Add;
        n.a = x.id;
        n.b = y.id;
        return make_raw(n);
    }

    Expr sub(Expr x, Expr y)
    {
        double cx, cy;
        if (is_const(x, &cx) && is_const(y, &cy))
        {
            return constant(cx - cy);
        }
        if (is_const(y, &cy) && exactly(cy, 0.0))
        {
            return x;
        }
        if (is_const(x, &cx) && exactly(cx, 0.0))
        {
            return neg(y);
        }
        if (x.id == y.id)
        {
            return constant(0.0);
        }
        Node n;
        n.op = Op::Sub;
        n.a = x.id;
        n.b = y.id;
        return make_raw(n);
    }

    Expr mul(Expr x, Expr y)
    {
        double cx, cy;
        if (is_const(x, &cx) && is_const(y, &cy))
        {
            return constant(cx * cy);
        }
        if (is_const(x, &cx))
        {
            if (exactly(cx, 0.0))
            {
                return constant(0.0);
            }
            if (exactly(cx, 1.0))
            {
                return y;
            }
            if (exactly(cx, -1.0))
            {
                return neg(y);
            }
        }
        if (is_const(y, &cy))
        {
            if (exactly(cy, 0.0))
            {
                return constant(0.0);
            }
            if (exactly(cy, 1.0))
            {
                return x;
            }
            if (exactly(cy, -1.0))
            {
                return neg(x);
            }
        }
        Node n;
        n.op = Op::Mul;
        n.a = x.id;
        n.b = y.id;
        return make_raw(n);
    }

    Expr div(Expr x, Expr y)
    {
        double cx, cy;
        if (is_const(x, &cx) && is_const(y, &cy))
        {
            return constant(cx / cy);
        }
        if (is_const(x, &cx) && exactly(cx, 0.0))
        {
            return constant(0.0);
        }
        if (is_const(y, &cy) && exactly(cy, 1.0))
        {
            return x;
        }
        Node n;
        n.op = Op::Div;
        n.a = x.id;
        n.b = y.id;
        return make_raw(n);
    }

    Expr neg(Expr x)
    {
        double cx;
        if (is_const(x, &cx))
        {
            return constant(-cx);
        }
        const auto &nx = g.nodes[x.id];
        if (nx.op == Op::Neg)
        {
            return Expr(&g, nx.a);
        }
        Node n;
        n.op = Op::Neg;
        n.a = x.id;
        return make_raw(n);
    }

    Expr unary(Op op, Expr x)
    {
        double cx;
        if (is_const(x, &cx))
        {
            switch (op)
            {
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

    Expr pow_const(Expr x, double p)
    {
        double cx;
        if (exactly(p, 0.0))
        {
            return constant(1.0);
        }
        if (exactly(p, 1.0))
        {
            return x;
        }
        if (exactly(p, 2.0))
        {
            return mul(x, x);
        }
        if (is_const(x, &cx))
        {
            return constant(std::pow(cx, p));
        }
        Node n;
        n.op = Op::PowConst;
        n.a = x.id;
        n.value = p;
        return make_raw(n);
    }
};

inline Expr ob_binary(OptimizingBuilder &b, Op op, Expr x, Expr y)
{
    switch (op)
    {
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

inline std::vector<NodeId> topo_used(const Graph &g, const std::vector<NodeId> &outputs)
{
    std::vector<char> seen(g.nodes.size(), false);
    std::vector<NodeId> order;
    std::function<void(NodeId)> dfs = [&](NodeId id)
    {
        if (id < 0 || seen[id])
        {
            return;
        }
        seen[id] = true;
        const auto &n = g.nodes[id];
        dfs(n.a);
        dfs(n.b);
        order.push_back(id);
    };
    for (auto o : outputs)
    {
        dfs(o);
    }
    return order;
}

inline GraphFunction optimize(const GraphFunction &f)
{
    auto order = topo_used(f.graph, f.outputs);
    OptimizingBuilder b;
    std::unordered_map<NodeId, Expr> map;

    for (NodeId old : order)
    {
        const auto &n = f.graph.nodes[old];
        Expr e;
        switch (n.op)
        {
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
    for (auto id : f.inputs)
    {
        if (map.count(id))
        {
            out.inputs.push_back(map.at(id).id);
        }
    }
    for (auto id : f.params)
    {
        if (map.count(id))
        {
            out.params.push_back(map.at(id).id);
        }
    }
    for (auto id : f.outputs)
    {
        out.outputs.push_back(map.at(id).id);
    }
    out.input_groups = f.input_groups;
    out.param_groups = f.param_groups;
    return out;
}

// -----------------------------------------------------------------------------
// Copy graph into OptimizingBuilder
// -----------------------------------------------------------------------------

inline std::vector<Expr> clone_nodes(const Graph &src, OptimizingBuilder &dst, const std::optional<std::string> &input_as_param = std::nullopt)
{
    std::vector<Expr> map(src.nodes.size());
    for (NodeId i = 0; i < static_cast<NodeId>(src.nodes.size()); ++i)
    {
        const auto &n = src.nodes[i];
        switch (n.op)
        {
            case Op::Constant:
                map[i] = dst.constant(n.value);
                break;
            case Op::Input:
                if (input_as_param && n.name == *input_as_param)
                {
                    map[i] = dst.param(n.name, n.index);
                }
                else
                {
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

inline GraphFunction forward_diff(const GraphFunction &f, const std::string &wrt_input_name = "x", const std::string &direction_name = "v")
{
    OptimizingBuilder b;
    auto primal = clone_nodes(f.graph, b);
    std::vector<Expr> tangent(f.graph.nodes.size());
    auto zero = b.constant(0.0);

    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
    {
        const auto &n = f.graph.nodes[i];
        switch (n.op)
        {
            case Op::Constant:
            case Op::Param:
                tangent[i] = zero;
                break;
            case Op::Input:
                if (n.name == wrt_input_name)
                {
                    tangent[i] = b.param(direction_name, n.index);
                }
                else
                {
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
            case Op::Div:
            {
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
            case Op::PowConst:
            {
                if (exactly(n.value, 0.0))
                {
                    tangent[i] = zero;
                }
                else
                {
                    auto c = b.constant(n.value);
                    auto p = b.pow_const(primal[n.a], n.value - 1.0);
                    tangent[i] = b.mul(b.mul(c, p), tangent[n.a]);
                }
                break;
            }
        }
    }

    GraphFunction out;
    for (auto id : f.inputs)
    {
        out.inputs.push_back(primal[id].id);
    }

    for (auto id : f.params)
    {
        out.params.push_back(primal[id].id);
    }

    for (auto [name, size] : f.param_groups)
    {
        add_unique_group(out.param_groups, name, size);
    }

    for (auto [name, size] : f.input_groups)
    {
        add_unique_group(out.input_groups, name, size);
    }

    int wrt_size = 0;
    for (auto [name, size] : f.input_groups)
    {
        if (name == wrt_input_name)
        {
            wrt_size = size;
        }
    }
    add_unique_group(out.param_groups, direction_name, wrt_size);

    for (auto id : f.outputs)
    {
        out.outputs.push_back(tangent[id].id);
    }
    out.graph = std::move(b.g);
    return optimize(out);
}

// -----------------------------------------------------------------------------
// Reverse-mode graph transform: VJP for vector-valued functions.
// Returns grad_x(lambda^T f(x)).
// -----------------------------------------------------------------------------

inline GraphFunction reverse_diff(const GraphFunction &f, const std::string &lambda_name = "lambda", const std::string &wrt_input_name = "x")
{
    OptimizingBuilder b;
    auto primal = clone_nodes(f.graph, b);
    auto zero = b.constant(0.0);

    std::vector<std::vector<Expr>> adj_terms(f.graph.nodes.size());
    auto add_adj = [&](NodeId id, Expr term)
    {
        if (id >= 0)
        {
            adj_terms[id].push_back(term);
        }
    };

    for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i)
    {
        add_adj(f.outputs[i], b.param(lambda_name, i));
    }

    auto sum_terms = [&](const std::vector<Expr> &terms) -> Expr
    {
        if (terms.empty())
        {
            return zero;
        }
        Expr s = terms[0];
        for (std::size_t i = 1; i < terms.size(); ++i)
        {
            s = b.add(s, terms[i]);
        }
        return s;
    };

    for (NodeId i = static_cast<NodeId>(f.graph.nodes.size()) - 1; i >= 0; --i)
    {
        if (adj_terms[i].empty())
        {
            continue;
        }
        const auto &n = f.graph.nodes[i];
        Expr adj = sum_terms(adj_terms[i]);
        switch (n.op)
        {
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
        if (i == 0)
        {
            break;
        }
    }

    GraphFunction out;

    int out_size = static_cast<int>(f.outputs.size());
    add_unique_group(out.param_groups, lambda_name, out_size);
    for (int i = 0; i < out_size; ++i)
    {
        // lambda seed nodes were created while seeding output adjoints.
        // They do not need to be duplicated here.
    }
    for (auto [name, size] : f.param_groups)
    {
        add_unique_group(out.param_groups, name, size);
    }
    for (auto [name, size] : f.input_groups)
    {
        add_unique_group(out.input_groups, name, size);
    }

    for (auto id : f.inputs)
    {
        out.inputs.push_back(primal[id].id);
    }
    for (auto id : f.params)
    {
        out.params.push_back(primal[id].id);
    }

    for (auto id : f.inputs)
    {
        const auto &n = f.graph.nodes[id];
        if (n.op == Op::Input && n.name == wrt_input_name)
        {
            out.outputs.push_back(sum_terms(adj_terms[id]).id);
        }
    }
    out.graph = std::move(b.g);
    return optimize(out);
}

// -----------------------------------------------------------------------------
// Evaluation environment and VM.
// -----------------------------------------------------------------------------

struct EvalEnv
{
    std::unordered_map<std::string, const double *> inputs;
    std::unordered_map<std::string, const double *> params;

    EvalEnv &input(const std::string &name, const double *data)
    {
        inputs[name] = data;
        return *this;
    }
    EvalEnv &param(const std::string &name, const double *data)
    {
        params[name] = data;
        return *this;
    }
};

struct VM
{
    GraphFunction f;

    explicit VM(GraphFunction fn)
        : f(optimize(fn))
    {
    }

    void evaluate(const EvalEnv &env, double *out) const
    {
        std::vector<double> mem(f.graph.nodes.size(), 0.0);
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            switch (n.op)
            {
                case Op::Constant:
                    mem[i] = n.value;
                    break;
                case Op::Input:
                {
                    auto it = env.inputs.find(n.name);
                    if (it == env.inputs.end())
                    {
                        throw std::runtime_error("missing input group: " + n.name);
                    }
                    mem[i] = it->second[n.index];
                    break;
                }
                case Op::Param:
                {
                    auto it = env.params.find(n.name);
                    if (it == env.params.end())
                    {
                        throw std::runtime_error("missing param group: " + n.name);
                    }
                    mem[i] = it->second[n.index];
                    break;
                }
                case Op::Add:
                    mem[i] = mem[n.a] + mem[n.b];
                    break;
                case Op::Sub:
                    mem[i] = mem[n.a] - mem[n.b];
                    break;
                case Op::Mul:
                    mem[i] = mem[n.a] * mem[n.b];
                    break;
                case Op::Div:
                    mem[i] = mem[n.a] / mem[n.b];
                    break;
                case Op::Neg:
                    mem[i] = -mem[n.a];
                    break;
                case Op::Sin:
                    mem[i] = std::sin(mem[n.a]);
                    break;
                case Op::Cos:
                    mem[i] = std::cos(mem[n.a]);
                    break;
                case Op::Tan:
                    mem[i] = std::tan(mem[n.a]);
                    break;
                case Op::Exp:
                    mem[i] = std::exp(mem[n.a]);
                    break;
                case Op::Log:
                    mem[i] = std::log(mem[n.a]);
                    break;
                case Op::PowConst:
                    mem[i] = std::pow(mem[n.a], n.value);
                    break;
            }
        }
        for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i)
        {
            out[i] = mem[f.outputs[i]];
        }
    }
};

// -----------------------------------------------------------------------------
// Staged VM: prepares all nodes independent of a chosen direction parameter group.
// -----------------------------------------------------------------------------

struct StagedVM
{
    GraphFunction f;
    std::string direction_name;
    std::vector<char> depends_on_direction;

    explicit StagedVM(GraphFunction fn, std::string direction = "v")
        : f(optimize(fn)),
          direction_name(std::move(direction))
    {
        analyze();
    }

    void analyze()
    {
        depends_on_direction.assign(f.graph.nodes.size(), false);
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            bool dep = false;
            if (n.op == Op::Param && n.name == direction_name)
            {
                dep = true;
            }
            if (n.a >= 0)
            {
                dep = dep || depends_on_direction[n.a];
            }
            if (n.b >= 0)
            {
                dep = dep || depends_on_direction[n.b];
            }
            depends_on_direction[i] = dep;
        }
    }

    struct Prepared
    {
        const StagedVM *vm = nullptr;
        std::vector<double> mem;

        void apply(const double *direction, double *out) const
        {
            EvalEnv env;
            env.params.emplace(vm->direction_name, direction);
            std::vector<double> local = mem;
            vm->eval_dependent(env, local, out);
        }
    };

    Prepared prepare(const EvalEnv &env) const
    {
        Prepared p;
        p.vm = this;
        p.mem.assign(f.graph.nodes.size(), 0.0);
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            if (depends_on_direction[i])
            {
                continue;
            }
            eval_node(i, env, p.mem);
        }
        return p;
    }

    void eval_dependent(const EvalEnv &env, std::vector<double> &mem, double *out) const
    {
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            if (!depends_on_direction[i])
            {
                continue;
            }
            eval_node(i, env, mem);
        }
        for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i)
        {
            out[i] = mem[f.outputs[i]];
        }
    }

    void eval_node(NodeId i, const EvalEnv &env, std::vector<double> &mem) const
    {
        const auto &n = f.graph.nodes[i];
        switch (n.op)
        {
            case Op::Constant:
                mem[i] = n.value;
                break;
            case Op::Input:
            {
                auto it = env.inputs.find(n.name);
                if (it == env.inputs.end())
                {
                    throw std::runtime_error("missing input group: " + n.name);
                }
                mem[i] = it->second[n.index];
                break;
            }
            case Op::Param:
            {
                auto it = env.params.find(n.name);
                if (it == env.params.end())
                {
                    throw std::runtime_error("missing param group: " + n.name);
                }
                mem[i] = it->second[n.index];
                break;
            }
            case Op::Add:
                mem[i] = mem[n.a] + mem[n.b];
                break;
            case Op::Sub:
                mem[i] = mem[n.a] - mem[n.b];
                break;
            case Op::Mul:
                mem[i] = mem[n.a] * mem[n.b];
                break;
            case Op::Div:
                mem[i] = mem[n.a] / mem[n.b];
                break;
            case Op::Neg:
                mem[i] = -mem[n.a];
                break;
            case Op::Sin:
                mem[i] = std::sin(mem[n.a]);
                break;
            case Op::Cos:
                mem[i] = std::cos(mem[n.a]);
                break;
            case Op::Tan:
                mem[i] = std::tan(mem[n.a]);
                break;
            case Op::Exp:
                mem[i] = std::exp(mem[n.a]);
                break;
            case Op::Log:
                mem[i] = std::log(mem[n.a]);
                break;
            case Op::PowConst:
                mem[i] = std::pow(mem[n.a], n.value);
                break;
        }
    }
};

// -----------------------------------------------------------------------------
// C emitter.
// -----------------------------------------------------------------------------

struct CEmitter
{
    static std::string number(double v)
    {
        std::ostringstream os;
        os << std::setprecision(17) << v;
        return os.str();
    }

    static std::string cname(const std::string &s)
    {
        std::string r;
        for (char c : s)
        {
            if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_')
            {
                r.push_back(c);
            }
            else
            {
                r.push_back('_');
            }
        }
        if (r.empty() || (r[0] >= '0' && r[0] <= '9'))
        {
            r = "_" + r;
        }
        return r;
    }

    static std::string cache_slot(const Node &n) { return cname(n.name) + "_" + std::to_string(n.index); }

    static std::string node_ref(const Node &n)
    {
        if (n.op == Op::Input)
        {
            return n.name + "[" + std::to_string(n.index) + "]";
        }
        if (n.op == Op::Param)
        {
            return n.name + "[" + std::to_string(n.index) + "]";
        }
        throw std::runtime_error("node_ref only supports input/param");
    }

    static std::string expr_rhs(const Graph &g, NodeId id, const std::function<std::string(NodeId)> &ref)
    {
        const auto &n = g.nodes[id];
        switch (n.op)
        {
            case Op::Constant:
                return number(n.value);
            case Op::Input:
                return node_ref(n);
            case Op::Param:
                return node_ref(n);
            case Op::Add:
                return ref(n.a) + " + " + ref(n.b);
            case Op::Sub:
                return ref(n.a) + " - " + ref(n.b);
            case Op::Mul:
                return ref(n.a) + " * " + ref(n.b);
            case Op::Div:
                return ref(n.a) + " / " + ref(n.b);
            case Op::Neg:
                return "-" + ref(n.a);
            case Op::Sin:
                return "sin(" + ref(n.a) + ")";
            case Op::Cos:
                return "cos(" + ref(n.a) + ")";
            case Op::Tan:
                return "tan(" + ref(n.a) + ")";
            case Op::Exp:
                return "exp(" + ref(n.a) + ")";
            case Op::Log:
                return "log(" + ref(n.a) + ")";
            case Op::PowConst:
                return "pow(" + ref(n.a) + ", " + number(n.value) + ")";
        }
        return "0.0";
    }

    static void emit_function(const GraphFunction &fn, const std::string &name, std::ostream &os)
    {
        auto f = optimize(fn);
        os << "#include <math.h>\n\n";
        os << "void " << name << "(\n";
        for (auto [n, s] : f.input_groups)
        {
            os << "    const double* " << n << ",\n";
        }
        for (auto [n, s] : f.param_groups)
        {
            os << "    const double* " << n << ",\n";
        }
        os << "    double* out\n) {\n";

        std::unordered_set<NodeId> output_set(f.outputs.begin(), f.outputs.end());
        auto ref = [&](NodeId id) -> std::string
        {
            const auto &n = f.graph.nodes[id];
            if (n.op == Op::Input || n.op == Op::Param)
            {
                return node_ref(n);
            }
            if (n.op == Op::Constant)
            {
                return number(n.value);
            }
            return "t" + std::to_string(id);
        };

        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            if (n.op == Op::Input || n.op == Op::Param || n.op == Op::Constant)
            {
                continue;
            }
            os << "    double t" << i << " = " << expr_rhs(f.graph, i, ref) << ";\n";
        }
        for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i)
        {
            os << "    out[" << i << "] = " << ref(f.outputs[i]) << ";\n";
        }
        os << "}\n";
    }

    static void emit_staged(const GraphFunction &fn, const std::string &basename, const std::string &direction_name, std::ostream &os)
    {
        StagedVM svm(fn, direction_name);
        const auto &f = svm.f;
        os << "#include <math.h>\n\n";
        os << "typedef struct {\n";
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            if ((n.op == Op::Input) || (n.op == Op::Param && n.name != direction_name))
            {
                os << "    double " << cache_slot(n) << ";\n";
            }
        }
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            if (!svm.depends_on_direction[i] && n.op != Op::Input && n.op != Op::Param && n.op != Op::Constant)
            {
                os << "    double t" << i << ";\n";
            }
        }
        os << "} " << basename << "_cache_t;\n\n";

        auto ref_prepare = [&](NodeId id) -> std::string
        {
            const auto &n = f.graph.nodes[id];
            if (n.op == Op::Input || n.op == Op::Param)
            {
                return node_ref(n);
            }
            if (n.op == Op::Constant)
            {
                return number(n.value);
            }
            return "cache->t" + std::to_string(id);
        };

        auto ref_apply = [&](NodeId id) -> std::string
        {
            const auto &n = f.graph.nodes[id];
            if (n.op == Op::Input)
            {
                return "cache->" + cache_slot(n);
            }
            if (n.op == Op::Param)
            {
                if (n.name == direction_name)
                {
                    return node_ref(n);
                }
                return "cache->" + cache_slot(n);
            }
            if (n.op == Op::Constant)
            {
                return number(n.value);
            }
            if (!svm.depends_on_direction[id])
            {
                return "cache->t" + std::to_string(id);
            }
            return "t" + std::to_string(id);
        };

        os << "void " << basename << "_prepare(\n";
        for (auto [n, s] : f.input_groups)
        {
            os << "    const double* " << n << ",\n";
        }
        for (auto [n, s] : f.param_groups)
        {
            if (n != direction_name)
            {
                os << "    const double* " << n << ",\n";
            }
        }
        os << "    " << basename << "_cache_t* cache\n) {\n";
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            if (n.op == Op::Input || (n.op == Op::Param && n.name != direction_name))
            {
                os << "    cache->" << cache_slot(n) << " = " << node_ref(n) << ";\n";
            }
        }
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            if (svm.depends_on_direction[i])
            {
                continue;
            }
            if (n.op == Op::Input || n.op == Op::Param || n.op == Op::Constant)
            {
                continue;
            }
            os << "    cache->t" << i << " = " << expr_rhs(f.graph, i, ref_prepare) << ";\n";
        }
        os << "}\n\n";

        os << "void " << basename << "_apply(\n";
        os << "    const " << basename << "_cache_t* cache,\n";
        os << "    const double* " << direction_name << ",\n";
        os << "    double* out\n) {\n";
        for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i)
        {
            const auto &n = f.graph.nodes[i];
            if (!svm.depends_on_direction[i])
            {
                continue;
            }
            if (n.op == Op::Input || n.op == Op::Param || n.op == Op::Constant)
            {
                continue;
            }
            os << "    double t" << i << " = " << expr_rhs(f.graph, i, ref_apply) << ";\n";
        }
        for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i)
        {
            os << "    out[" << i << "] = " << ref_apply(f.outputs[i]) << ";\n";
        }
        os << "}\n";
    }
};

inline std::string to_c(const GraphFunction &f, const std::string &name)
{
    std::ostringstream os;
    CEmitter::emit_function(f, name, os);
    return os.str();
}

inline std::string to_staged_c(const GraphFunction &f, const std::string &basename, const std::string &direction_name = "v")
{
    std::ostringstream os;
    CEmitter::emit_staged(f, basename, direction_name, os);
    return os.str();
}
// -----------------------------------------------------------------------------
// Structural sparsity via graph reachability on an already-differentiated graph.
//
// Philosophy:
//   Build the highest-order graph you need (Grad, HVP, …), then ask it
//   "which outputs structurally depend on which components of group G?"
//   The optimizer has already folded genuine zeros to Op::Constant(0.0),
//   so reachability on the optimized graph is exact structural sparsity.
//
// Primary primitive:
//   structural_sparsity(f, name)
//     — name is an Input group  → use for Jacobian  (pass F,   probe "x")
//     — name is a Param  group  → use for Hessian   (pass HVP, probe "v")
//
// Convenience wrappers:
//   jacobian_sparsity(F,   "x")           → SparsityPattern
//   hessian_sparsity (HVP, "v")           → lower-triangular (row,col) pairs
//   hessian_sparsity_full(HVP, "v")       → full symmetric   (row,col) pairs
//
// SparsityPattern post-processing (chainable on the returned pattern):
//   .to_pairs()          — raw (output, var) pairs
//   .contract_outputs()  — (output, var) → (row, col) deduplicated
//                          valid when output_i IS x_i (gradient functions)
//
// Static helpers on SparsityPattern:
//   SparsityPattern::lower_triangular(pairs)  — fold to i >= j
//   SparsityPattern::symmetrize(pairs)        — add (j,i) for every (i,j)
// -----------------------------------------------------------------------------

struct SparsityEntry
{
    int output; // which output of the function
    int var;    // which component of the probed group
};

struct SparsityPattern
{
    std::string variable_name;
    int output_size = 0;
    int variable_size = 0;

    std::vector<SparsityEntry> entries;

    std::size_t nnz() const { return entries.size(); }
    bool empty() const { return entries.empty(); }

    // Raw (output, var) pairs.
    std::vector<std::pair<int, int>> to_pairs() const
    {
        std::vector<std::pair<int, int>> out;
        out.reserve(entries.size());
        for (auto &e : entries)
        {
            out.push_back({e.output, e.var});
        }
        return out;
    }

    // Reinterpret (output, var) as (row, col) and deduplicate.
    // Valid when the function's output index i IS x_i (gradient / HVP outputs).
    std::vector<std::pair<int, int>> contract_outputs() const
    {
        std::vector<std::pair<int, int>> out;
        out.reserve(entries.size());
        for (auto &e : entries)
        {
            out.push_back({e.output, e.var});
        }
        std::sort(out.begin(), out.end());
        out.erase(std::unique(out.begin(), out.end()), out.end());
        return out;
    }

    // Fold (i,j) and (j,i) to the lower triangle: i >= j.
    static std::vector<std::pair<int, int>> lower_triangular(std::vector<std::pair<int, int>> pairs)
    {
        for (auto &p : pairs)
        {
            if (p.first < p.second)
            {
                std::swap(p.first, p.second);
            }
        }
        std::sort(pairs.begin(), pairs.end());
        pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
        return pairs;
    }

    // Add (j,i) for every (i,j), keep unique → full symmetric set.
    static std::vector<std::pair<int, int>> symmetrize(std::vector<std::pair<int, int>> pairs)
    {
        const std::size_t n = pairs.size();
        pairs.reserve(n * 2);
        for (std::size_t k = 0; k < n; ++k)
        {
            pairs.push_back({pairs[k].second, pairs[k].first});
        }
        std::sort(pairs.begin(), pairs.end());
        pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
        return pairs;
    }

    std::string str() const
    {
        std::ostringstream os;
        os << "SparsityPattern("
           << "outputs=" << output_size << ", " << variable_name << "[" << variable_size << "]"
           << ", nnz=" << nnz() << ")\n";
        for (auto &e : entries)
        {
            os << "  out[" << e.output << "]  " << variable_name << "[" << e.var << "]\n";
        }
        return os.str();
    }
};

inline std::ostream &operator<<(std::ostream &os, const SparsityPattern &p)
{
    return os << p.str();
}

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

inline int input_group_size(const GraphFunction &f, const std::string &name)
{
    for (auto [n, size] : f.input_groups)
    {
        if (n == name)
        {
            return size;
        }
    }
    return 0;
}

inline int param_group_size(const GraphFunction &f, const std::string &name)
{
    for (auto [n, size] : f.param_groups)
    {
        if (n == name)
        {
            return size;
        }
    }
    return 0;
}

// -----------------------------------------------------------------------------
// Core primitive: reachability-based structural sparsity.
//
// Optimizes f, then marks which output nodes can reach each component of
// the named group (Input or Param).  A node that the optimizer collapsed
// to Op::Constant is unreachable from every variable — exactly right.
// -----------------------------------------------------------------------------

inline SparsityPattern structural_sparsity(const GraphFunction &fn, const std::string &wrt_name)
{
    auto f = optimize(fn);

    // Determine node type and size from the name alone.
    int n = input_group_size(f, wrt_name);

    bool wrt_param = false;
    if (n <= 0)
    {
        n = param_group_size(f, wrt_name);
        wrt_param = true;
    }
    if (n <= 0)
    {
        throw std::runtime_error("structural_sparsity: unknown group: " + wrt_name);
    }


    const int N = static_cast<int>(f.graph.nodes.size());

    // deps[i*n + v] = true  iff  node i structurally depends on wrt_name[v].
    // Filled in topological (post) order: children before parents.
    std::vector<bool> deps(N * n, false);

    for (NodeId i : topo_used(f.graph, f.outputs))
    {
        const auto &nd = f.graph.nodes[i];

        // Seed: this node is one component of the probed group.
        if (!wrt_param && nd.op == Op::Input && nd.name == wrt_name)
        {
            deps[i * n + nd.index] = true;
            continue;
        }
        if (wrt_param && nd.op == Op::Param && nd.name == wrt_name)
        {
            deps[i * n + nd.index] = true;
            continue;
        }

        // Union children's dependency sets into this node.
        if (nd.a >= 0)
        {
            for (int v = 0; v < n; ++v)
            {
                deps[i * n + v] = deps[i * n + v] || deps[nd.a * n + v];
            }
        }
        if (nd.b >= 0)
        {
            for (int v = 0; v < n; ++v)
            {
                deps[i * n + v] = deps[i * n + v] || deps[nd.b * n + v];
            }
        }
    }

    SparsityPattern pat;
    pat.variable_name = wrt_name;
    pat.output_size = static_cast<int>(f.outputs.size());
    pat.variable_size = n;

    for (int out = 0; out < static_cast<int>(f.outputs.size()); ++out)
    {
        const NodeId oid = f.outputs[out];
        for (int v = 0; v < n; ++v)
        {
            if (deps[oid * n + v])
            {
                pat.entries.push_back({out, v});
            }
        }
    }
    return pat;
}

// -----------------------------------------------------------------------------
// Convenience wrappers — pass the graph you already built.
//
//   F   = function_from(...)
//   Grad = reverse_diff(F, "lambda", "x")
//   HVP  = forward_diff(Grad, "x", "v")
//
//   jacobian_sparsity(F)           probe Input "x"  → (output, col) pattern
//   hessian_sparsity(HVP)          probe Param "v"  → lower-triangle (row,col)
//   hessian_sparsity_full(HVP)     probe Param "v"  → full symmetric  (row,col)
// -----------------------------------------------------------------------------

// Jacobian: pass F, probe the input group.
inline SparsityPattern jacobian_sparsity(const GraphFunction &F, const std::string &wrt = "x")
{
    return structural_sparsity(F, wrt);
}

// Hessian lower triangle: pass HVP, probe the direction param.
// Output i of HVP == grad component i == x_i, so contract_outputs()
// gives (H_row, H_col) = (output, v_index) directly.
inline std::vector<std::pair<int, int>> hessian_sparsity(const GraphFunction &HVP, const std::string &v_name = "v")
{
    return SparsityPattern::lower_triangular(structural_sparsity(HVP, v_name).contract_outputs());
}

// Hessian full symmetric: same but both triangles.
inline std::vector<std::pair<int, int>> hessian_sparsity_full(const GraphFunction &HVP, const std::string &v_name = "v")
{
    return SparsityPattern::symmetrize(structural_sparsity(HVP, v_name).contract_outputs());
}

} // namespace advec

#endif // AD_ADVEC_H
