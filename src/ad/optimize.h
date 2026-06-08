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

#ifndef MOO_AD_OPTIMIZE_H
#define MOO_AD_OPTIMIZE_H

#include "core.h"

namespace ad {

struct NodeKey {
    Op op{};
    NodeId a = invalid_node;
    NodeId b = invalid_node;
    double value = 0.0;
    std::string name;
    int index = -1;

    bool operator==(const NodeKey &o) const;
};

struct NodeKeyHash {
    std::size_t operator()(const NodeKey &k) const;
};

struct OptimizingBuilder {
    Graph g;
    std::unordered_map<NodeKey, NodeId, NodeKeyHash> cse;

    Expr make_raw(const Node &n);
    Expr constant(double v);
    Expr input(const std::string &name, int index);
    Expr param(const std::string &name, int index);
    bool is_const(Expr e, double *value = nullptr) const;
    Expr add(Expr x, Expr y);
    Expr sub(Expr x, Expr y);
    Expr mul(Expr x, Expr y);
    Expr div(Expr x, Expr y);
    Expr neg(Expr x);
    Expr unary(Op op, Expr x);
    Expr pow_const(Expr x, double p);
};

Expr ob_binary(OptimizingBuilder &b, Op op, Expr x, Expr y);
std::vector<NodeId> topo_used(const Graph &g, const std::vector<NodeId> &outputs);
GraphFunction optimize(const GraphFunction &f);

} // namespace ad

#endif // MOO_AD_OPTIMIZE_H
