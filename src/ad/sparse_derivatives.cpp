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

#include "sparse_derivatives.h"

namespace ad {

std::vector<int> greedy_column_coloring(int column_count, const std::vector<std::pair<int, int>> &pattern) {
    std::vector<std::unordered_set<int>> rows_by_col(column_count);
    std::unordered_map<int, std::vector<int>> cols_by_row;
    for (auto [row, col] : pattern) {
        if (col >= 0 && col < column_count) {
            cols_by_row[row].push_back(col);
        }
    }
    std::vector<std::unordered_set<int>> adjacency(column_count);
    for (auto &[row, cols] : cols_by_row) {
        (void)row;
        std::sort(cols.begin(), cols.end());
        cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        for (std::size_t i = 0; i < cols.size(); ++i) {
            for (std::size_t j = i + 1; j < cols.size(); ++j) {
                adjacency[cols[i]].insert(cols[j]);
                adjacency[cols[j]].insert(cols[i]);
            }
        }
    }

    std::vector<int> order(column_count);
    for (int i = 0; i < column_count; ++i) {
        order[i] = i;
    }
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        if (adjacency[a].size() != adjacency[b].size()) {
            return adjacency[a].size() > adjacency[b].size();
        }
        return a < b;
    });

    std::vector<int> color(column_count, -1);
    for (int col : order) {
        std::unordered_set<int> used;
        for (int other : adjacency[col]) {
            if (color[other] >= 0) {
                used.insert(color[other]);
            }
        }
        int c = 0;
        while (used.count(c)) {
            ++c;
        }
        color[col] = c;
    }
    return color;
}

SparseDerivativePlan derivative_plan(int column_count, const std::vector<std::pair<int, int>> &pattern) {
    SparseDerivativePlan plan;
    plan.pattern = pattern;
    if (pattern.empty()) {
        plan.column_color.assign(static_cast<std::size_t>(std::max(column_count, 0)), -1);
        plan.color_count = 0;
        plan.strategy = "direct";
        return plan;
    }
    plan.column_color = greedy_column_coloring(column_count, pattern);
    plan.color_count = 0;
    for (int c : plan.column_color) {
        plan.color_count = std::max(plan.color_count, c + 1);
    }
    plan.strategy = (plan.color_count > 0 && plan.color_count < static_cast<int>(pattern.size())) ? "direct" : "direct";
    return plan;
}

int graph_function_group_size(const GraphFunction &f, const std::string &name) {
    int n = f.input_group_size(name);
    if (n >= 0) {
        return n;
    }
    n = f.param_group_size(name);
    if (n >= 0) {
        return n;
    }
    return 0;
}

bool has_coeff(const LinearForm &f) {
    return !f.coeff.empty();
}

LinearForm make_constant_form(Expr e) {
    LinearForm out;
    out.constant = e;
    return out;
}

LinearForm linear_add(OptimizingBuilder &b, const LinearForm &a, const LinearForm &c, double sign) {
    LinearForm out;
    out.constant = sign > 0.0 ? b.add(a.constant, c.constant) : b.sub(a.constant, c.constant);
    out.coeff = a.coeff;
    for (auto &[idx, expr] : c.coeff) {
        auto it = out.coeff.find(idx);
        Expr term = sign > 0.0 ? expr : b.neg(expr);
        if (it == out.coeff.end()) {
            out.coeff.emplace(idx, term);
        } else {
            it->second = b.add(it->second, term);
        }
    }
    return out;
}

LinearForm linear_scale(OptimizingBuilder &b, const LinearForm &a, Expr scale) {
    LinearForm out;
    out.constant = b.mul(a.constant, scale);
    for (auto &[idx, expr] : a.coeff) {
        out.coeff.emplace(idx, b.mul(expr, scale));
    }
    return out;
}

GraphFunction sparse_coefficients(const GraphFunction &fn, const std::string &seed_name, const std::vector<std::pair<int, int>> &entries) {
    auto f = optimize(fn);
    const int seed_size = param_group_size(f, seed_name);
    if (seed_size <= 0) {
        throw std::runtime_error("sparse_coefficients: unknown seed parameter group: " + seed_name);
    }

    OptimizingBuilder b;
    std::vector<LinearForm> form(f.graph.nodes.size());
    std::vector<Expr> cloned(f.graph.nodes.size());
    auto zero = b.constant(0.0);
    auto one = b.constant(1.0);

    auto require_constant = [&](const LinearForm &lf, const char *op) -> Expr {
        if (has_coeff(lf)) {
            throw std::runtime_error(std::string("sparse_coefficients: nonlinear seed use in ") + op);
        }
        return lf.constant;
    };

    std::vector<NodeId> roots = f.outputs;
    roots.insert(roots.end(), f.inputs.begin(), f.inputs.end());
    roots.insert(roots.end(), f.params.begin(), f.params.end());
    for (NodeId i : topo_used(f.graph, roots)) {
        const auto &n = f.graph.nodes[i];
        switch (n.op) {
            case Op::Constant:
                cloned[i] = b.constant(n.value);
                form[i] = make_constant_form(cloned[i]);
                break;
            case Op::Input:
                cloned[i] = b.input(n.name, n.index);
                form[i] = make_constant_form(cloned[i]);
                break;
            case Op::Param:
                if (n.name == seed_name) {
                    cloned[i] = zero;
                    form[i] = make_constant_form(zero);
                    form[i].coeff[n.index] = one;
                } else {
                    cloned[i] = b.param(n.name, n.index);
                    form[i] = make_constant_form(cloned[i]);
                }
                break;
            case Op::Add:
                form[i] = linear_add(b, form[n.a], form[n.b], 1.0);
                cloned[i] = form[i].constant;
                break;
            case Op::Sub:
                form[i] = linear_add(b, form[n.a], form[n.b], -1.0);
                cloned[i] = form[i].constant;
                break;
            case Op::Neg:
                form[i] = linear_scale(b, form[n.a], b.constant(-1.0));
                cloned[i] = form[i].constant;
                break;
            case Op::Mul: {
                bool la = has_coeff(form[n.a]);
                bool lb = has_coeff(form[n.b]);
                if (la && lb) {
                    throw std::runtime_error("sparse_coefficients: nonlinear seed use in multiplication");
                }
                if (la) {
                    form[i] = linear_scale(b, form[n.a], form[n.b].constant);
                } else if (lb) {
                    form[i] = linear_scale(b, form[n.b], form[n.a].constant);
                } else {
                    form[i] = make_constant_form(b.mul(form[n.a].constant, form[n.b].constant));
                }
                cloned[i] = form[i].constant;
                break;
            }
            case Op::Div: {
                if (has_coeff(form[n.b])) {
                    throw std::runtime_error("sparse_coefficients: nonlinear seed use in division denominator");
                }
                Expr denom = form[n.b].constant;
                form[i].constant = b.div(form[n.a].constant, denom);
                for (auto &[idx, expr] : form[n.a].coeff) {
                    form[i].coeff.emplace(idx, b.div(expr, denom));
                }
                cloned[i] = form[i].constant;
                break;
            }
            case Op::Sin:
                form[i] = make_constant_form(b.unary(Op::Sin, require_constant(form[n.a], "sin")));
                cloned[i] = form[i].constant;
                break;
            case Op::Cos:
                form[i] = make_constant_form(b.unary(Op::Cos, require_constant(form[n.a], "cos")));
                cloned[i] = form[i].constant;
                break;
            case Op::Tan:
                form[i] = make_constant_form(b.unary(Op::Tan, require_constant(form[n.a], "tan")));
                cloned[i] = form[i].constant;
                break;
            case Op::Exp:
                form[i] = make_constant_form(b.unary(Op::Exp, require_constant(form[n.a], "exp")));
                cloned[i] = form[i].constant;
                break;
            case Op::Log:
                form[i] = make_constant_form(b.unary(Op::Log, require_constant(form[n.a], "log")));
                cloned[i] = form[i].constant;
                break;
            case Op::PowConst:
                form[i] = make_constant_form(b.pow_const(require_constant(form[n.a], "pow_const"), n.value));
                cloned[i] = form[i].constant;
                break;
        }
    }

    GraphFunction out;
    out.input_groups = f.input_groups;
    for (auto id : f.inputs) {
        out.inputs.push_back(cloned[id].id);
    }
    for (auto [name, size] : f.param_groups) {
        if (name != seed_name) {
            add_unique_group(out.param_groups, name, size);
        }
    }
    for (auto id : f.params) {
        const auto &n = f.graph.nodes[id];
        if (n.op == Op::Param && n.name != seed_name) {
            out.params.push_back(cloned[id].id);
        }
    }
    for (auto [row, col] : entries) {
        if (row < 0 || row >= static_cast<int>(f.outputs.size()) || col < 0 || col >= seed_size) {
            throw std::runtime_error("sparse_coefficients: requested entry out of range");
        }
        const auto &coeff = form[f.outputs[row]].coeff;
        auto it = coeff.find(col);
        out.outputs.push_back(it == coeff.end() ? zero.id : it->second.id);
    }
    out.graph = std::move(b.g);
    return optimize(out);
}

GraphFunction sparse_jacobian_function(const GraphFunction &F,
                                              const std::string &wrt,
                                              const std::vector<std::pair<int, int>> &entries,
                                              const std::string &direction_name) {
    return sparse_coefficients(forward_diff(F, wrt, direction_name), direction_name, entries);
}

GraphFunction sparse_hessian_function(const GraphFunction &HVP, const std::string &direction_name, const std::vector<std::pair<int, int>> &entries) {
    return sparse_coefficients(HVP, direction_name, entries);
}


} // namespace ad
