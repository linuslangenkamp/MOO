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

#include "sparsity.h"

namespace ad {

std::size_t SparsityPattern::nnz() const {
    return entries.size();
}

bool SparsityPattern::empty() const {
    return entries.empty();
}

std::vector<std::pair<int, int>> SparsityPattern::to_pairs() const {
    std::vector<std::pair<int, int>> out;
    out.reserve(entries.size());
    for (auto &e : entries) {
        out.push_back({e.output, e.var});
    }
    return out;
}

std::vector<std::pair<int, int>> SparsityPattern::contract_outputs() const {
    std::vector<std::pair<int, int>> out;
    out.reserve(entries.size());
    for (auto &e : entries) {
        out.push_back({e.output, e.var});
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::vector<std::pair<int, int>> SparsityPattern::lower_triangular(std::vector<std::pair<int, int>> pairs) {
    for (auto &p : pairs) {
        if (p.first < p.second) {
            std::swap(p.first, p.second);
        }
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

std::vector<std::pair<int, int>> SparsityPattern::symmetrize(std::vector<std::pair<int, int>> pairs) {
    const std::size_t n = pairs.size();
    pairs.reserve(n * 2);
    for (std::size_t k = 0; k < n; ++k) {
        pairs.push_back({pairs[k].second, pairs[k].first});
    }
    std::sort(pairs.begin(), pairs.end());
    pairs.erase(std::unique(pairs.begin(), pairs.end()), pairs.end());
    return pairs;
}

std::string SparsityPattern::str() const {
    std::ostringstream os;
    os << "SparsityPattern("
       << "outputs=" << output_size << ", " << variable_name << "[" << variable_size << "]"
       << ", nnz=" << nnz() << ")\n";
    for (auto &e : entries) {
        os << "  out[" << e.output << "]  " << variable_name << "[" << e.var << "]\n";
    }
    return os.str();
}

std::ostream &operator<<(std::ostream &os, const SparsityPattern &p) {
    return os << p.str();
}

SparsityPattern structural_sparsity(const GraphFunction &fn, const std::string &wrt_name) {
    auto f = optimize(fn);
    int n = input_group_size(f, wrt_name);

    bool wrt_param = false;
    if (n <= 0) {
        n = param_group_size(f, wrt_name);
        wrt_param = true;
    }
    if (n <= 0) {
        throw std::runtime_error("structural_sparsity: unknown group: " + wrt_name);
    }

    const int N = static_cast<int>(f.graph.nodes.size());
    std::vector<bool> deps(N * n, false);

    for (NodeId i : topo_used(f.graph, f.outputs)) {
        const auto &nd = f.graph.nodes[i];
        if (!wrt_param && nd.op == Op::Input && nd.name == wrt_name) {
            deps[i * n + nd.index] = true;
            continue;
        }
        if (wrt_param && nd.op == Op::Param && nd.name == wrt_name) {
            deps[i * n + nd.index] = true;
            continue;
        }
        if (nd.a >= 0) {
            for (int v = 0; v < n; ++v) {
                deps[i * n + v] = deps[i * n + v] || deps[nd.a * n + v];
            }
        }
        if (nd.b >= 0) {
            for (int v = 0; v < n; ++v) {
                deps[i * n + v] = deps[i * n + v] || deps[nd.b * n + v];
            }
        }
    }

    SparsityPattern pat;
    pat.variable_name = wrt_name;
    pat.output_size = static_cast<int>(f.outputs.size());
    pat.variable_size = n;

    for (int out = 0; out < static_cast<int>(f.outputs.size()); ++out) {
        const NodeId oid = f.outputs[out];
        for (int v = 0; v < n; ++v) {
            if (deps[oid * n + v]) {
                pat.entries.push_back({out, v});
            }
        }
    }
    return pat;
}

SparsityPattern jacobian_sparsity(const GraphFunction &F, const std::string &wrt) {
    return structural_sparsity(F, wrt);
}

std::vector<std::pair<int, int>> hessian_sparsity(const GraphFunction &HVP, const std::string &v_name) {
    return SparsityPattern::lower_triangular(structural_sparsity(HVP, v_name).contract_outputs());
}

std::vector<std::pair<int, int>> hessian_sparsity_full(const GraphFunction &HVP, const std::string &v_name) {
    return SparsityPattern::symmetrize(structural_sparsity(HVP, v_name).contract_outputs());
}

ExactDerivativePlan exact_derivative_plan(const GraphFunction &F,
                                          const std::string &wrt,
                                          const std::string &direction_name,
                                          const std::string &lambda_name) {
    ExactDerivativePlan plan;
    const int column_count = graph_function_group_size(F, wrt);
    if (column_count <= 0) {
        throw std::runtime_error("exact_derivative_plan: unknown differentiation group: " + wrt);
    }

    plan.jvp = forward_diff(F, wrt, direction_name);
    plan.grad = reverse_diff(F, lambda_name, wrt);
    plan.hvp = forward_diff(plan.grad, wrt, direction_name);
    plan.jacobian_sparsity = jacobian_sparsity(F, wrt).to_pairs();
    plan.hessian_sparsity = hessian_sparsity(plan.hvp, direction_name);
    plan.hessian_full_sparsity = hessian_sparsity_full(plan.hvp, direction_name);

    auto jac_plan = derivative_plan(column_count, plan.jacobian_sparsity);
    auto hes_plan = derivative_plan(column_count, plan.hessian_full_sparsity);
    plan.jacobian_colors = std::move(jac_plan.column_color);
    plan.hessian_colors = std::move(hes_plan.column_color);
    plan.jacobian_color_count = jac_plan.color_count;
    plan.hessian_color_count = hes_plan.color_count;
    return plan;
}

} // namespace ad
