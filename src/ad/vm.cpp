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

#include "vm.h"

namespace ad {

EvalEnv &EvalEnv::input(const std::string &name, const double *data) {
    inputs[name] = data;
    return *this;
}

EvalEnv &EvalEnv::param(const std::string &name, const double *data) {
    params[name] = data;
    return *this;
}

VM::VM(GraphFunction fn)
    : f(optimize(fn)) {}

void VM::evaluate(const EvalEnv &env, double *out) const {
    std::vector<double> mem(f.graph.nodes.size(), 0.0);
    std::vector<char> seen(f.graph.nodes.size(), false);
    if (f.has_vector_structure()) {
        auto values = eval_vector_node(f.output_vector, env, mem, seen);
        if (static_cast<int>(values.size()) == f.output_size()) {
            for (int i = 0; i < static_cast<int>(values.size()); ++i) {
                out[i] = values[static_cast<std::size_t>(i)];
            }
            return;
        }
    }
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        eval_scalar_node(i, env, mem, seen);
    }
    for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i) {
        out[i] = mem[f.outputs[i]];
    }
}

double VM::eval_scalar_node(NodeId i, const EvalEnv &env, std::vector<double> &mem, std::vector<char> &seen) const {
    if (i < 0) {
        return 0.0;
    }
    if (seen[static_cast<std::size_t>(i)]) {
        return mem[static_cast<std::size_t>(i)];
    }
    const auto &n = f.graph.nodes[static_cast<std::size_t>(i)];
    switch (n.op) {
        case Op::Constant:
            mem[static_cast<std::size_t>(i)] = n.value;
            break;
        case Op::Input: {
            auto it = env.inputs.find(n.name);
            if (it == env.inputs.end()) {
                throw std::runtime_error("missing input group: " + n.name);
            }
            mem[static_cast<std::size_t>(i)] = it->second[n.index];
            break;
        }
        case Op::Param: {
            auto it = env.params.find(n.name);
            if (it == env.params.end()) {
                throw std::runtime_error("missing param group: " + n.name);
            }
            mem[static_cast<std::size_t>(i)] = it->second[n.index];
            break;
        }
        case Op::Add:
            mem[static_cast<std::size_t>(i)] = eval_scalar_node(n.a, env, mem, seen) + eval_scalar_node(n.b, env, mem, seen);
            break;
        case Op::Sub:
            mem[static_cast<std::size_t>(i)] = eval_scalar_node(n.a, env, mem, seen) - eval_scalar_node(n.b, env, mem, seen);
            break;
        case Op::Mul:
            mem[static_cast<std::size_t>(i)] = eval_scalar_node(n.a, env, mem, seen) * eval_scalar_node(n.b, env, mem, seen);
            break;
        case Op::Div:
            mem[static_cast<std::size_t>(i)] = eval_scalar_node(n.a, env, mem, seen) / eval_scalar_node(n.b, env, mem, seen);
            break;
        case Op::Neg:
            mem[static_cast<std::size_t>(i)] = -eval_scalar_node(n.a, env, mem, seen);
            break;
        case Op::Sin:
            mem[static_cast<std::size_t>(i)] = std::sin(eval_scalar_node(n.a, env, mem, seen));
            break;
        case Op::Cos:
            mem[static_cast<std::size_t>(i)] = std::cos(eval_scalar_node(n.a, env, mem, seen));
            break;
        case Op::Tan:
            mem[static_cast<std::size_t>(i)] = std::tan(eval_scalar_node(n.a, env, mem, seen));
            break;
        case Op::Exp:
            mem[static_cast<std::size_t>(i)] = std::exp(eval_scalar_node(n.a, env, mem, seen));
            break;
        case Op::Log:
            mem[static_cast<std::size_t>(i)] = std::log(eval_scalar_node(n.a, env, mem, seen));
            break;
        case Op::PowConst:
            mem[static_cast<std::size_t>(i)] = std::pow(eval_scalar_node(n.a, env, mem, seen), n.value);
            break;
    }
    seen[static_cast<std::size_t>(i)] = true;
    return mem[static_cast<std::size_t>(i)];
}

std::vector<double> VM::eval_vector_node(int vector_id, const EvalEnv &env, std::vector<double> &mem, std::vector<char> &seen) const {
    if (vector_id < 0 || vector_id >= static_cast<int>(f.vector_nodes.size())) {
        throw std::runtime_error("invalid vector node id");
    }
    const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
    std::vector<double> out(static_cast<std::size_t>(node.size), 0.0);
    switch (node.op) {
        case VectorOp::Values:
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = eval_scalar_node(node.values[static_cast<std::size_t>(i)], env, mem, seen);
            }
            break;
        case VectorOp::Add: {
            auto lhs = eval_vector_node(node.a, env, mem, seen);
            auto rhs = eval_vector_node(node.b, env, mem, seen);
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = lhs[static_cast<std::size_t>(i)] + rhs[static_cast<std::size_t>(i)];
            }
            break;
        }
        case VectorOp::Sub: {
            auto lhs = eval_vector_node(node.a, env, mem, seen);
            auto rhs = eval_vector_node(node.b, env, mem, seen);
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = lhs[static_cast<std::size_t>(i)] - rhs[static_cast<std::size_t>(i)];
            }
            break;
        }
        case VectorOp::Mul: {
            auto lhs = eval_vector_node(node.a, env, mem, seen);
            auto rhs = eval_vector_node(node.b, env, mem, seen);
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = lhs[static_cast<std::size_t>(i)] * rhs[static_cast<std::size_t>(i)];
            }
            break;
        }
        case VectorOp::Scale: {
            auto rhs = eval_vector_node(node.a, env, mem, seen);
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = node.scale * rhs[static_cast<std::size_t>(i)];
            }
            break;
        }
        case VectorOp::PowConst: {
            auto rhs = eval_vector_node(node.a, env, mem, seen);
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = std::pow(rhs[static_cast<std::size_t>(i)], node.power);
            }
            break;
        }
        case VectorOp::DenseMatVec: {
            auto rhs = eval_vector_node(node.a, env, mem, seen);
            for (int row = 0; row < node.matrix.rows; ++row) {
                double acc = 0.0;
                for (int col = 0; col < node.matrix.cols; ++col) {
                    acc += node.matrix(row, col) * rhs[static_cast<std::size_t>(col)];
                }
                out[static_cast<std::size_t>(row)] = acc;
            }
            break;
        }
        case VectorOp::SparseMatVec: {
            auto rhs = eval_vector_node(node.a, env, mem, seen);
            for (int k = 0; k < node.sparse_matrix.nnz(); ++k) {
                int row = node.sparse_matrix.row_indices[static_cast<std::size_t>(k)];
                int col = node.sparse_matrix.col_indices[static_cast<std::size_t>(k)];
                out[static_cast<std::size_t>(row)] += node.sparse_matrix.values[static_cast<std::size_t>(k)] * rhs[static_cast<std::size_t>(col)];
            }
            break;
        }
        case VectorOp::KronEyeMatVec: {
            auto rhs = eval_vector_node(node.a, env, mem, seen);
            for (int row = 0; row < node.matrix.rows; ++row) {
                for (int inner = 0; inner < node.eye_size; ++inner) {
                    double acc = 0.0;
                    for (int col = 0; col < node.matrix.cols; ++col) {
                        acc += node.matrix(row, col) * rhs[static_cast<std::size_t>(col * node.eye_size + inner)];
                    }
                    out[static_cast<std::size_t>(row * node.eye_size + inner)] = acc;
                }
            }
            break;
        }
        case VectorOp::Slice: {
            auto base = eval_vector_node(node.a, env, mem, seen);
            for (int i = 0; i < node.size; ++i) {
                out[static_cast<std::size_t>(i)] = base[static_cast<std::size_t>(node.start + i * node.stride)];
            }
            break;
        }
        case VectorOp::Concat: {
            auto lhs = eval_vector_node(node.a, env, mem, seen);
            auto rhs = eval_vector_node(node.b, env, mem, seen);
            for (int i = 0; i < static_cast<int>(lhs.size()); ++i) {
                out[static_cast<std::size_t>(i)] = lhs[static_cast<std::size_t>(i)];
            }
            for (int i = 0; i < static_cast<int>(rhs.size()); ++i) {
                out[static_cast<std::size_t>(lhs.size() + i)] = rhs[static_cast<std::size_t>(i)];
            }
            break;
        }
    }
    return out;
}

StagedVM::StagedVM(GraphFunction fn, std::string direction)
    : f(optimize(fn)),
      direction_name(std::move(direction)) {
    analyze();
}

void StagedVM::analyze() {
    depends_on_direction.assign(f.graph.nodes.size(), false);
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        bool dep = false;
        if (n.op == Op::Param && n.name == direction_name) {
            dep = true;
        }
        if (n.a >= 0) {
            dep = dep || depends_on_direction[n.a];
        }
        if (n.b >= 0) {
            dep = dep || depends_on_direction[n.b];
        }
        depends_on_direction[i] = dep;
    }
}

void StagedVM::Prepared::apply(const double *direction, double *out) const {
    EvalEnv env;
    env.params.emplace(vm->direction_name, direction);
    std::vector<double> local = mem;
    vm->eval_dependent(env, local, out);
}

StagedVM::Prepared StagedVM::prepare(const EvalEnv &env) const {
    Prepared p;
    p.vm = this;
    p.mem.assign(f.graph.nodes.size(), 0.0);
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        if (depends_on_direction[i]) {
            continue;
        }
        eval_node(i, env, p.mem);
    }
    return p;
}

void StagedVM::eval_dependent(const EvalEnv &env, std::vector<double> &mem, double *out) const {
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        if (!depends_on_direction[i]) {
            continue;
        }
        eval_node(i, env, mem);
    }
    for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i) {
        out[i] = mem[f.outputs[i]];
    }
}

void StagedVM::eval_node(NodeId i, const EvalEnv &env, std::vector<double> &mem) const {
    const auto &n = f.graph.nodes[i];
    switch (n.op) {
        case Op::Constant:
            mem[i] = n.value;
            break;
        case Op::Input: {
            auto it = env.inputs.find(n.name);
            if (it == env.inputs.end()) {
                throw std::runtime_error("missing input group: " + n.name);
            }
            mem[i] = it->second[n.index];
            break;
        }
        case Op::Param: {
            auto it = env.params.find(n.name);
            if (it == env.params.end()) {
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

} // namespace ad
