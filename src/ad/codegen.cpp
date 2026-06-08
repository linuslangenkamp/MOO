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

#include "codegen.h"

namespace ad {

std::string CEmitter::number(double v) {
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

std::string CEmitter::cname(const std::string &s) {
    std::string r;
    for (char c : s) {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '_') {
            r.push_back(c);
        } else {
            r.push_back('_');
        }
    }
    if (r.empty() || (r[0] >= '0' && r[0] <= '9')) {
        r = "_" + r;
    }
    return r;
}

std::string CEmitter::cache_slot(const Node &n) { return cname(n.name) + "_" + std::to_string(n.index); }

std::string CEmitter::node_ref(const Node &n) {
    if (n.op == Op::Input) {
        return n.name + "[" + std::to_string(n.index) + "]";
    }
    if (n.op == Op::Param) {
        return n.name + "[" + std::to_string(n.index) + "]";
    }
    throw std::runtime_error("node_ref only supports input/param");
}

std::string CEmitter::expr_rhs(const Graph &g, NodeId id, const std::function<std::string(NodeId)> &ref) {
    const auto &n = g.nodes[id];
    switch (n.op) {
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

bool CEmitter::emit_vector_function(const GraphFunction &f, const std::string &name, std::ostream &os) {
    if (!f.has_vector_structure() || f.output_vector < 0 || f.output_vector >= static_cast<int>(f.vector_nodes.size())) {
        return false;
    }
    const auto &output_node = f.vector_nodes[static_cast<std::size_t>(f.output_vector)];
    if (output_node.size != f.output_size()) {
        return false;
    }

    os << "void " << name << "(\n";
    for (auto [n, s] : f.input_groups) {
        os << "    const double* " << n << ",\n";
    }
    for (auto [n, s] : f.param_groups) {
        os << "    const double* " << n << ",\n";
    }
    os << "    double* out\n) {\n";

    std::unordered_set<NodeId> scalar_emitted;
    std::unordered_set<int> vector_emitted;

    auto direct_group_name = [&](const FunctionVectorNode &node) -> std::string {
        if (node.op != VectorOp::Values || (!node.is_input_group && !node.is_param_group) || static_cast<int>(node.values.size()) != node.size) {
            return "";
        }
        std::string group_name;
        for (int i = 0; i < node.size; ++i) {
            const auto &scalar = f.graph.nodes[static_cast<std::size_t>(node.values[static_cast<std::size_t>(i)])];
            if (node.is_input_group && scalar.op != Op::Input) {
                return "";
            }
            if (node.is_param_group && scalar.op != Op::Param) {
                return "";
            }
            if (scalar.index != i) {
                return "";
            }
            if (i == 0) {
                group_name = scalar.name;
            } else if (scalar.name != group_name) {
                return "";
            }
        }
        return group_name;
    };

    auto vector_ref = [&](int vector_id, const std::string &index) -> std::string {
        const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
        auto group = direct_group_name(node);
        if (!group.empty()) {
            return group + "[" + index + "]";
        }
        return "vec" + std::to_string(vector_id) + "[" + index + "]";
    };

    auto scalar_ref = [&](NodeId id) -> std::string {
        const auto &n = f.graph.nodes[static_cast<std::size_t>(id)];
        if (n.op == Op::Input || n.op == Op::Param) {
            return node_ref(n);
        }
        if (n.op == Op::Constant) {
            return number(n.value);
        }
        return "t" + std::to_string(id);
    };

    std::function<void(NodeId)> emit_scalar = [&](NodeId id) {
        if (id < 0 || scalar_emitted.count(id)) {
            return;
        }
        const auto &n = f.graph.nodes[static_cast<std::size_t>(id)];
        if (n.a >= 0) {
            emit_scalar(n.a);
        }
        if (n.b >= 0) {
            emit_scalar(n.b);
        }
        if (n.op != Op::Input && n.op != Op::Param && n.op != Op::Constant) {
            os << "    double t" << id << " = " << expr_rhs(f.graph, id, scalar_ref) << ";\n";
        }
        scalar_emitted.insert(id);
    };

    std::function<void(int)> emit_vector = [&](int vector_id) {
        if (vector_id < 0 || vector_id >= static_cast<int>(f.vector_nodes.size())) {
            throw std::runtime_error("invalid vector node id");
        }
        if (vector_emitted.count(vector_id)) {
            return;
        }
        const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
        if (node.size < 0) {
            throw std::runtime_error("invalid vector node size");
        }
        if (node.op == VectorOp::Add || node.op == VectorOp::Sub || node.op == VectorOp::Concat) {
            emit_vector(node.a);
            emit_vector(node.b);
        } else if (node.op == VectorOp::Scale || node.op == VectorOp::DenseMatVec || node.op == VectorOp::SparseMatVec || node.op == VectorOp::KronEyeMatVec ||
                   node.op == VectorOp::Slice) {
            emit_vector(node.a);
        } else if (node.op == VectorOp::Values) {
            for (NodeId value : node.values) {
                emit_scalar(value);
            }
        }

        if (!direct_group_name(node).empty()) {
            vector_emitted.insert(vector_id);
            return;
        }

        os << "    double vec" << vector_id << "[" << std::max(node.size, 1) << "] = {0};\n";
        switch (node.op) {
            case VectorOp::Values:
                for (int i = 0; i < node.size; ++i) {
                    os << "    vec" << vector_id << "[" << i << "] = " << scalar_ref(node.values[static_cast<std::size_t>(i)]) << ";\n";
                }
                break;
            case VectorOp::Add:
                os << "    for (int i = 0; i < " << node.size << "; ++i) { vec" << vector_id << "[i] = " << vector_ref(node.a, "i") << " + " << vector_ref(node.b, "i")
                   << "; }\n";
                break;
            case VectorOp::Sub:
                os << "    for (int i = 0; i < " << node.size << "; ++i) { vec" << vector_id << "[i] = " << vector_ref(node.a, "i") << " - " << vector_ref(node.b, "i")
                   << "; }\n";
                break;
            case VectorOp::Scale:
                os << "    for (int i = 0; i < " << node.size << "; ++i) { vec" << vector_id << "[i] = " << number(node.scale) << " * " << vector_ref(node.a, "i")
                   << "; }\n";
                break;
            case VectorOp::DenseMatVec:
                os << "    static const double mat" << vector_id << "[" << std::max(node.matrix.rows * node.matrix.cols, 1) << "] = {";
                for (int i = 0; i < node.matrix.rows * node.matrix.cols; ++i) {
                    if (i > 0) {
                        os << ", ";
                    }
                    os << number(node.matrix.values[static_cast<std::size_t>(i)]);
                }
                if (node.matrix.rows * node.matrix.cols == 0) {
                    os << "0";
                }
                os << "};\n";
                os << "    for (int row = 0; row < " << node.matrix.rows << "; ++row) {\n";
                os << "        double acc = 0.0;\n";
                os << "        for (int col = 0; col < " << node.matrix.cols << "; ++col) { acc += mat" << vector_id << "[row * " << node.matrix.cols << " + col] * "
                   << vector_ref(node.a, "col") << "; }\n";
                os << "        vec" << vector_id << "[row] = acc;\n";
                os << "    }\n";
                break;
            case VectorOp::SparseMatVec:
                os << "    static const int sp_row" << vector_id << "[" << std::max(node.sparse_matrix.nnz(), 1) << "] = {";
                for (int i = 0; i < node.sparse_matrix.nnz(); ++i) {
                    if (i > 0) {
                        os << ", ";
                    }
                    os << node.sparse_matrix.row_indices[static_cast<std::size_t>(i)];
                }
                if (node.sparse_matrix.nnz() == 0) {
                    os << "0";
                }
                os << "};\n";
                os << "    static const int sp_col" << vector_id << "[" << std::max(node.sparse_matrix.nnz(), 1) << "] = {";
                for (int i = 0; i < node.sparse_matrix.nnz(); ++i) {
                    if (i > 0) {
                        os << ", ";
                    }
                    os << node.sparse_matrix.col_indices[static_cast<std::size_t>(i)];
                }
                if (node.sparse_matrix.nnz() == 0) {
                    os << "0";
                }
                os << "};\n";
                os << "    static const double sp_val" << vector_id << "[" << std::max(node.sparse_matrix.nnz(), 1) << "] = {";
                for (int i = 0; i < node.sparse_matrix.nnz(); ++i) {
                    if (i > 0) {
                        os << ", ";
                    }
                    os << number(node.sparse_matrix.values[static_cast<std::size_t>(i)]);
                }
                if (node.sparse_matrix.nnz() == 0) {
                    os << "0";
                }
                os << "};\n";
                os << "    for (int k = 0; k < " << node.sparse_matrix.nnz() << "; ++k) { vec" << vector_id << "[sp_row" << vector_id << "[k]] += sp_val" << vector_id
                   << "[k] * " << vector_ref(node.a, "sp_col" + std::to_string(vector_id) + "[k]") << "; }\n";
                break;
            case VectorOp::KronEyeMatVec:
                os << "    static const double kron_mat" << vector_id << "[" << std::max(node.matrix.rows * node.matrix.cols, 1) << "] = {";
                for (int i = 0; i < node.matrix.rows * node.matrix.cols; ++i) {
                    if (i > 0) {
                        os << ", ";
                    }
                    os << number(node.matrix.values[static_cast<std::size_t>(i)]);
                }
                if (node.matrix.rows * node.matrix.cols == 0) {
                    os << "0";
                }
                os << "};\n";
                os << "    for (int row = 0; row < " << node.matrix.rows << "; ++row) {\n";
                os << "        for (int inner = 0; inner < " << node.eye_size << "; ++inner) {\n";
                os << "            double acc = 0.0;\n";
                os << "            for (int col = 0; col < " << node.matrix.cols << "; ++col) { acc += kron_mat" << vector_id << "[row * " << node.matrix.cols << " + col] * "
                   << vector_ref(node.a, "col * " + std::to_string(node.eye_size) + " + inner") << "; }\n";
                os << "            vec" << vector_id << "[row * " << node.eye_size << " + inner] = acc;\n";
                os << "        }\n";
                os << "    }\n";
                break;
            case VectorOp::Slice:
                os << "    for (int i = 0; i < " << node.size << "; ++i) { vec" << vector_id << "[i] = " << vector_ref(node.a, std::to_string(node.start) + " + i * " + std::to_string(node.stride))
                   << "; }\n";
                break;
            case VectorOp::Concat: {
                const auto lhs_size = f.vector_nodes[static_cast<std::size_t>(node.a)].size;
                const auto rhs_size = f.vector_nodes[static_cast<std::size_t>(node.b)].size;
                os << "    for (int i = 0; i < " << lhs_size << "; ++i) { vec" << vector_id << "[i] = " << vector_ref(node.a, "i") << "; }\n";
                os << "    for (int i = 0; i < " << rhs_size << "; ++i) { vec" << vector_id << "[" << lhs_size << " + i] = " << vector_ref(node.b, "i") << "; }\n";
                break;
            }
        }
        vector_emitted.insert(vector_id);
    };

    emit_vector(f.output_vector);
    os << "    for (int i = 0; i < " << f.output_size() << "; ++i) { out[i] = " << vector_ref(f.output_vector, "i") << "; }\n";
    os << "}\n";
    return true;
}

void CEmitter::emit_function(const GraphFunction &fn, const std::string &name, std::ostream &os) {
    auto f = optimize(fn);
    if (emit_vector_function(f, name, os)) {
        return;
    }
    os << "void " << name << "(\n";
    for (auto [n, s] : f.input_groups) {
        os << "    const double* " << n << ",\n";
    }
    for (auto [n, s] : f.param_groups) {
        os << "    const double* " << n << ",\n";
    }
    os << "    double* out\n) {\n";

    std::unordered_set<NodeId> output_set(f.outputs.begin(), f.outputs.end());
    auto ref = [&](NodeId id) -> std::string {
        const auto &n = f.graph.nodes[id];
        if (n.op == Op::Input || n.op == Op::Param) {
            return node_ref(n);
        }
        if (n.op == Op::Constant) {
            return number(n.value);
        }
        return "t" + std::to_string(id);
    };

    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        if (n.op == Op::Input || n.op == Op::Param || n.op == Op::Constant) {
            continue;
        }
        os << "    double t" << i << " = " << expr_rhs(f.graph, i, ref) << ";\n";
    }
    for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i) {
        os << "    out[" << i << "] = " << ref(f.outputs[i]) << ";\n";
    }
    os << "}\n";
}

void CEmitter::emit_staged(const GraphFunction &fn, const std::string &basename, const std::string &direction_name, std::ostream &os) {
    StagedVM svm(fn, direction_name);
    const auto &f = svm.f;
    os << "typedef struct {\n";
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        if ((n.op == Op::Input) || (n.op == Op::Param && n.name != direction_name)) {
            os << "    double " << cache_slot(n) << ";\n";
        }
    }
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        if (!svm.depends_on_direction[i] && n.op != Op::Input && n.op != Op::Param && n.op != Op::Constant) {
            os << "    double t" << i << ";\n";
        }
    }
    os << "} " << basename << "_cache_t;\n\n";

    auto ref_prepare = [&](NodeId id) -> std::string {
        const auto &n = f.graph.nodes[id];
        if (n.op == Op::Input || n.op == Op::Param) {
            return node_ref(n);
        }
        if (n.op == Op::Constant) {
            return number(n.value);
        }
        return "cache->t" + std::to_string(id);
    };

    auto ref_apply = [&](NodeId id) -> std::string {
        const auto &n = f.graph.nodes[id];
        if (n.op == Op::Input) {
            return "cache->" + cache_slot(n);
        }
        if (n.op == Op::Param) {
            if (n.name == direction_name) {
                return node_ref(n);
            }
            return "cache->" + cache_slot(n);
        }
        if (n.op == Op::Constant) {
            return number(n.value);
        }
        if (!svm.depends_on_direction[id]) {
            return "cache->t" + std::to_string(id);
        }
        return "t" + std::to_string(id);
    };

    os << "void " << basename << "_prepare(\n";
    for (auto [n, s] : f.input_groups) {
        os << "    const double* " << n << ",\n";
    }
    for (auto [n, s] : f.param_groups) {
        if (n != direction_name) {
            os << "    const double* " << n << ",\n";
        }
    }
    os << "    " << basename << "_cache_t* cache\n) {\n";
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        if (n.op == Op::Input || (n.op == Op::Param && n.name != direction_name)) {
            os << "    cache->" << cache_slot(n) << " = " << node_ref(n) << ";\n";
        }
    }
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        if (svm.depends_on_direction[i]) {
            continue;
        }
        if (n.op == Op::Input || n.op == Op::Param || n.op == Op::Constant) {
            continue;
        }
        os << "    cache->t" << i << " = " << expr_rhs(f.graph, i, ref_prepare) << ";\n";
    }
    os << "}\n\n";

    os << "void " << basename << "_apply(\n";
    os << "    const " << basename << "_cache_t* cache,\n";
    os << "    const double* " << direction_name << ",\n";
    os << "    double* out\n) {\n";
    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        const auto &n = f.graph.nodes[i];
        if (!svm.depends_on_direction[i]) {
            continue;
        }
        if (n.op == Op::Input || n.op == Op::Param || n.op == Op::Constant) {
            continue;
        }
        os << "    double t" << i << " = " << expr_rhs(f.graph, i, ref_apply) << ";\n";
    }
    for (int i = 0; i < static_cast<int>(f.outputs.size()); ++i) {
        os << "    out[" << i << "] = " << ref_apply(f.outputs[i]) << ";\n";
    }
    os << "}\n";
}


std::string to_c(const GraphFunction &f, const std::string &name) {
    std::ostringstream os;
    CEmitter::emit_function(f, name, os);
    return os.str();
}

std::string to_staged_c(const GraphFunction &f, const std::string &basename, const std::string &direction_name) {
    std::ostringstream os;
    CEmitter::emit_staged(f, basename, direction_name, os);
    return os.str();
}

std::string to_sparse_coefficients_c(const GraphFunction &linear_f, const std::string &seed_name, const std::vector<std::pair<int, int>> &entries, const std::string &name) {
    return to_c(sparse_coefficients(linear_f, seed_name, entries), name);
}

std::string to_sparse_jacobian_c(const GraphFunction &F,
                                 const std::string &wrt,
                                 const std::vector<std::pair<int, int>> &entries,
                                 const std::string &name,
                                 const std::string &direction_name) {
    return to_c(sparse_jacobian_function(F, wrt, entries, direction_name), name);
}

std::string to_sparse_hessian_c(const GraphFunction &HVP, const std::string &direction_name, const std::vector<std::pair<int, int>> &entries, const std::string &name) {
    return to_c(sparse_hessian_function(HVP, direction_name, entries), name);
}

} // namespace ad
