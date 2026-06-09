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

#include "codegen.h"

namespace ad {

namespace {

std::vector<int> scalar_use_counts(const Graph &g) {
    std::vector<int> counts(g.nodes.size(), 0);
    for (const auto &node : g.nodes) {
        if (node.a >= 0) {
            counts[static_cast<std::size_t>(node.a)] += 1;
        }
        if (node.b >= 0) {
            counts[static_cast<std::size_t>(node.b)] += 1;
        }
    }
    return counts;
}

bool cheap_inline_op(Op op) {
    return op == Op::Add || op == Op::Sub || op == Op::Mul || op == Op::Neg;
}

int inline_cost(const Graph &g, NodeId id, const std::vector<int> &use_counts, int limit) {
    if (id < 0) {
        return 0;
    }
    const auto &node = g.nodes[static_cast<std::size_t>(id)];
    if (node.op == Op::Input || node.op == Op::Param || node.op == Op::Constant) {
        return 0;
    }
    if (use_counts[static_cast<std::size_t>(id)] > 1 || !cheap_inline_op(node.op)) {
        return limit + 1;
    }
    int cost = 1;
    cost += inline_cost(g, node.a, use_counts, limit - cost);
    if (cost > limit) {
        return cost;
    }
    cost += inline_cost(g, node.b, use_counts, limit - cost);
    return cost;
}

bool should_inline_node(const Graph &g, NodeId id, const std::vector<int> &use_counts) {
    return inline_cost(g, id, use_counts, 24) <= 24;
}

std::string scalar_expr_ref(const Graph &g, NodeId id, const std::vector<int> &use_counts);

std::string inline_expr_rhs(const Graph &g, NodeId id, const std::vector<int> &use_counts) {
    const auto &node = g.nodes[static_cast<std::size_t>(id)];
    auto ref = [&](NodeId child) { return scalar_expr_ref(g, child, use_counts); };
    switch (node.op) {
        case Op::Constant:
            return CEmitter::number(node.value);
        case Op::Input:
        case Op::Param:
            return CEmitter::node_ref(node);
        case Op::Add:
            return "(" + ref(node.a) + " + " + ref(node.b) + ")";
        case Op::Sub:
            return "(" + ref(node.a) + " - " + ref(node.b) + ")";
        case Op::Mul:
            return "(" + ref(node.a) + " * " + ref(node.b) + ")";
        case Op::Div:
            return "(" + ref(node.a) + " / " + ref(node.b) + ")";
        case Op::Neg:
            return "(-" + ref(node.a) + ")";
        case Op::Sin:
            return "sin(" + ref(node.a) + ")";
        case Op::Cos:
            return "cos(" + ref(node.a) + ")";
        case Op::Tan:
            return "tan(" + ref(node.a) + ")";
        case Op::Exp:
            return "exp(" + ref(node.a) + ")";
        case Op::Log:
            return "log(" + ref(node.a) + ")";
        case Op::PowConst:
            return "pow(" + ref(node.a) + ", " + CEmitter::number(node.value) + ")";
    }
    return "0.0";
}

std::string scalar_expr_ref(const Graph &g, NodeId id, const std::vector<int> &use_counts) {
    const auto &node = g.nodes[static_cast<std::size_t>(id)];
    if (node.op == Op::Input || node.op == Op::Param) {
        return CEmitter::node_ref(node);
    }
    if (node.op == Op::Constant) {
        return CEmitter::number(node.value);
    }
    if (should_inline_node(g, id, use_counts)) {
        return inline_expr_rhs(g, id, use_counts);
    }
    return "t" + std::to_string(id);
}

bool should_emit_scalar_temp(const Graph &g, NodeId id, const std::vector<int> &use_counts) {
    const auto &node = g.nodes[static_cast<std::size_t>(id)];
    if (node.op == Op::Input || node.op == Op::Param || node.op == Op::Constant) {
        return false;
    }
    return !should_inline_node(g, id, use_counts);
}

} // namespace

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

namespace {

void emit_vector_helpers(std::ostream &os) {
    os << "#ifndef MOO_AD_VECTOR_HELPERS_DEFINED\n";
    os << "#define MOO_AD_VECTOR_HELPERS_DEFINED\n";
    os << "static void moo_ad_vec_add(int n, const double* a, const double* b, double* out) {\n";
    os << "    for (int i = 0; i < n; ++i) { out[i] = a[i] + b[i]; }\n";
    os << "}\n";
    os << "static void moo_ad_vec_sub(int n, const double* a, const double* b, double* out) {\n";
    os << "    for (int i = 0; i < n; ++i) { out[i] = a[i] - b[i]; }\n";
    os << "}\n";
    os << "static void moo_ad_vec_mul(int n, const double* a, const double* b, double* out) {\n";
    os << "    for (int i = 0; i < n; ++i) { out[i] = a[i] * b[i]; }\n";
    os << "}\n";
    os << "static void moo_ad_vec_scale(int n, double factor, const double* x, double* out) {\n";
    os << "    for (int i = 0; i < n; ++i) { out[i] = factor * x[i]; }\n";
    os << "}\n";
    os << "static void moo_ad_vec_pow_const(int n, const double* x, double power, double* out) {\n";
    os << "    for (int i = 0; i < n; ++i) { out[i] = pow(x[i], power); }\n";
    os << "}\n";
    os << "static void moo_ad_vec_slice(int n, const double* x, int start, int stride, double* out) {\n";
    os << "    for (int i = 0; i < n; ++i) { out[i] = x[start + i * stride]; }\n";
    os << "}\n";
    os << "static void moo_ad_vec_concat(int lhs_n, int rhs_n, const double* lhs, const double* rhs, double* out) {\n";
    os << "    for (int i = 0; i < lhs_n; ++i) { out[i] = lhs[i]; }\n";
    os << "    for (int i = 0; i < rhs_n; ++i) { out[lhs_n + i] = rhs[i]; }\n";
    os << "}\n";
    os << "static void moo_ad_dense_matvec(int rows, int cols, const double* mat, const double* x, double* out) {\n";
    os << "    for (int row = 0; row < rows; ++row) {\n";
    os << "        double acc = 0.0;\n";
    os << "        for (int col = 0; col < cols; ++col) { acc += mat[row * cols + col] * x[col]; }\n";
    os << "        out[row] = acc;\n";
    os << "    }\n";
    os << "}\n";
    os << "static void moo_ad_sparse_matvec(int nnz, const int* row, const int* col, const double* val, const double* x, double* out) {\n";
    os << "    for (int k = 0; k < nnz; ++k) { out[row[k]] += val[k] * x[col[k]]; }\n";
    os << "}\n";
    os << "static void moo_ad_kron_eye_matvec(int rows, int cols, int eye_size, const double* mat, const double* x, double* out) {\n";
    os << "    for (int row = 0; row < rows; ++row) {\n";
    os << "        for (int inner = 0; inner < eye_size; ++inner) {\n";
    os << "            double acc = 0.0;\n";
    os << "            for (int col = 0; col < cols; ++col) { acc += mat[row * cols + col] * x[col * eye_size + inner]; }\n";
    os << "            out[row * eye_size + inner] = acc;\n";
    os << "        }\n";
    os << "    }\n";
    os << "}\n";
    os << "#endif\n\n";
}

std::string double_array_c(const std::string &indent, const std::string &name, const std::vector<double> &values) {
    std::ostringstream os;
    os << indent << "static const double " << name << "[" << std::max<std::size_t>(values.size(), 1) << "] = {";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            os << ", ";
        }
        os << CEmitter::number(values[i]);
    }
    if (values.empty()) {
        os << "0";
    }
    os << "};\n";
    return os.str();
}

std::string int_array_c(const std::string &indent, const std::string &name, const std::vector<int> &values) {
    std::ostringstream os;
    os << indent << "static const int " << name << "[" << std::max<std::size_t>(values.size(), 1) << "] = {";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i > 0) {
            os << ", ";
        }
        os << values[i];
    }
    if (values.empty()) {
        os << "0";
    }
    os << "};\n";
    return os.str();
}

struct CoeffKey {
    int vector_id = -1;
    int element = -1;
    int seed_col = -1;

    bool operator==(const CoeffKey &other) const { return vector_id == other.vector_id && element == other.element && seed_col == other.seed_col; }
};

struct CoeffKeyHash {
    std::size_t operator()(const CoeffKey &key) const {
        std::size_t h = std::hash<int>{}(key.vector_id);
        h ^= std::hash<int>{}(key.element) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        h ^= std::hash<int>{}(key.seed_col) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
        return h;
    }
};

struct SeedLinearAnalysis {
    OptimizingBuilder builder;
    std::vector<LinearForm> forms;
    std::vector<Expr> cloned;
    Expr zero;
    std::unordered_map<CoeffKey, LinearForm, CoeffKeyHash> vector_form_cache;
};

struct MixedCoeff {
    double constant = 0.0;
    Expr dynamic;
};

struct MixedCoeffAnalysis {
    SeedLinearAnalysis &seed;
    std::unordered_map<CoeffKey, MixedCoeff, CoeffKeyHash> cache;
};

bool expr_constant(const Graph &graph, Expr expr, double *value = nullptr) {
    if (!expr || expr.id < 0 || expr.id >= static_cast<NodeId>(graph.nodes.size())) {
        return false;
    }
    const auto &node = graph.nodes[static_cast<std::size_t>(expr.id)];
    if (node.op != Op::Constant) {
        return false;
    }
    if (value) {
        *value = node.value;
    }
    return true;
}

Expr require_seed_constant(SeedLinearAnalysis &analysis, const LinearForm &form, const char *op) {
    if (has_coeff(form)) {
        throw std::runtime_error(std::string("vector sparse coefficients: nonlinear seed use in ") + op);
    }
    return form.constant;
}

SeedLinearAnalysis analyze_seed_linear_forms(const GraphFunction &f, const std::string &seed_name) {
    SeedLinearAnalysis analysis;
    analysis.forms.resize(f.graph.nodes.size());
    analysis.cloned.resize(f.graph.nodes.size());
    analysis.zero = analysis.builder.constant(0.0);
    Expr one = analysis.builder.constant(1.0);

    std::vector<NodeId> roots = f.outputs;
    roots.insert(roots.end(), f.inputs.begin(), f.inputs.end());
    roots.insert(roots.end(), f.params.begin(), f.params.end());
    if (f.has_vector_structure()) {
        std::vector<char> seen(f.vector_nodes.size(), false);
        std::function<void(int)> add_vector_roots = [&](int vector_id) {
            if (vector_id < 0 || vector_id >= static_cast<int>(f.vector_nodes.size()) || seen[static_cast<std::size_t>(vector_id)]) {
                return;
            }
            seen[static_cast<std::size_t>(vector_id)] = true;
            const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
            roots.insert(roots.end(), node.values.begin(), node.values.end());
            if (node.op == VectorOp::Add || node.op == VectorOp::Sub || node.op == VectorOp::Mul || node.op == VectorOp::Concat) {
                add_vector_roots(node.a);
                add_vector_roots(node.b);
            } else if (node.op == VectorOp::Scale || node.op == VectorOp::PowConst || node.op == VectorOp::DenseMatVec || node.op == VectorOp::SparseMatVec ||
                       node.op == VectorOp::KronEyeMatVec || node.op == VectorOp::Slice) {
                add_vector_roots(node.a);
            }
        };
        add_vector_roots(f.output_vector);
        add_vector_roots(f.input_vector);
        for (int vector_id : f.param_vectors) {
            add_vector_roots(vector_id);
        }
    }
    for (NodeId i : topo_used(f.graph, roots)) {
        const auto &n = f.graph.nodes[static_cast<std::size_t>(i)];
        switch (n.op) {
            case Op::Constant:
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.builder.constant(n.value);
                analysis.forms[static_cast<std::size_t>(i)] = make_constant_form(analysis.cloned[static_cast<std::size_t>(i)]);
                break;
            case Op::Input:
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.builder.input(n.name, n.index);
                analysis.forms[static_cast<std::size_t>(i)] = make_constant_form(analysis.cloned[static_cast<std::size_t>(i)]);
                break;
            case Op::Param:
                if (n.name == seed_name) {
                    analysis.cloned[static_cast<std::size_t>(i)] = analysis.zero;
                    analysis.forms[static_cast<std::size_t>(i)] = make_constant_form(analysis.zero);
                    analysis.forms[static_cast<std::size_t>(i)].coeff[n.index] = one;
                } else {
                    analysis.cloned[static_cast<std::size_t>(i)] = analysis.builder.param(n.name, n.index);
                    analysis.forms[static_cast<std::size_t>(i)] = make_constant_form(analysis.cloned[static_cast<std::size_t>(i)]);
                }
                break;
            case Op::Add:
                analysis.forms[static_cast<std::size_t>(i)] =
                    linear_add(analysis.builder, analysis.forms[static_cast<std::size_t>(n.a)], analysis.forms[static_cast<std::size_t>(n.b)], 1.0);
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Sub:
                analysis.forms[static_cast<std::size_t>(i)] =
                    linear_add(analysis.builder, analysis.forms[static_cast<std::size_t>(n.a)], analysis.forms[static_cast<std::size_t>(n.b)], -1.0);
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Neg:
                analysis.forms[static_cast<std::size_t>(i)] = linear_scale(analysis.builder, analysis.forms[static_cast<std::size_t>(n.a)], analysis.builder.constant(-1.0));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Mul: {
                bool lhs_linear = has_coeff(analysis.forms[static_cast<std::size_t>(n.a)]);
                bool rhs_linear = has_coeff(analysis.forms[static_cast<std::size_t>(n.b)]);
                if (lhs_linear && rhs_linear) {
                    throw std::runtime_error("vector sparse coefficients: nonlinear seed use in multiplication");
                }
                if (lhs_linear) {
                    analysis.forms[static_cast<std::size_t>(i)] =
                        linear_scale(analysis.builder, analysis.forms[static_cast<std::size_t>(n.a)], analysis.forms[static_cast<std::size_t>(n.b)].constant);
                } else if (rhs_linear) {
                    analysis.forms[static_cast<std::size_t>(i)] =
                        linear_scale(analysis.builder, analysis.forms[static_cast<std::size_t>(n.b)], analysis.forms[static_cast<std::size_t>(n.a)].constant);
                } else {
                    analysis.forms[static_cast<std::size_t>(i)] =
                        make_constant_form(analysis.builder.mul(analysis.forms[static_cast<std::size_t>(n.a)].constant, analysis.forms[static_cast<std::size_t>(n.b)].constant));
                }
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            }
            case Op::Div: {
                if (has_coeff(analysis.forms[static_cast<std::size_t>(n.b)])) {
                    throw std::runtime_error("vector sparse coefficients: nonlinear seed use in division denominator");
                }
                Expr denom = analysis.forms[static_cast<std::size_t>(n.b)].constant;
                LinearForm form;
                form.constant = analysis.builder.div(analysis.forms[static_cast<std::size_t>(n.a)].constant, denom);
                for (auto &[idx, expr] : analysis.forms[static_cast<std::size_t>(n.a)].coeff) {
                    form.coeff.emplace(idx, analysis.builder.div(expr, denom));
                }
                analysis.forms[static_cast<std::size_t>(i)] = std::move(form);
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            }
            case Op::Sin:
                analysis.forms[static_cast<std::size_t>(i)] =
                    make_constant_form(analysis.builder.unary(Op::Sin, require_seed_constant(analysis, analysis.forms[static_cast<std::size_t>(n.a)], "sin")));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Cos:
                analysis.forms[static_cast<std::size_t>(i)] =
                    make_constant_form(analysis.builder.unary(Op::Cos, require_seed_constant(analysis, analysis.forms[static_cast<std::size_t>(n.a)], "cos")));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Tan:
                analysis.forms[static_cast<std::size_t>(i)] =
                    make_constant_form(analysis.builder.unary(Op::Tan, require_seed_constant(analysis, analysis.forms[static_cast<std::size_t>(n.a)], "tan")));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Exp:
                analysis.forms[static_cast<std::size_t>(i)] =
                    make_constant_form(analysis.builder.unary(Op::Exp, require_seed_constant(analysis, analysis.forms[static_cast<std::size_t>(n.a)], "exp")));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::Log:
                analysis.forms[static_cast<std::size_t>(i)] =
                    make_constant_form(analysis.builder.unary(Op::Log, require_seed_constant(analysis, analysis.forms[static_cast<std::size_t>(n.a)], "log")));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
            case Op::PowConst:
                analysis.forms[static_cast<std::size_t>(i)] = make_constant_form(
                    analysis.builder.pow_const(require_seed_constant(analysis, analysis.forms[static_cast<std::size_t>(n.a)], "pow_const"), n.value));
                analysis.cloned[static_cast<std::size_t>(i)] = analysis.forms[static_cast<std::size_t>(i)].constant;
                break;
        }
    }
    return analysis;
}

LinearForm require_vector_seed_constant(SeedLinearAnalysis &analysis, const LinearForm &form, const char *op) {
    if (has_coeff(form)) {
        throw std::runtime_error(std::string("vector sparse coefficients: nonlinear seed use in vector ") + op);
    }
    return form;
}

LinearForm vector_linear_form(const GraphFunction &f, SeedLinearAnalysis &analysis, int vector_id, int element) {
    if (vector_id < 0 || vector_id >= static_cast<int>(f.vector_nodes.size())) {
        throw std::runtime_error("vector sparse coefficients: invalid vector node id");
    }
    const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
    if (element < 0 || element >= node.size) {
        throw std::runtime_error("vector sparse coefficients: vector element out of bounds");
    }
    CoeffKey key{vector_id, element, -1};
    auto cached = analysis.vector_form_cache.find(key);
    if (cached != analysis.vector_form_cache.end()) {
        return cached->second;
    }

    auto add_scaled = [&](LinearForm acc, double scale, const LinearForm &term) {
        if (exactly(scale, 0.0)) {
            return acc;
        }
        auto scaled = exactly(scale, 1.0) ? term : linear_scale(analysis.builder, term, analysis.builder.constant(scale));
        return linear_add(analysis.builder, acc, scaled);
    };

    LinearForm result = make_constant_form(analysis.zero);
    switch (node.op) {
        case VectorOp::Values:
            result = analysis.forms[static_cast<std::size_t>(node.values[static_cast<std::size_t>(element)])];
            break;
        case VectorOp::Add:
            result = linear_add(analysis.builder, vector_linear_form(f, analysis, node.a, element), vector_linear_form(f, analysis, node.b, element));
            break;
        case VectorOp::Sub:
            result = linear_add(analysis.builder, vector_linear_form(f, analysis, node.a, element), vector_linear_form(f, analysis, node.b, element), -1.0);
            break;
        case VectorOp::Mul: {
            auto lhs = vector_linear_form(f, analysis, node.a, element);
            auto rhs = vector_linear_form(f, analysis, node.b, element);
            const bool lhs_linear = has_coeff(lhs);
            const bool rhs_linear = has_coeff(rhs);
            if (lhs_linear && rhs_linear) {
                throw std::runtime_error("vector sparse coefficients: nonlinear seed use in vector multiplication");
            }
            if (lhs_linear) {
                result = linear_scale(analysis.builder, lhs, rhs.constant);
            } else if (rhs_linear) {
                result = linear_scale(analysis.builder, rhs, lhs.constant);
            } else {
                result = make_constant_form(analysis.builder.mul(lhs.constant, rhs.constant));
            }
            break;
        }
        case VectorOp::Scale:
            result = linear_scale(analysis.builder, vector_linear_form(f, analysis, node.a, element), analysis.builder.constant(node.scale));
            break;
        case VectorOp::PowConst: {
            auto rhs = require_vector_seed_constant(analysis, vector_linear_form(f, analysis, node.a, element), "pow_const");
            result = make_constant_form(analysis.builder.pow_const(rhs.constant, node.power));
            break;
        }
        case VectorOp::Slice:
            result = vector_linear_form(f, analysis, node.a, node.start + element * node.stride);
            break;
        case VectorOp::Concat: {
            const int lhs_size = f.vector_nodes[static_cast<std::size_t>(node.a)].size;
            if (element < lhs_size) {
                result = vector_linear_form(f, analysis, node.a, element);
            } else {
                result = vector_linear_form(f, analysis, node.b, element - lhs_size);
            }
            break;
        }
        case VectorOp::DenseMatVec:
            for (int col = 0; col < node.matrix.cols; ++col) {
                result = add_scaled(result, node.matrix(element, col), vector_linear_form(f, analysis, node.a, col));
            }
            break;
        case VectorOp::SparseMatVec:
            for (int k = 0; k < node.sparse_matrix.nnz(); ++k) {
                if (node.sparse_matrix.row_indices[static_cast<std::size_t>(k)] == element) {
                    result = add_scaled(result, node.sparse_matrix.values[static_cast<std::size_t>(k)],
                                        vector_linear_form(f, analysis, node.a, node.sparse_matrix.col_indices[static_cast<std::size_t>(k)]));
                }
            }
            break;
        case VectorOp::KronEyeMatVec: {
            const int row = element / node.eye_size;
            const int inner = element % node.eye_size;
            for (int col = 0; col < node.matrix.cols; ++col) {
                result = add_scaled(result, node.matrix(row, col), vector_linear_form(f, analysis, node.a, col * node.eye_size + inner));
            }
            break;
        }
    }
    analysis.vector_form_cache.emplace(key, result);
    return result;
}

bool mixed_has_dynamic(const MixedCoeff &coeff) { return static_cast<bool>(coeff.dynamic); }

MixedCoeff mixed_from_expr(SeedLinearAnalysis &analysis, Expr expr) {
    double value = 0.0;
    if (expr_constant(analysis.builder.g, expr, &value)) {
        return MixedCoeff{value, Expr{}};
    }
    return MixedCoeff{0.0, expr};
}

MixedCoeff mixed_vector_coeff(const GraphFunction &f, MixedCoeffAnalysis &analysis, int vector_id, int element, int seed_col) {
    if (vector_id < 0 || vector_id >= static_cast<int>(f.vector_nodes.size())) {
        throw std::runtime_error("vector sparse coefficients: invalid vector node id");
    }
    const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
    if (element < 0 || element >= node.size) {
        throw std::runtime_error("vector sparse coefficients: vector element out of bounds");
    }
    CoeffKey key{vector_id, element, seed_col};
    auto cached = analysis.cache.find(key);
    if (cached != analysis.cache.end()) {
        return cached->second;
    }

    (void)node;
    MixedCoeff result;
    auto form = vector_linear_form(f, analysis.seed, vector_id, element);
    auto it = form.coeff.find(seed_col);
    result = it == form.coeff.end() ? MixedCoeff{} : mixed_from_expr(analysis.seed, it->second);
    analysis.cache.emplace(key, result);
    return result;
}

std::optional<std::string> try_vector_sparse_coefficients_c(const GraphFunction &linear_f,
                                                            const std::string &seed_name,
                                                            const std::vector<std::pair<int, int>> &entries,
                                                            const std::string &name) {
    auto f = optimize(linear_f);
    if (!f.has_vector_structure() || f.output_vector < 0 || f.output_vector >= static_cast<int>(f.vector_nodes.size()) ||
        f.vector_nodes[static_cast<std::size_t>(f.output_vector)].size != f.output_size()) {
        return std::nullopt;
    }

    const int seed_size = param_group_size(f, seed_name);
    if (seed_size <= 0) {
        return std::nullopt;
    }

    SeedLinearAnalysis analysis = analyze_seed_linear_forms(f, seed_name);
    MixedCoeffAnalysis mixed_analysis{analysis, {}};
    std::vector<int> constant_indices;
    std::vector<double> constant_values;
    std::vector<std::tuple<int, double, Expr>> dynamic_values;
    constant_indices.reserve(entries.size());
    constant_values.reserve(entries.size());

    for (std::size_t out_i = 0; out_i < entries.size(); ++out_i) {
        const auto [row, col] = entries[out_i];
        if (row < 0 || row >= f.output_size() || col < 0 || col >= seed_size) {
            throw std::runtime_error("vector sparse coefficients: requested entry out of range");
        }
        MixedCoeff coeff = mixed_vector_coeff(f, mixed_analysis, f.output_vector, row, col);
        if (mixed_has_dynamic(coeff)) {
            dynamic_values.push_back({static_cast<int>(out_i), coeff.constant, coeff.dynamic});
        } else {
            if (!exactly(coeff.constant, 0.0)) {
                constant_indices.push_back(static_cast<int>(out_i));
                constant_values.push_back(coeff.constant);
            }
        }
    }

    std::ostringstream os;
    os << "void " << name << "(\n";
    for (auto [group_name, size] : f.input_groups) {
        os << "    const double* " << group_name << ",\n";
    }
    for (auto [group_name, size] : f.param_groups) {
        if (group_name != seed_name) {
            os << "    const double* " << group_name << ",\n";
        }
    }
    os << "    double* out\n) {\n";
    os << "    for (int i = 0; i < " << entries.size() << "; ++i) { out[i] = 0.0; }\n";
    if (!constant_indices.empty()) {
        os << int_array_c("    ", "constant_index", constant_indices);
        os << double_array_c("    ", "constant_value", constant_values);
        os << "    for (int i = 0; i < " << constant_indices.size() << "; ++i) { out[constant_index[i]] = constant_value[i]; }\n";
    }

    std::vector<NodeId> dynamic_roots;
    dynamic_roots.reserve(dynamic_values.size());
    for (auto &[idx, constant, expr] : dynamic_values) {
        (void)idx;
        (void)constant;
        dynamic_roots.push_back(expr.id);
    }
    const auto use_counts = scalar_use_counts(analysis.builder.g);
    auto ref = [&](NodeId id) -> std::string {
        return scalar_expr_ref(analysis.builder.g, id, use_counts);
    };
    for (NodeId id : topo_used(analysis.builder.g, dynamic_roots)) {
        if (!should_emit_scalar_temp(analysis.builder.g, id, use_counts)) {
            continue;
        }
        os << "    double t" << id << " = " << CEmitter::expr_rhs(analysis.builder.g, id, ref) << ";\n";
    }
    for (auto &[idx, constant, expr] : dynamic_values) {
        if (exactly(constant, 0.0)) {
            os << "    out[" << idx << "] = " << ref(expr.id) << ";\n";
        } else {
            os << "    out[" << idx << "] = " << CEmitter::number(constant) << " + " << ref(expr.id) << ";\n";
        }
    }
    os << "}\n";
    return os.str();
}

bool vector_sparse_error_can_fallback(const std::runtime_error &err) {
    const std::string message = err.what();
    const std::string prefix = "vector sparse coefficients:";
    if (message.rfind(prefix, 0) != 0) {
        return false;
    }
    return message.find("requested entry out of range") == std::string::npos;
}

} // namespace

bool CEmitter::emit_vector_function(const GraphFunction &f, const std::string &name, std::ostream &os) {
    if (!f.has_vector_structure() || f.output_vector < 0 || f.output_vector >= static_cast<int>(f.vector_nodes.size())) {
        return false;
    }
    const auto &output_node = f.vector_nodes[static_cast<std::size_t>(f.output_vector)];
    if (output_node.size != f.output_size()) {
        return false;
    }

    std::function<bool(int)> uses_vector_helper = [&](int vector_id) -> bool {
        const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
        switch (node.op) {
            case VectorOp::Values:
                return false;
            case VectorOp::Add:
            case VectorOp::Sub:
            case VectorOp::Mul:
            case VectorOp::Concat:
            case VectorOp::Scale:
            case VectorOp::PowConst:
            case VectorOp::DenseMatVec:
            case VectorOp::SparseMatVec:
            case VectorOp::KronEyeMatVec:
            case VectorOp::Slice:
                return true;
        }
        return false;
    };
    if (uses_vector_helper(f.output_vector)) {
        emit_vector_helpers(os);
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
    const auto use_counts = scalar_use_counts(f.graph);

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

    auto vector_ptr = [&](int vector_id) -> std::string {
        const auto &node = f.vector_nodes[static_cast<std::size_t>(vector_id)];
        auto group = direct_group_name(node);
        if (!group.empty()) {
            return group;
        }
        return "vec" + std::to_string(vector_id);
    };

    auto scalar_ref = [&](NodeId id) -> std::string {
        return scalar_expr_ref(f.graph, id, use_counts);
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
        if (should_emit_scalar_temp(f.graph, id, use_counts)) {
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
        if (node.op == VectorOp::Add || node.op == VectorOp::Sub || node.op == VectorOp::Mul || node.op == VectorOp::Concat) {
            emit_vector(node.a);
            emit_vector(node.b);
        } else if (node.op == VectorOp::Scale || node.op == VectorOp::PowConst || node.op == VectorOp::DenseMatVec || node.op == VectorOp::SparseMatVec ||
                   node.op == VectorOp::KronEyeMatVec || node.op == VectorOp::Slice) {
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
                os << "    moo_ad_vec_add(" << node.size << ", " << vector_ptr(node.a) << ", " << vector_ptr(node.b) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::Sub:
                os << "    moo_ad_vec_sub(" << node.size << ", " << vector_ptr(node.a) << ", " << vector_ptr(node.b) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::Mul:
                os << "    moo_ad_vec_mul(" << node.size << ", " << vector_ptr(node.a) << ", " << vector_ptr(node.b) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::Scale:
                os << "    moo_ad_vec_scale(" << node.size << ", " << number(node.scale) << ", " << vector_ptr(node.a) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::PowConst:
                os << "    moo_ad_vec_pow_const(" << node.size << ", " << vector_ptr(node.a) << ", " << number(node.power) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::DenseMatVec:
                os << double_array_c("    ", "mat" + std::to_string(vector_id), node.matrix.values);
                os << "    moo_ad_dense_matvec(" << node.matrix.rows << ", " << node.matrix.cols << ", mat" << vector_id << ", " << vector_ptr(node.a) << ", vec" << vector_id
                   << ");\n";
                break;
            case VectorOp::SparseMatVec:
                os << int_array_c("    ", "sp_row" + std::to_string(vector_id), node.sparse_matrix.row_indices);
                os << int_array_c("    ", "sp_col" + std::to_string(vector_id), node.sparse_matrix.col_indices);
                os << double_array_c("    ", "sp_val" + std::to_string(vector_id), node.sparse_matrix.values);
                os << "    moo_ad_sparse_matvec(" << node.sparse_matrix.nnz() << ", sp_row" << vector_id << ", sp_col" << vector_id << ", sp_val" << vector_id << ", "
                   << vector_ptr(node.a) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::KronEyeMatVec:
                os << double_array_c("    ", "kron_mat" + std::to_string(vector_id), node.matrix.values);
                os << "    moo_ad_kron_eye_matvec(" << node.matrix.rows << ", " << node.matrix.cols << ", " << node.eye_size << ", kron_mat" << vector_id << ", "
                   << vector_ptr(node.a) << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::Slice:
                os << "    moo_ad_vec_slice(" << node.size << ", " << vector_ptr(node.a) << ", " << node.start << ", " << node.stride << ", vec" << vector_id << ");\n";
                break;
            case VectorOp::Concat: {
                const auto lhs_size = f.vector_nodes[static_cast<std::size_t>(node.a)].size;
                const auto rhs_size = f.vector_nodes[static_cast<std::size_t>(node.b)].size;
                os << "    moo_ad_vec_concat(" << lhs_size << ", " << rhs_size << ", " << vector_ptr(node.a) << ", " << vector_ptr(node.b) << ", vec" << vector_id << ");\n";
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

    const auto use_counts = scalar_use_counts(f.graph);
    auto ref = [&](NodeId id) -> std::string {
        return scalar_expr_ref(f.graph, id, use_counts);
    };

    for (NodeId i = 0; i < static_cast<NodeId>(f.graph.nodes.size()); ++i) {
        if (!should_emit_scalar_temp(f.graph, i, use_counts)) {
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
    try {
        if (auto vector_code = try_vector_sparse_coefficients_c(linear_f, seed_name, entries, name)) {
            return *vector_code;
        }
    } catch (const std::runtime_error &err) {
        if (!vector_sparse_error_can_fallback(err)) {
            throw;
        }
    }
    return to_c(sparse_coefficients(linear_f, seed_name, entries), name);
}

std::string to_sparse_jacobian_c(const GraphFunction &F,
                                 const std::string &wrt,
                                 const std::vector<std::pair<int, int>> &entries,
                                 const std::string &name,
                                 const std::string &direction_name) {
    return to_sparse_coefficients_c(forward_diff(F, wrt, direction_name), direction_name, entries, name);
}

std::string to_sparse_hessian_c(const GraphFunction &HVP, const std::string &direction_name, const std::vector<std::pair<int, int>> &entries, const std::string &name) {
    return to_sparse_coefficients_c(HVP, direction_name, entries, name);
}

std::string derivative_callback_mode(const std::string &strategy, const std::vector<std::pair<int, int>> &pairs, const std::vector<int> &colors) {
    (void)pairs;
    (void)colors;
    if (strategy == "colored") {
        return "colored";
    }
    return "direct";
}

std::vector<std::string> derivative_section_keys(const std::string &prefix, const std::string &jac_mode, const std::string &hes_mode) {
    const std::string stem = prefix.empty() ? "" : prefix + "_";
    return {
        stem + "VALUE",
        stem + (jac_mode == "direct" ? "JAC" : "JVP"),
        stem + (hes_mode == "direct" ? "HES" : "HVP"),
    };
}

namespace {

std::string render_colored_fill(const std::vector<std::pair<int, int>> &pairs,
                                const std::vector<int> &colors,
                                const std::string &call,
                                const std::vector<int> &buf_indices,
                                const std::string &indent) {
    if (pairs.empty()) {
        return indent + "(void)out;";
    }

    std::vector<int> indices = buf_indices;
    if (indices.empty()) {
        indices.resize(pairs.size());
        for (int i = 0; i < static_cast<int>(indices.size()); ++i) {
            indices[static_cast<std::size_t>(i)] = i;
        }
    }
    if (indices.size() != pairs.size()) {
        throw std::runtime_error("render_colored_fill: buf_indices must be empty or match sparsity size");
    }

    std::map<int, std::vector<std::tuple<int, int, int>>> by_color;
    for (std::size_t i = 0; i < pairs.size(); ++i) {
        const auto [row, col] = pairs[i];
        const int color = (col >= 0 && col < static_cast<int>(colors.size())) ? colors[static_cast<std::size_t>(col)] : col;
        by_color[color].push_back({indices[i], row, col});
    }

    std::vector<int> color_cols;
    std::vector<int> color_offsets{0};
    std::vector<int> scatter_buf;
    std::vector<int> scatter_row;
    std::vector<int> scatter_offsets{0};
    for (const auto &[color, entries] : by_color) {
        (void)color;
        std::vector<int> cols;
        cols.reserve(entries.size());
        for (const auto &[buf, row, col] : entries) {
            (void)buf;
            (void)row;
            cols.push_back(col);
        }
        std::sort(cols.begin(), cols.end());
        cols.erase(std::unique(cols.begin(), cols.end()), cols.end());
        color_cols.insert(color_cols.end(), cols.begin(), cols.end());
        color_offsets.push_back(static_cast<int>(color_cols.size()));
        for (const auto &[buf, row, col] : entries) {
            (void)col;
            scatter_buf.push_back(buf);
            scatter_row.push_back(row);
        }
        scatter_offsets.push_back(static_cast<int>(scatter_buf.size()));
    }

    std::ostringstream os;
    os << int_array_c(indent, "color_offsets", color_offsets);
    os << int_array_c(indent, "color_cols", color_cols);
    os << int_array_c(indent, "scatter_offsets", scatter_offsets);
    os << int_array_c(indent, "scatter_buf", scatter_buf);
    os << int_array_c(indent, "scatter_row", scatter_row);
    os << indent << "for (int color = 0; color < " << std::max<int>(static_cast<int>(color_offsets.size()) - 1, 0) << "; ++color) {\n";
    os << indent << "    for (int i = color_offsets[color]; i < color_offsets[color + 1]; ++i) { v[color_cols[i]] = 1.0; }\n";
    os << indent << "    " << call << "\n";
    os << indent << "    for (int i = scatter_offsets[color]; i < scatter_offsets[color + 1]; ++i) { out[scatter_buf[i]] = tmp[scatter_row[i]]; }\n";
    os << indent << "    for (int i = color_offsets[color]; i < color_offsets[color + 1]; ++i) { v[color_cols[i]] = 0.0; }\n";
    os << indent << "}";
    return os.str();
}

} // namespace

std::string render_jacobian_callback_body(const std::string &mode,
                                          const std::string &direct_call,
                                          const std::string &input_size,
                                          const std::string &output_size,
                                          const std::vector<std::pair<int, int>> &pairs,
                                          const std::vector<int> &colors,
                                          const std::string &colored_call) {
    if (mode == "direct") {
        return "    " + direct_call;
    }
    std::ostringstream os;
    os << "    f64 v[" << input_size << "] = {0};\n";
    os << "    f64 tmp[" << output_size << "] = {0};\n";
    os << render_colored_fill(pairs, colors, colored_call, {}, "    ");
    return os.str();
}

std::string render_hessian_callback_body(const std::string &mode,
                                         const std::string &direct_body,
                                         const std::string &input_size,
                                         const std::string &tmp_size,
                                         const std::vector<std::pair<int, int>> &pairs,
                                         const std::vector<int> &colors,
                                         const std::string &prepare_body,
                                         const std::string &apply_call,
                                         const std::vector<int> &buf_indices) {
    if (mode == "direct") {
        return direct_body;
    }
    std::ostringstream os;
    os << "    f64 v[" << input_size << "] = {0};\n";
    os << "    f64 tmp[" << tmp_size << "] = {0};\n";
    os << prepare_body;
    if (!prepare_body.empty() && prepare_body.back() != '\n') {
        os << "\n";
    }
    os << render_colored_fill(pairs, colors, apply_call, buf_indices, "    ");
    return os.str();
}

ExactDerivativeCode emit_exact_derivative_code(const GraphFunction &F,
                                               const std::string &wrt,
                                               const std::string &direction_name,
                                               const std::string &lambda_name,
                                               const std::string &value_name,
                                               const std::string &jvp_name,
                                               const std::string &hvp_name) {
    auto plan = exact_derivative_plan(F, wrt, direction_name, lambda_name);

    ExactDerivativeCode code;
    code.value = to_c(F, value_name);
    code.jvp = to_c(plan.jvp, jvp_name);
    code.hvp = to_c(plan.hvp, hvp_name);
    code.jacobian = to_sparse_jacobian_c(F, wrt, plan.jacobian_sparsity, jvp_name + "_sparse", direction_name);
    code.hessian = to_sparse_hessian_c(plan.hvp, direction_name, plan.hessian_sparsity, hvp_name + "_sparse");
    code.jacobian_sparsity = std::move(plan.jacobian_sparsity);
    code.hessian_sparsity = std::move(plan.hessian_sparsity);
    code.jacobian_colors = std::move(plan.jacobian_colors);
    code.hessian_colors = std::move(plan.hessian_colors);
    code.jacobian_color_count = plan.jacobian_color_count;
    code.hessian_color_count = plan.hessian_color_count;
    code.value_bytes = code.value.size();
    code.jvp_bytes = code.jvp.size();
    code.hvp_bytes = code.hvp.size();
    code.jacobian_bytes = code.jacobian.size();
    code.hessian_bytes = code.hessian.size();
    return code;
}

} // namespace ad
