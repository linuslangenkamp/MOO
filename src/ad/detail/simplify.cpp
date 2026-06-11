// SPDX-License-Identifier: LGPL-3.0-or-later
#include "simplify.h"

#include "../function.h"
#include "../optimize.h"
#include "function_core.h"
#include "mapping.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>

namespace ad {
namespace detail {
namespace {

void require_valid(const Expr &expr) {
    if (!expr.valid()) {
        throw std::runtime_error("invalid scalar expression");
    }
}

void require_valid(const Vec &vec) {
    if (!vec.valid()) {
        throw std::runtime_error("invalid vector expression");
    }
}

void require_same_size(const Vec &lhs, const Vec &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("vector sizes differ");
    }
}

void validate_indices_in_range(const std::vector<int> &indices, int upper_bound, const char *message) {
    for (int index : indices) {
        if (index < 0 || index >= upper_bound) {
            throw std::runtime_error(message);
        }
    }
}

bool scalar_constant_value(const std::shared_ptr<const ScalarNode> &node, double *value = nullptr) {
    if (!node || node->kind != GraphNodeKind::ScalarConstant) {
        return false;
    }
    if (value) {
        *value = node->value;
    }
    return true;
}

bool scalar_constant_value(const Expr &expr, double *value = nullptr) {
    return scalar_constant_value(scalar_node(expr), value);
}

bool scalar_is_zero(const Expr &expr) {
    double value = 0.0;
    return scalar_constant_value(expr, &value) && value == 0.0;
}

bool scalar_is_one(const Expr &expr) {
    double value = 0.0;
    return scalar_constant_value(expr, &value) && value == 1.0;
}

bool vec_constant_values(const Vec &vec, const std::vector<double> **values = nullptr) {
    const auto &node = vec_node(vec);
    if (!node || node->kind != GraphNodeKind::VectorConstant) {
        return false;
    }
    if (values) {
        *values = &node->constants;
    }
    return true;
}

bool values_all_equal(const std::vector<double> &values, double expected) {
    return std::all_of(values.begin(), values.end(), [expected](double value) {
        return value == expected;
    });
}

bool vec_is_zero(const Vec &vec) {
    const std::vector<double> *values = nullptr;
    return vec_constant_values(vec, &values) && values_all_equal(*values, 0.0);
}

bool dense_is_zero(const DenseMatrix &matrix) {
    return std::all_of(matrix.values.begin(), matrix.values.end(), [](double value) {
        return value == 0.0;
    });
}

bool dense_is_identity(const DenseMatrix &matrix) {
    if (matrix.rows != matrix.cols) {
        return false;
    }
    for (int row = 0; row < matrix.rows; ++row) {
        for (int col = 0; col < matrix.cols; ++col) {
            const double expected = row == col ? 1.0 : 0.0;
            if (matrix(row, col) != expected) {
                return false;
            }
        }
    }
    return true;
}

bool sparse_is_zero(const SparseMatrix &matrix) {
    return std::all_of(matrix.values.begin(), matrix.values.end(), [](double value) {
        return value == 0.0;
    });
}

double eval_scalar_unary(ScalarUnaryOp op, double value, double exponent) {
    switch (op) {
        case ScalarUnaryOp::Neg:
            return -value;
        case ScalarUnaryOp::Sin:
            return std::sin(value);
        case ScalarUnaryOp::Cos:
            return std::cos(value);
        case ScalarUnaryOp::Tan:
            return std::tan(value);
        case ScalarUnaryOp::Exp:
            return std::exp(value);
        case ScalarUnaryOp::Log:
            return std::log(value);
        case ScalarUnaryOp::PowConst:
            return std::pow(value, exponent);
        case ScalarUnaryOp::Abs:
            return std::fabs(value);
        case ScalarUnaryOp::Sqrt:
            return std::sqrt(value);
        case ScalarUnaryOp::Asin:
            return std::asin(value);
        case ScalarUnaryOp::Acos:
            return std::acos(value);
        case ScalarUnaryOp::Atan:
            return std::atan(value);
        case ScalarUnaryOp::Sinh:
            return std::sinh(value);
        case ScalarUnaryOp::Cosh:
            return std::cosh(value);
        case ScalarUnaryOp::Tanh:
            return std::tanh(value);
        case ScalarUnaryOp::Log10:
            return std::log10(value);
        case ScalarUnaryOp::Sigmoid:
            return 1.0 / (1.0 + std::exp(-value));
    }
    throw std::runtime_error("unsupported scalar unary op");
}

double eval_scalar_binary(ScalarBinaryOp op, double lhs, double rhs) {
    switch (op) {
        case ScalarBinaryOp::Add:
            return lhs + rhs;
        case ScalarBinaryOp::Sub:
            return lhs - rhs;
        case ScalarBinaryOp::Mul:
            return lhs * rhs;
        case ScalarBinaryOp::Div:
            return lhs / rhs;
        case ScalarBinaryOp::Pow:
            return std::pow(lhs, rhs);
        case ScalarBinaryOp::Min:
            return std::fmin(lhs, rhs);
        case ScalarBinaryOp::Max:
            return std::fmax(lhs, rhs);
    }
    throw std::runtime_error("unsupported scalar binary op");
}

double eval_vec_unary(VecUnaryOp op, double value) {
    switch (op) {
        case VecUnaryOp::Neg:
            return -value;
        case VecUnaryOp::Sin:
            return std::sin(value);
        case VecUnaryOp::Cos:
            return std::cos(value);
        case VecUnaryOp::Tan:
            return std::tan(value);
        case VecUnaryOp::Exp:
            return std::exp(value);
        case VecUnaryOp::Log:
            return std::log(value);
        case VecUnaryOp::Sigmoid:
            return 1.0 / (1.0 + std::exp(-value));
        case VecUnaryOp::Abs:
            return std::fabs(value);
        case VecUnaryOp::Sqrt:
            return std::sqrt(value);
        case VecUnaryOp::Asin:
            return std::asin(value);
        case VecUnaryOp::Acos:
            return std::acos(value);
        case VecUnaryOp::Atan:
            return std::atan(value);
        case VecUnaryOp::Sinh:
            return std::sinh(value);
        case VecUnaryOp::Cosh:
            return std::cosh(value);
        case VecUnaryOp::Tanh:
            return std::tanh(value);
        case VecUnaryOp::Log10:
            return std::log10(value);
    }
    throw std::runtime_error("unsupported vector unary op");
}

double eval_vec_binary(VecBinaryOp op, double lhs, double rhs) {
    switch (op) {
        case VecBinaryOp::Add:
            return lhs + rhs;
        case VecBinaryOp::Sub:
            return lhs - rhs;
        case VecBinaryOp::Mul:
            return lhs * rhs;
        case VecBinaryOp::Div:
            return lhs / rhs;
        case VecBinaryOp::Pow:
            return std::pow(lhs, rhs);
        case VecBinaryOp::Min:
            return std::fmin(lhs, rhs);
        case VecBinaryOp::Max:
            return std::fmax(lhs, rhs);
    }
    throw std::runtime_error("unsupported vector binary op");
}

bool identity_indices(const std::vector<int> &indices, int size) {
    if (static_cast<int>(indices.size()) != size) {
        return false;
    }
    for (int i = 0; i < size; ++i) {
        if (indices[static_cast<std::size_t>(i)] != i) {
            return false;
        }
    }
    return true;
}

struct SimplifyContext {
    std::map<NodeId, Expr> scalar_cache;
    std::map<NodeId, Vec> vec_cache;
};

Expr simplify_scalar_node(const std::shared_ptr<const ScalarNode> &node, SimplifyContext &context);
Vec simplify_vec_node(const std::shared_ptr<const VecNode> &node, SimplifyContext &context);

Expr simplify_scalar_node(const std::shared_ptr<const ScalarNode> &node, SimplifyContext &context) {
    if (!node) {
        throw std::runtime_error("cannot simplify invalid scalar node");
    }

    const auto cached = context.scalar_cache.find(node->id);
    if (cached != context.scalar_cache.end()) {
        return cached->second;
    }

    Expr result;
    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
        case GraphNodeKind::ScalarVariable:
        case GraphNodeKind::ScalarParameter:
            result = expr_from_node(node);
            break;
        case GraphNodeKind::ScalarUnary: {
            const Expr child = simplify_scalar_node(node->lhs, context);
            result = node->unary == ScalarUnaryOp::PowConst ? make_scalar_pow_const(child, node->value)
                                                            : make_scalar_unary(node->unary, child);
            break;
        }
        case GraphNodeKind::ScalarBinary:
            result = make_scalar_binary(node->binary,
                                        simplify_scalar_node(node->lhs, context),
                                        simplify_scalar_node(node->rhs, context));
            break;
        case GraphNodeKind::VectorElement:
            result = make_vector_element(simplify_vec_node(node->vec, context), node->index);
            break;
        case GraphNodeKind::Sum:
            result = make_sum(simplify_vec_node(node->vec, context));
            break;
        case GraphNodeKind::Dot:
            result = make_dot(simplify_vec_node(node->vec_lhs, context),
                              simplify_vec_node(node->vec_rhs, context));
            break;
        default:
            throw std::runtime_error("unsupported scalar node while simplifying");
    }

    context.scalar_cache.emplace(node->id, result);
    return result;
}

Vec simplify_vec_node(const std::shared_ptr<const VecNode> &node, SimplifyContext &context) {
    if (!node) {
        throw std::runtime_error("cannot simplify invalid vector node");
    }

    const auto cached = context.vec_cache.find(node->id);
    if (cached != context.vec_cache.end()) {
        return cached->second;
    }

    Vec result;
    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
        case GraphNodeKind::VectorParameter:
        case GraphNodeKind::VectorConstant:
            result = vec_from_node(node);
            break;
        case GraphNodeKind::VectorFromElements: {
            std::vector<Expr> elements;
            elements.reserve(node->elements.size());
            for (const auto &element : node->elements) {
                elements.push_back(simplify_scalar_node(element, context));
            }
            result = make_vec_from_elements(std::move(elements));
            break;
        }
        case GraphNodeKind::VectorUnary:
            result = make_vec_unary(node->unary, simplify_vec_node(node->lhs, context));
            break;
        case GraphNodeKind::VectorBinary:
            result = make_vec_binary(node->binary,
                                     simplify_vec_node(node->lhs, context),
                                     simplify_vec_node(node->rhs, context));
            break;
        case GraphNodeKind::VectorScale:
            result = make_vec_scale(simplify_scalar_node(node->scale, context),
                                    simplify_vec_node(node->lhs, context));
            break;
        case GraphNodeKind::DenseMatVec:
            result = make_dense_matvec(node->dense, simplify_vec_node(node->lhs, context));
            break;
        case GraphNodeKind::SparseMatVec:
            result = make_sparse_matvec(node->sparse, simplify_vec_node(node->lhs, context));
            break;
        case GraphNodeKind::Slice:
            result = make_slice(simplify_vec_node(node->lhs, context), node->start, node->size);
            break;
        case GraphNodeKind::ScatterSlice:
            result = make_scatter_slice(simplify_vec_node(node->lhs, context), node->start, node->size);
            break;
        case GraphNodeKind::Gather:
            result = make_gather(simplify_vec_node(node->lhs, context), node->indices);
            break;
        case GraphNodeKind::ScatterAdd:
            result = make_scatter_add(simplify_vec_node(node->lhs, context), node->indices, node->size);
            break;
        case GraphNodeKind::Concat:
            result = make_concat(simplify_vec_node(node->lhs, context),
                                 simplify_vec_node(node->rhs, context));
            break;
        case GraphNodeKind::FunctionCall: {
            auto rebuilt = std::make_shared<VecNode>();
            rebuilt->kind = GraphNodeKind::FunctionCall;
            rebuilt->size = node->size;
            rebuilt->function = node->function;
            rebuilt->arguments.reserve(node->arguments.size());
            for (const auto &argument : node->arguments) {
                rebuilt->arguments.push_back(vec_node(simplify_vec_node(argument, context)));
            }
            result = vec_from_node(rebuilt);
            break;
        }
        case GraphNodeKind::MappedFunctionCall: {
            std::vector<MappedBindingNode> bindings = node->mapped_bindings;
            for (MappedBindingNode &binding : bindings) {
                binding.source = vec_node(simplify_vec_node(binding.source, context));
            }
            result = mapped_call_from_bindings(node->function, std::move(bindings), node->reps, node->mapped_output);
            break;
        }
        default:
            throw std::runtime_error("unsupported vector node while simplifying");
    }

    context.vec_cache.emplace(node->id, result);
    return result;
}

} // namespace

Expr make_scalar_constant(double value) {
    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::ScalarConstant;
    node->value = value;
    return expr_from_node(node);
}

Expr make_scalar_variable(const std::string &label) {
    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::ScalarVariable;
    node->var = make_var(label, -1, next_symbol_group_id());
    return expr_from_node(node);
}

Expr make_scalar_parameter(const std::string &label) {
    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::ScalarParameter;
    node->param = make_param(label, -1, next_symbol_group_id());
    return expr_from_node(node);
}

Expr make_scalar_unary(ScalarUnaryOp op, const Expr &arg) {
    require_valid(arg);

    if (op == ScalarUnaryOp::Neg) {
        double value = 0.0;
        if (scalar_constant_value(arg, &value)) {
            return make_scalar_constant(-value);
        }
        const auto &arg_node = scalar_node(arg);
        if (arg_node->kind == GraphNodeKind::ScalarUnary && arg_node->unary == ScalarUnaryOp::Neg) {
            return expr_from_node(arg_node->lhs);
        }
    } else {
        double value = 0.0;
        if (scalar_constant_value(arg, &value)) {
            return make_scalar_constant(eval_scalar_unary(op, value, 0.0));
        }
    }

    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::ScalarUnary;
    node->unary = op;
    node->lhs = scalar_node(arg);
    return expr_from_node(node);
}

Expr make_scalar_pow_const(const Expr &arg, double exponent) {
    require_valid(arg);
    if (exponent == 0.0) {
        return make_scalar_constant(1.0);
    }
    if (exponent == 1.0) {
        return arg;
    }

    double value = 0.0;
    if (scalar_constant_value(arg, &value)) {
        return make_scalar_constant(std::pow(value, exponent));
    }

    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::ScalarUnary;
    node->unary = ScalarUnaryOp::PowConst;
    node->lhs = scalar_node(arg);
    node->value = exponent;
    return expr_from_node(node);
}

Expr make_scalar_binary(ScalarBinaryOp op, const Expr &lhs, const Expr &rhs) {
    require_valid(lhs);
    require_valid(rhs);

    double lhs_value = 0.0;
    double rhs_value = 0.0;
    const bool lhs_constant = scalar_constant_value(lhs, &lhs_value);
    const bool rhs_constant = scalar_constant_value(rhs, &rhs_value);
    if (lhs_constant && rhs_constant) {
        return make_scalar_constant(eval_scalar_binary(op, lhs_value, rhs_value));
    }

    switch (op) {
        case ScalarBinaryOp::Add:
            if (lhs_constant && lhs_value == 0.0) {
                return rhs;
            }
            if (rhs_constant && rhs_value == 0.0) {
                return lhs;
            }
            break;
        case ScalarBinaryOp::Sub:
            if (rhs_constant && rhs_value == 0.0) {
                return lhs;
            }
            if (lhs_constant && lhs_value == 0.0) {
                return make_scalar_unary(ScalarUnaryOp::Neg, rhs);
            }
            break;
        case ScalarBinaryOp::Mul:
            if ((lhs_constant && lhs_value == 0.0) || (rhs_constant && rhs_value == 0.0)) {
                return make_scalar_constant(0.0);
            }
            if (lhs_constant && lhs_value == 1.0) {
                return rhs;
            }
            if (rhs_constant && rhs_value == 1.0) {
                return lhs;
            }
            if (lhs_constant && lhs_value == -1.0) {
                return make_scalar_unary(ScalarUnaryOp::Neg, rhs);
            }
            if (rhs_constant && rhs_value == -1.0) {
                return make_scalar_unary(ScalarUnaryOp::Neg, lhs);
            }
            break;
        case ScalarBinaryOp::Div:
            if (lhs_constant && lhs_value == 0.0) {
                return make_scalar_constant(0.0);
            }
            if (rhs_constant && rhs_value == 1.0) {
                return lhs;
            }
            if (rhs_constant && rhs_value == -1.0) {
                return make_scalar_unary(ScalarUnaryOp::Neg, lhs);
            }
            break;
        case ScalarBinaryOp::Pow:
            if (rhs_constant && rhs_value == 0.0) {
                return make_scalar_constant(1.0);
            }
            if (rhs_constant && rhs_value == 1.0) {
                return lhs;
            }
            if (lhs_constant && lhs_value == 1.0) {
                return make_scalar_constant(1.0);
            }
            break;
        case ScalarBinaryOp::Min:
            if (lhs.node_id() == rhs.node_id()) {
                return lhs;
            }
            break;
        case ScalarBinaryOp::Max:
            if (lhs.node_id() == rhs.node_id()) {
                return lhs;
            }
            break;
    }

    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::ScalarBinary;
    node->binary = op;
    node->lhs = scalar_node(lhs);
    node->rhs = scalar_node(rhs);
    return expr_from_node(node);
}

Expr make_vector_element(const Vec &vec, int index) {
    require_valid(vec);
    if (index < 0 || index >= vec.size()) {
        throw std::out_of_range("vector element index out of range");
    }

    const auto &node = vec_node(vec);
    if (node->kind == GraphNodeKind::VectorConstant) {
        return make_scalar_constant(node->constants[static_cast<std::size_t>(index)]);
    }
    if (node->kind == GraphNodeKind::VectorFromElements) {
        return expr_from_node(node->elements[static_cast<std::size_t>(index)]);
    }
    if (node->kind == GraphNodeKind::Slice) {
        return make_vector_element(vec_from_node(node->lhs), node->start + index);
    }
    if (node->kind == GraphNodeKind::Gather) {
        return make_vector_element(vec_from_node(node->lhs), node->indices[static_cast<std::size_t>(index)]);
    }

    auto out = std::make_shared<ScalarNode>();
    out->kind = GraphNodeKind::VectorElement;
    out->vec = node;
    out->index = index;
    return expr_from_node(out);
}

Expr make_sum(const Vec &values) {
    require_valid(values);
    if (vec_is_zero(values)) {
        return make_scalar_constant(0.0);
    }

    const std::vector<double> *constants = nullptr;
    if (vec_constant_values(values, &constants)) {
        double total = 0.0;
        for (double value : *constants) {
            total += value;
        }
        return make_scalar_constant(total);
    }

    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::Sum;
    node->vec = vec_node(values);
    return expr_from_node(node);
}

Expr make_dot(const Vec &lhs, const Vec &rhs) {
    require_same_size(lhs, rhs);
    if (vec_is_zero(lhs) || vec_is_zero(rhs)) {
        return make_scalar_constant(0.0);
    }

    const std::vector<double> *lhs_constants = nullptr;
    const std::vector<double> *rhs_constants = nullptr;
    if (vec_constant_values(lhs, &lhs_constants) && vec_constant_values(rhs, &rhs_constants)) {
        double total = 0.0;
        for (int i = 0; i < lhs.size(); ++i) {
            total += (*lhs_constants)[static_cast<std::size_t>(i)] * (*rhs_constants)[static_cast<std::size_t>(i)];
        }
        return make_scalar_constant(total);
    }

    auto node = std::make_shared<ScalarNode>();
    node->kind = GraphNodeKind::Dot;
    node->vec_lhs = vec_node(lhs);
    node->vec_rhs = vec_node(rhs);
    return expr_from_node(node);
}

Vec make_vec_from_elements(std::vector<Expr> elements) {
    std::vector<double> constants;
    constants.reserve(elements.size());
    bool all_constant = true;
    for (const Expr &element : elements) {
        if (!element.valid()) {
            throw std::runtime_error("invalid scalar element in vector expression");
        }
        double value = 0.0;
        if (scalar_constant_value(element, &value)) {
            constants.push_back(value);
        } else {
            all_constant = false;
        }
    }
    if (all_constant) {
        return make_vec_constant(std::move(constants));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorFromElements;
    node->size = static_cast<int>(elements.size());
    node->elements.reserve(elements.size());
    for (const Expr &element : elements) {
        node->elements.push_back(scalar_node(element));
    }
    return vec_from_node(node);
}

Vec make_vec_variable(const std::string &label, int size) {
    if (size < 0) {
        throw std::runtime_error("vector variable size must be non-negative");
    }
    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorVariable;
    node->size = size;
    node->label = label;
    node->vars.reserve(static_cast<std::size_t>(size));
    const SymbolGroupId group_id = next_symbol_group_id();
    for (int i = 0; i < size; ++i) {
        node->vars.push_back(make_var(label, i, group_id));
    }
    return vec_from_node(node);
}

Vec make_vec_parameter(const std::string &label, int size) {
    if (size < 0) {
        throw std::runtime_error("vector parameter size must be non-negative");
    }
    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorParameter;
    node->size = size;
    node->label = label;
    node->params.reserve(static_cast<std::size_t>(size));
    const SymbolGroupId group_id = next_symbol_group_id();
    for (int i = 0; i < size; ++i) {
        node->params.push_back(make_param(label, i, group_id));
    }
    return vec_from_node(node);
}

Vec make_vec_constant(std::vector<double> values) {
    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorConstant;
    node->size = static_cast<int>(values.size());
    node->constants = std::move(values);
    return vec_from_node(node);
}

Vec make_vec_unary(VecUnaryOp op, const Vec &arg) {
    require_valid(arg);
    if (op == VecUnaryOp::Neg) {
        const auto &arg_node = vec_node(arg);
        if (arg_node->kind == GraphNodeKind::VectorUnary && arg_node->unary == VecUnaryOp::Neg) {
            return vec_from_node(arg_node->lhs);
        }
    }
    const std::vector<double> *constants = nullptr;
    if (vec_constant_values(arg, &constants)) {
        std::vector<double> values;
        values.reserve(constants->size());
        for (double value : *constants) {
            values.push_back(eval_vec_unary(op, value));
        }
        return make_vec_constant(std::move(values));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorUnary;
    node->size = arg.size();
    node->unary = op;
    node->lhs = vec_node(arg);
    return vec_from_node(node);
}

Vec make_vec_binary(VecBinaryOp op, const Vec &lhs, const Vec &rhs) {
    require_same_size(lhs, rhs);

    const std::vector<double> *lhs_constants = nullptr;
    const std::vector<double> *rhs_constants = nullptr;
    const bool lhs_constant = vec_constant_values(lhs, &lhs_constants);
    const bool rhs_constant = vec_constant_values(rhs, &rhs_constants);
    if (lhs_constant && rhs_constant) {
        std::vector<double> values;
        values.reserve(lhs_constants->size());
        for (int i = 0; i < lhs.size(); ++i) {
            values.push_back(eval_vec_binary(op,
                                             (*lhs_constants)[static_cast<std::size_t>(i)],
                                             (*rhs_constants)[static_cast<std::size_t>(i)]));
        }
        return make_vec_constant(std::move(values));
    }

    switch (op) {
        case VecBinaryOp::Add:
            if (lhs_constant && values_all_equal(*lhs_constants, 0.0)) {
                return rhs;
            }
            if (rhs_constant && values_all_equal(*rhs_constants, 0.0)) {
                return lhs;
            }
            break;
        case VecBinaryOp::Sub:
            if (rhs_constant && values_all_equal(*rhs_constants, 0.0)) {
                return lhs;
            }
            if (lhs_constant && values_all_equal(*lhs_constants, 0.0)) {
                return make_vec_scale(make_scalar_constant(-1.0), rhs);
            }
            break;
        case VecBinaryOp::Mul:
            if ((lhs_constant && values_all_equal(*lhs_constants, 0.0)) ||
                (rhs_constant && values_all_equal(*rhs_constants, 0.0))) {
                return make_vec_constant(std::vector<double>(static_cast<std::size_t>(lhs.size()), 0.0));
            }
            if (lhs_constant && values_all_equal(*lhs_constants, 1.0)) {
                return rhs;
            }
            if (rhs_constant && values_all_equal(*rhs_constants, 1.0)) {
                return lhs;
            }
            break;
        case VecBinaryOp::Div:
            if (lhs_constant && values_all_equal(*lhs_constants, 0.0)) {
                return make_vec_constant(std::vector<double>(static_cast<std::size_t>(lhs.size()), 0.0));
            }
            if (rhs_constant && values_all_equal(*rhs_constants, 1.0)) {
                return lhs;
            }
            break;
        case VecBinaryOp::Pow:
            if (rhs_constant && values_all_equal(*rhs_constants, 0.0)) {
                return make_vec_constant(std::vector<double>(static_cast<std::size_t>(lhs.size()), 1.0));
            }
            if (rhs_constant && values_all_equal(*rhs_constants, 1.0)) {
                return lhs;
            }
            if (lhs_constant && values_all_equal(*lhs_constants, 1.0)) {
                return make_vec_constant(std::vector<double>(static_cast<std::size_t>(lhs.size()), 1.0));
            }
            break;
        case VecBinaryOp::Min:
            if (lhs.node_id() == rhs.node_id()) {
                return lhs;
            }
            break;
        case VecBinaryOp::Max:
            if (lhs.node_id() == rhs.node_id()) {
                return lhs;
            }
            break;
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorBinary;
    node->size = lhs.size();
    node->binary = op;
    node->lhs = vec_node(lhs);
    node->rhs = vec_node(rhs);
    return vec_from_node(node);
}

Vec make_vec_scale(const Expr &scale, const Vec &vec) {
    require_valid(scale);
    require_valid(vec);
    if (scalar_is_zero(scale) || vec_is_zero(vec)) {
        return make_vec_constant(std::vector<double>(static_cast<std::size_t>(vec.size()), 0.0));
    }
    if (scalar_is_one(scale)) {
        return vec;
    }
    double scale_value = 0.0;
    const std::vector<double> *constants = nullptr;
    if (scalar_constant_value(scale, &scale_value) && vec_constant_values(vec, &constants)) {
        std::vector<double> values;
        values.reserve(constants->size());
        for (double value : *constants) {
            values.push_back(scale_value * value);
        }
        return make_vec_constant(std::move(values));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::VectorScale;
    node->size = vec.size();
    node->scale = scalar_node(scale);
    node->lhs = vec_node(vec);
    return vec_from_node(node);
}

Vec make_dense_matvec(const DenseMatrix &matrix, const Vec &rhs) {
    require_valid(rhs);
    if (matrix.cols != rhs.size()) {
        throw std::runtime_error("dense matrix/vector dimensions differ");
    }
    if (matrix.rows == 0) {
        return make_vec_constant({});
    }
    if (vec_is_zero(rhs) || dense_is_zero(matrix)) {
        return make_vec_constant(std::vector<double>(static_cast<std::size_t>(matrix.rows), 0.0));
    }
    if (dense_is_identity(matrix)) {
        return rhs;
    }

    const std::vector<double> *constants = nullptr;
    if (vec_constant_values(rhs, &constants)) {
        std::vector<double> values(static_cast<std::size_t>(matrix.rows), 0.0);
        for (int row = 0; row < matrix.rows; ++row) {
            for (int col = 0; col < matrix.cols; ++col) {
                values[static_cast<std::size_t>(row)] += matrix(row, col) * (*constants)[static_cast<std::size_t>(col)];
            }
        }
        return make_vec_constant(std::move(values));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::DenseMatVec;
    node->size = matrix.rows;
    node->dense = matrix;
    node->lhs = vec_node(rhs);
    return vec_from_node(node);
}

Vec make_sparse_matvec(const SparseMatrix &matrix, const Vec &rhs) {
    require_valid(rhs);
    if (matrix.cols != rhs.size()) {
        throw std::runtime_error("sparse matrix/vector dimensions differ");
    }
    if (matrix.rows == 0) {
        return make_vec_constant({});
    }
    if (vec_is_zero(rhs) || sparse_is_zero(matrix)) {
        return make_vec_constant(std::vector<double>(static_cast<std::size_t>(matrix.rows), 0.0));
    }

    const std::vector<double> *constants = nullptr;
    if (vec_constant_values(rhs, &constants)) {
        std::vector<double> values(static_cast<std::size_t>(matrix.rows), 0.0);
        for (std::size_t k = 0; k < matrix.values.size(); ++k) {
            const double value = matrix.values[k];
            if (value == 0.0) {
                continue;
            }
            values[static_cast<std::size_t>(matrix.row[k])] += value * (*constants)[static_cast<std::size_t>(matrix.col[k])];
        }
        return make_vec_constant(std::move(values));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::SparseMatVec;
    node->size = matrix.rows;
    node->sparse = matrix;
    node->lhs = vec_node(rhs);
    return vec_from_node(node);
}

Vec make_slice(const Vec &source, int start, int length) {
    require_valid(source);
    if (start < 0 || length < 0 || start + length > source.size()) {
        throw std::runtime_error("invalid vector slice");
    }
    if (start == 0 && length == source.size()) {
        return source;
    }

    const auto &source_node = vec_node(source);
    if (source_node->kind == GraphNodeKind::VectorConstant) {
        return make_vec_constant(std::vector<double>(source_node->constants.begin() + start,
                                                     source_node->constants.begin() + start + length));
    }
    if (source_node->kind == GraphNodeKind::Slice) {
        return make_slice(vec_from_node(source_node->lhs), source_node->start + start, length);
    }
    if (source_node->kind == GraphNodeKind::Gather) {
        std::vector<int> indices;
        indices.reserve(static_cast<std::size_t>(length));
        for (int i = 0; i < length; ++i) {
            indices.push_back(source_node->indices[static_cast<std::size_t>(start + i)]);
        }
        return make_gather(vec_from_node(source_node->lhs), std::move(indices));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::Slice;
    node->size = length;
    node->lhs = source_node;
    node->start = start;
    return vec_from_node(node);
}

Vec make_scatter_slice(const Vec &values, int start, int output_size) {
    require_valid(values);
    if (start < 0 || output_size < 0 || start + values.size() > output_size) {
        throw std::runtime_error("invalid scatter-slice bounds");
    }
    if (start == 0 && values.size() == output_size) {
        return values;
    }
    if (vec_is_zero(values)) {
        return make_vec_constant(std::vector<double>(static_cast<std::size_t>(output_size), 0.0));
    }

    const std::vector<double> *constants = nullptr;
    if (vec_constant_values(values, &constants)) {
        std::vector<double> out(static_cast<std::size_t>(output_size), 0.0);
        for (int i = 0; i < values.size(); ++i) {
            out[static_cast<std::size_t>(start + i)] = (*constants)[static_cast<std::size_t>(i)];
        }
        return make_vec_constant(std::move(out));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::ScatterSlice;
    node->size = output_size;
    node->start = start;
    node->lhs = vec_node(values);
    return vec_from_node(node);
}

Vec make_concat(const Vec &lhs, const Vec &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    if (lhs.size() == 0) {
        return rhs;
    }
    if (rhs.size() == 0) {
        return lhs;
    }

    const std::vector<double> *lhs_constants = nullptr;
    const std::vector<double> *rhs_constants = nullptr;
    if (vec_constant_values(lhs, &lhs_constants) && vec_constant_values(rhs, &rhs_constants)) {
        std::vector<double> values;
        values.reserve(lhs_constants->size() + rhs_constants->size());
        values.insert(values.end(), lhs_constants->begin(), lhs_constants->end());
        values.insert(values.end(), rhs_constants->begin(), rhs_constants->end());
        return make_vec_constant(std::move(values));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::Concat;
    node->size = lhs.size() + rhs.size();
    node->lhs = vec_node(lhs);
    node->rhs = vec_node(rhs);
    return vec_from_node(node);
}

Vec make_gather(const Vec &source, std::vector<int> indices) {
    require_valid(source);
    validate_indices_in_range(indices, source.size(), "gather index is out of source vector range");
    if (identity_indices(indices, source.size())) {
        return source;
    }

    const auto &source_node = vec_node(source);
    if (source_node->kind == GraphNodeKind::VectorConstant) {
        std::vector<double> values;
        values.reserve(indices.size());
        for (int index : indices) {
            values.push_back(source_node->constants[static_cast<std::size_t>(index)]);
        }
        return make_vec_constant(std::move(values));
    }
    if (source_node->kind == GraphNodeKind::Gather) {
        std::vector<int> composed;
        composed.reserve(indices.size());
        for (int index : indices) {
            composed.push_back(source_node->indices[static_cast<std::size_t>(index)]);
        }
        return make_gather(vec_from_node(source_node->lhs), std::move(composed));
    }
    if (source_node->kind == GraphNodeKind::Slice) {
        std::vector<int> composed;
        composed.reserve(indices.size());
        for (int index : indices) {
            composed.push_back(source_node->start + index);
        }
        return make_gather(vec_from_node(source_node->lhs), std::move(composed));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::Gather;
    node->size = static_cast<int>(indices.size());
    node->lhs = source_node;
    node->indices = std::move(indices);
    return vec_from_node(node);
}

Vec make_scatter_add(const Vec &values, std::vector<int> indices, int output_size) {
    require_valid(values);
    if (output_size < 0) {
        throw std::runtime_error("scatter_add output size must be non-negative");
    }
    if (static_cast<int>(indices.size()) != values.size()) {
        throw std::runtime_error("scatter_add index count must match value vector size");
    }
    validate_indices_in_range(indices, output_size, "scatter_add index is out of output vector range");
    if (identity_indices(indices, output_size) && values.size() == output_size) {
        return values;
    }
    if (vec_is_zero(values)) {
        return make_vec_constant(std::vector<double>(static_cast<std::size_t>(output_size), 0.0));
    }

    const std::vector<double> *constants = nullptr;
    if (vec_constant_values(values, &constants)) {
        std::vector<double> out(static_cast<std::size_t>(output_size), 0.0);
        for (int i = 0; i < values.size(); ++i) {
            out[static_cast<std::size_t>(indices[static_cast<std::size_t>(i)])] += (*constants)[static_cast<std::size_t>(i)];
        }
        return make_vec_constant(std::move(out));
    }

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::ScatterAdd;
    node->size = output_size;
    node->lhs = vec_node(values);
    node->indices = std::move(indices);
    return vec_from_node(node);
}

Expr simplify_expr(const Expr &expr) {
    require_valid(expr);
    SimplifyContext context;
    return simplify_scalar_node(scalar_node(expr), context);
}

Vec simplify_vec(const Vec &vec) {
    require_valid(vec);
    SimplifyContext context;
    return simplify_vec_node(vec_node(vec), context);
}

} // namespace detail

Expr optimize(const Expr &expr) {
    return detail::simplify_expr(expr);
}

Vec optimize(const Vec &vec) {
    return detail::simplify_vec(vec);
}

Function optimize(const Function &function) {
    if (!function.valid()) {
        throw std::runtime_error("cannot optimize invalid function");
    }
    return Function(function.inputs(), function.parameters(), optimize(function.output()));
}

} // namespace ad
