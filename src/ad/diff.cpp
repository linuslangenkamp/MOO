// SPDX-License-Identifier: LGPL-3.0-or-later
#include "expr.h"

#include "function.h"
#include "vec.h"
#include "detail/function_core.h"
#include "detail/mapping.h"
#include "detail/node.h"
#include "detail/simplify.h"

#include <limits>
#include <map>
#include <cmath>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ad {
namespace {

struct ForwardContext {
    Vars wrt;
    Vec seed;
    std::map<VarId, int> seed_index;
};

struct ReverseContext {
    Vars wrt;
    std::map<VarId, int> adjoint_index;
};

void require_valid_seed(const Vars &wrt, const Vec &seed) {
    if (!seed.valid()) {
        throw std::runtime_error("forward seed must be a valid vector expression");
    }
    if (seed.size() != wrt.size()) {
        throw std::runtime_error("forward seed size must match differentiation variable layout");
    }
}

ForwardContext make_context(const Vars &wrt, const Vec &seed) {
    require_valid_seed(wrt, seed);

    ForwardContext context;
    context.wrt = wrt;
    context.seed = seed;

    std::set<VarId> seen;
    for (int i = 0; i < wrt.size(); ++i) {
        const Var &var = wrt[i];
        if (!var.valid()) {
            throw std::runtime_error("forward differentiation variable must be valid");
        }
        if (!seen.insert(var.id()).second) {
            throw std::runtime_error("forward differentiation variable layout contains duplicate variable IDs");
        }
        context.seed_index.emplace(var.id(), i);
    }

    return context;
}

ReverseContext make_reverse_context(const Vars &wrt) {
    ReverseContext context;
    context.wrt = wrt;

    std::set<VarId> seen;
    for (int i = 0; i < wrt.size(); ++i) {
        const Var &var = wrt[i];
        if (!var.valid()) {
            throw std::runtime_error("reverse differentiation variable must be valid");
        }
        if (!seen.insert(var.id()).second) {
            throw std::runtime_error("reverse differentiation variable layout contains duplicate variable IDs");
        }
        context.adjoint_index.emplace(var.id(), i);
    }

    return context;
}

Expr zero_scalar() {
    return constant(0.0);
}

Expr one_scalar() {
    return constant(1.0);
}

Vec zero_vec(int size) {
    return vec_constant(std::vector<double>(static_cast<std::size_t>(size), 0.0));
}

Vec one_vec(int size) {
    return vec_constant(std::vector<double>(static_cast<std::size_t>(size), 1.0));
}

bool is_zero(const Expr &expr) {
    double value = 0.0;
    return expr.is_constant(&value) && value == 0.0;
}

bool is_zero(const Vec &vec_expr) {
    if (!vec_expr.valid() || vec_expr.node_kind() != GraphNodeKind::VectorConstant) {
        return false;
    }
    const auto &node = detail::vec_node(vec_expr);
    for (double value : node->constants) {
        if (value != 0.0) {
            return false;
        }
    }
    return true;
}

Expr add(const Expr &lhs, const Expr &rhs) {
    if (is_zero(lhs)) {
        return rhs;
    }
    if (is_zero(rhs)) {
        return lhs;
    }
    return lhs + rhs;
}

Expr sub(const Expr &lhs, const Expr &rhs) {
    if (is_zero(rhs)) {
        return lhs;
    }
    return lhs - rhs;
}

Expr mul(const Expr &lhs, const Expr &rhs) {
    if (is_zero(lhs) || is_zero(rhs)) {
        return zero_scalar();
    }
    return lhs * rhs;
}

Expr div(const Expr &lhs, const Expr &rhs) {
    if (is_zero(lhs)) {
        return zero_scalar();
    }
    return lhs / rhs;
}

Vec add(const Vec &lhs, const Vec &rhs) {
    if (is_zero(lhs)) {
        return rhs;
    }
    if (is_zero(rhs)) {
        return lhs;
    }
    return lhs + rhs;
}

Vec sub(const Vec &lhs, const Vec &rhs) {
    if (is_zero(rhs)) {
        return lhs;
    }
    return lhs - rhs;
}

Vec mul(const Vec &lhs, const Vec &rhs) {
    if (is_zero(lhs) || is_zero(rhs)) {
        return zero_vec(lhs.size());
    }
    return lhs * rhs;
}

Vec div(const Vec &lhs, const Vec &rhs) {
    if (is_zero(lhs)) {
        return zero_vec(lhs.size());
    }
    return lhs / rhs;
}

Vec scale(const Expr &lhs, const Vec &rhs) {
    if (is_zero(lhs) || is_zero(rhs)) {
        return zero_vec(rhs.size());
    }
    return lhs * rhs;
}

Vec matvec(const DenseMatrix &matrix, const Vec &rhs) {
    if (is_zero(rhs)) {
        return zero_vec(matrix.rows);
    }
    return matrix * rhs;
}

Vec matvec(const SparseMatrix &matrix, const Vec &rhs) {
    if (is_zero(rhs)) {
        return zero_vec(matrix.rows);
    }
    return matrix * rhs;
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
    return SparseMatrix(matrix.cols, matrix.rows, matrix.col, matrix.row, matrix.values);
}

Expr as_expr(const std::shared_ptr<const detail::ScalarNode> &node) {
    return detail::expr_from_node(node);
}

Vec as_vec(const std::shared_ptr<const detail::VecNode> &node) {
    return detail::vec_from_node(node);
}

int checked_int(long long value, const char *message) {
    if (value < static_cast<long long>(std::numeric_limits<int>::min()) ||
        value > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(message);
    }
    return static_cast<int>(value);
}

std::size_t checked_table_size(int reps, int local_size, const char *message) {
    if (reps < 0 || local_size < 0) {
        throw std::runtime_error(message);
    }
    const long long value = static_cast<long long>(reps) * static_cast<long long>(local_size);
    checked_int(value, message);
    return static_cast<std::size_t>(value);
}

Expr seed_component(const ForwardContext &context, int index) {
    return context.seed[index];
}

Expr forward_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node, const ForwardContext &context);
Vec forward_vec_node(const std::shared_ptr<const detail::VecNode> &node, const ForwardContext &context);
Vec reverse_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node, const Expr &adjoint, const ReverseContext &context);
Vec reverse_vec_node(const std::shared_ptr<const detail::VecNode> &node, const Vec &adjoint, const ReverseContext &context);

Vec scatter_slice(const Vec &values, int start, int output_size) {
    return detail::make_scatter_slice(values, start, output_size);
}

std::vector<Vec> function_call_arguments(const std::shared_ptr<const detail::VecNode> &node) {
    std::vector<Vec> arguments;
    arguments.reserve(node->arguments.size());
    for (const auto &argument : node->arguments) {
        arguments.push_back(as_vec(argument));
    }
    return arguments;
}

detail::MappedBindingNode mapped_binding(const Vec &local_input,
                                         const std::shared_ptr<const detail::VecNode> &source,
                                         int reps,
                                         int local_size,
                                         const std::vector<int> &indices,
                                         MapKind map_kind = MapKind::ExplicitIndices,
                                         int base = 0,
                                         int rep_stride = 0,
                                         int component_stride = 1,
                                         int shift = 0,
                                         std::vector<int> offsets = {}) {
    detail::MappedBindingNode binding;
    binding.local_input = detail::vec_node(local_input);
    binding.source = source;
    binding.reps = reps;
    binding.local_size = local_size;
    binding.indices = indices;
    binding.map_kind = map_kind;
    binding.base = base;
    binding.rep_stride = rep_stride;
    binding.component_stride = component_stride;
    binding.shift = shift;
    binding.offsets = std::move(offsets);
    return binding;
}

detail::MappedBindingNode mapped_binding(const Vec &local_input,
                                         const Vec &source,
                                         int reps,
                                         int local_size,
                                         const std::vector<int> &indices,
                                         MapKind map_kind,
                                         int base,
                                         int rep_stride,
                                         int component_stride,
                                         int shift,
                                         std::vector<int> offsets) {
    return mapped_binding(local_input,
                          detail::vec_node(source),
                          reps,
                          local_size,
                          indices,
                          map_kind,
                          base,
                          rep_stride,
                          component_stride,
                          shift,
                          std::move(offsets));
}

detail::MappedBindingNode mapped_binding_like(const Vec &local_input,
                                              const std::shared_ptr<const detail::VecNode> &source,
                                              const detail::MappedBindingNode &like) {
    return mapped_binding(local_input,
                          source,
                          like.reps,
                          like.local_size,
                          like.indices,
                          like.map_kind,
                          like.base,
                          like.rep_stride,
                          like.component_stride,
                          like.shift,
                          like.offsets);
}

detail::MappedBindingNode mapped_binding_like(const Vec &local_input,
                                              const Vec &source,
                                              const detail::MappedBindingNode &like) {
    return mapped_binding_like(local_input, detail::vec_node(source), like);
}

Vec weighted_repeated_lambda(const Vec &lambda, int reps, int local_output_size, const std::vector<double> &weights) {
    if (static_cast<int>(weights.size()) != reps) {
        throw std::runtime_error("mapped weighted reverse weight count must match reps");
    }
    std::vector<int> rows;
    std::vector<int> cols;
    std::vector<double> values;
    rows.reserve(checked_table_size(reps, local_output_size, "mapped weighted reverse size overflow"));
    cols.reserve(rows.capacity());
    values.reserve(rows.capacity());
    for (int rep = 0; rep < reps; ++rep) {
        const double weight = weights[static_cast<std::size_t>(rep)];
        if (weight == 0.0) {
            continue;
        }
        for (int component = 0; component < local_output_size; ++component) {
            rows.push_back(checked_int(static_cast<long long>(rep) * static_cast<long long>(local_output_size) + component,
                                       "mapped weighted reverse row index overflow"));
            cols.push_back(component);
            values.push_back(weight);
        }
    }
    return SparseMatrix(checked_int(static_cast<long long>(reps) * static_cast<long long>(local_output_size),
                                    "mapped weighted reverse source size overflow"),
                        local_output_size,
                        std::move(rows),
                        std::move(cols),
                        std::move(values)) *
           lambda;
}

std::vector<int> repeated_local_adjoint_indices(int reps, int local_input_size, int function_input_size, int input_offset) {
    std::vector<int> indices;
    indices.reserve(checked_table_size(reps, local_input_size, "mapped reverse local adjoint index table size overflow"));
    for (int rep = 0; rep < reps; ++rep) {
        for (int component = 0; component < local_input_size; ++component) {
            const long long base = static_cast<long long>(rep) * static_cast<long long>(function_input_size);
            indices.push_back(checked_int(base + input_offset + component, "mapped reverse local adjoint index overflow"));
        }
    }
    return indices;
}

Vec forward_vector_variable(const detail::VecNode &node, const ForwardContext &context) {
    if (node.size == 0) {
        return zero_vec(0);
    }

    std::vector<int> indices(static_cast<std::size_t>(node.size), -1);
    int matched = 0;
    for (int i = 0; i < node.size; ++i) {
        const auto found = context.seed_index.find(node.vars[static_cast<std::size_t>(i)].id());
        if (found != context.seed_index.end()) {
            indices[static_cast<std::size_t>(i)] = found->second;
            ++matched;
        }
    }

    if (matched == 0) {
        return zero_vec(node.size);
    }

    if (matched == node.size) {
        const int start = indices.front();
        bool contiguous = true;
        for (int i = 0; i < node.size; ++i) {
            if (indices[static_cast<std::size_t>(i)] != start + i) {
                contiguous = false;
                break;
            }
        }
        if (contiguous) {
            if (start == 0 && node.size == context.seed.size()) {
                return context.seed;
            }
            return context.seed.slice(start, node.size);
        }
    }

    std::vector<Expr> elements;
    elements.reserve(static_cast<std::size_t>(node.size));
    for (int index : indices) {
        elements.push_back(index >= 0 ? seed_component(context, index) : zero_scalar());
    }
    return vec(std::move(elements));
}

Expr forward_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node, const ForwardContext &context) {
    if (!node) {
        throw std::runtime_error("invalid scalar graph node while building forward derivative");
    }

    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
        case GraphNodeKind::ScalarParameter:
            return zero_scalar();
        case GraphNodeKind::ScalarVariable: {
            const auto found = context.seed_index.find(node->var.id());
            return found == context.seed_index.end() ? zero_scalar() : seed_component(context, found->second);
        }
        case GraphNodeKind::ScalarUnary: {
            const Expr arg = as_expr(node->lhs);
            const Expr darg = forward_scalar_node(node->lhs, context);
            switch (node->unary) {
                case detail::ScalarUnaryOp::Neg:
                    return -darg;
                case detail::ScalarUnaryOp::Sin:
                    return mul(cos(arg), darg);
                case detail::ScalarUnaryOp::Cos:
                    return mul(-sin(arg), darg);
                case detail::ScalarUnaryOp::Tan:
                    return mul(add(one_scalar(), tan(arg) * tan(arg)), darg);
                case detail::ScalarUnaryOp::Exp:
                    return mul(exp(arg), darg);
                case detail::ScalarUnaryOp::Log:
                    return div(darg, arg);
                case detail::ScalarUnaryOp::PowConst:
                    return node->value == 0.0 ? zero_scalar() : mul(mul(constant(node->value), pow(arg, node->value - 1.0)), darg);
                case detail::ScalarUnaryOp::Abs:
                    throw std::runtime_error("forward derivative of nonsmooth abs is not supported");
                case detail::ScalarUnaryOp::Sqrt:
                    return div(darg, mul(constant(2.0), sqrt(arg)));
                case detail::ScalarUnaryOp::Asin:
                    return div(darg, sqrt(sub(one_scalar(), arg * arg)));
                case detail::ScalarUnaryOp::Acos:
                    return -div(darg, sqrt(sub(one_scalar(), arg * arg)));
                case detail::ScalarUnaryOp::Atan:
                    return div(darg, add(one_scalar(), arg * arg));
                case detail::ScalarUnaryOp::Sinh:
                    return mul(cosh(arg), darg);
                case detail::ScalarUnaryOp::Cosh:
                    return mul(sinh(arg), darg);
                case detail::ScalarUnaryOp::Tanh: {
                    const Expr t = tanh(arg);
                    return mul(sub(one_scalar(), t * t), darg);
                }
                case detail::ScalarUnaryOp::Log10:
                    return div(darg, mul(arg, constant(std::log(10.0))));
                case detail::ScalarUnaryOp::Sigmoid: {
                    const Expr sig = sigmoid(arg);
                    return mul(mul(sig, sub(one_scalar(), sig)), darg);
                }
            }
            break;
        }
        case GraphNodeKind::ScalarBinary: {
            const Expr lhs = as_expr(node->lhs);
            const Expr rhs = as_expr(node->rhs);
            const Expr dlhs = forward_scalar_node(node->lhs, context);
            const Expr drhs = forward_scalar_node(node->rhs, context);
            switch (node->binary) {
                case detail::ScalarBinaryOp::Add:
                    return add(dlhs, drhs);
                case detail::ScalarBinaryOp::Sub:
                    return sub(dlhs, drhs);
                case detail::ScalarBinaryOp::Mul:
                    return add(mul(dlhs, rhs), mul(lhs, drhs));
                case detail::ScalarBinaryOp::Div:
                    return div(sub(mul(dlhs, rhs), mul(lhs, drhs)), rhs * rhs);
                case detail::ScalarBinaryOp::Pow: {
                    const Expr value = pow(lhs, rhs);
                    return mul(value, add(mul(drhs, log(lhs)), div(mul(rhs, dlhs), lhs)));
                }
                case detail::ScalarBinaryOp::Min:
                    throw std::runtime_error("forward derivative of nonsmooth min is not supported");
                case detail::ScalarBinaryOp::Max:
                    throw std::runtime_error("forward derivative of nonsmooth max is not supported");
            }
            break;
        }
        case GraphNodeKind::VectorElement:
            return forward_vec_node(node->vec, context)[node->index];
        case GraphNodeKind::Sum:
            return sum(forward_vec_node(node->vec, context));
        case GraphNodeKind::Dot: {
            const Vec lhs = as_vec(node->vec_lhs);
            const Vec rhs = as_vec(node->vec_rhs);
            const Vec dlhs = forward_vec_node(node->vec_lhs, context);
            const Vec drhs = forward_vec_node(node->vec_rhs, context);
            return add(is_zero(dlhs) ? zero_scalar() : dot(dlhs, rhs), is_zero(drhs) ? zero_scalar() : dot(lhs, drhs));
        }
        default:
            throw std::runtime_error("unsupported scalar graph node while building forward derivative");
    }

    throw std::runtime_error("unhandled scalar graph node while building forward derivative");
}

Vec forward_through_function_call(const std::shared_ptr<const detail::VecNode> &node, const ForwardContext &context) {
    if (!node->function) {
        throw std::runtime_error("function call node is missing callee while building forward derivative");
    }

    std::vector<Vec> arguments = function_call_arguments(node);
    arguments.reserve(arguments.size() + node->arguments.size());
    for (const auto &argument : node->arguments) {
        arguments.push_back(forward_vec_node(argument, context));
    }
    return detail::function_from_core(node->function).forward_function().call(std::move(arguments));
}

Vec forward_mapped_function_call(const std::shared_ptr<const detail::VecNode> &node, const ForwardContext &context) {
    if (!node->function) {
        throw std::runtime_error("mapped function call node is missing callee while building forward derivative");
    }

    Function transformed = detail::function_from_core(node->function).forward_function();
    const std::shared_ptr<const detail::FunctionCore> &transformed_core = detail::function_core(transformed);
    const int input_count = node->function->info.input_count;
    if (static_cast<int>(node->mapped_bindings.size()) != input_count) {
        throw std::runtime_error("mapped function call binding count does not match callee input count");
    }
    if (transformed_core->info.input_count != input_count * 2) {
        throw std::runtime_error("mapped forward transformed function has unexpected input count");
    }

    std::vector<detail::MappedBindingNode> bindings;
    bindings.reserve(static_cast<std::size_t>(input_count * 2));
    for (int i = 0; i < input_count; ++i) {
        const detail::MappedBindingNode &binding = node->mapped_bindings[static_cast<std::size_t>(i)];
        bindings.push_back(mapped_binding_like(transformed_core->inputs[static_cast<std::size_t>(i)],
                                               binding.source,
                                               binding));
    }
    for (int i = 0; i < input_count; ++i) {
        const detail::MappedBindingNode &binding = node->mapped_bindings[static_cast<std::size_t>(i)];
        Vec source_tangent = forward_vec_node(binding.source, context);
        bindings.push_back(mapped_binding_like(transformed_core->inputs[static_cast<std::size_t>(input_count + i)],
                                               source_tangent,
                                               binding));
    }

    return detail::mapped_call_from_bindings(transformed_core, std::move(bindings), node->reps, node->mapped_output);
}

Vec forward_vec_node(const std::shared_ptr<const detail::VecNode> &node, const ForwardContext &context) {
    if (!node) {
        throw std::runtime_error("invalid vector graph node while building forward derivative");
    }

    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
            return forward_vector_variable(*node, context);
        case GraphNodeKind::VectorParameter:
        case GraphNodeKind::VectorConstant:
            return zero_vec(node->size);
        case GraphNodeKind::VectorFromElements: {
            std::vector<Expr> elements;
            elements.reserve(node->elements.size());
            for (const auto &element : node->elements) {
                elements.push_back(forward_scalar_node(element, context));
            }
            return vec(std::move(elements));
        }
        case GraphNodeKind::VectorUnary: {
            const Vec arg = as_vec(node->lhs);
            const Vec darg = forward_vec_node(node->lhs, context);
            if (is_zero(darg)) {
                return zero_vec(node->size);
            }
            switch (node->unary) {
                case detail::VecUnaryOp::Neg:
                    return -darg;
                case detail::VecUnaryOp::Sin:
                    return mul(cos(arg), darg);
                case detail::VecUnaryOp::Cos:
                    return scale(constant(-1.0), mul(sin(arg), darg));
                case detail::VecUnaryOp::Tan:
                    return mul(one_vec(node->size) + tan(arg) * tan(arg), darg);
                case detail::VecUnaryOp::Exp:
                    return mul(exp(arg), darg);
                case detail::VecUnaryOp::Log:
                    return div(darg, arg);
                case detail::VecUnaryOp::Sigmoid: {
                    const Vec sig = sigmoid(arg);
                    return mul(mul(sig, one_vec(node->size) - sig), darg);
                }
                case detail::VecUnaryOp::Abs:
                    throw std::runtime_error("forward derivative of nonsmooth abs is not supported");
                case detail::VecUnaryOp::Sqrt:
                    return div(darg, scale(constant(2.0), sqrt(arg)));
                case detail::VecUnaryOp::Asin:
                    return div(darg, sqrt(one_vec(node->size) - arg * arg));
                case detail::VecUnaryOp::Acos:
                    return scale(constant(-1.0), div(darg, sqrt(one_vec(node->size) - arg * arg)));
                case detail::VecUnaryOp::Atan:
                    return div(darg, one_vec(node->size) + arg * arg);
                case detail::VecUnaryOp::Sinh:
                    return mul(cosh(arg), darg);
                case detail::VecUnaryOp::Cosh:
                    return mul(sinh(arg), darg);
                case detail::VecUnaryOp::Tanh: {
                    const Vec t = tanh(arg);
                    return mul(one_vec(node->size) - t * t, darg);
                }
                case detail::VecUnaryOp::Log10:
                    return div(darg, scale(constant(std::log(10.0)), arg));
            }
            break;
        }
        case GraphNodeKind::VectorBinary: {
            const Vec lhs = as_vec(node->lhs);
            const Vec rhs = as_vec(node->rhs);
            const Vec dlhs = forward_vec_node(node->lhs, context);
            const Vec drhs = forward_vec_node(node->rhs, context);
            switch (node->binary) {
                case detail::VecBinaryOp::Add:
                    return add(dlhs, drhs);
                case detail::VecBinaryOp::Sub:
                    return sub(dlhs, drhs);
                case detail::VecBinaryOp::Mul:
                    return add(mul(dlhs, rhs), mul(lhs, drhs));
                case detail::VecBinaryOp::Div:
                    return div(sub(mul(dlhs, rhs), mul(lhs, drhs)), rhs * rhs);
                case detail::VecBinaryOp::Pow: {
                    const Vec value = pow(lhs, rhs);
                    return mul(value, add(mul(drhs, log(lhs)), div(mul(rhs, dlhs), lhs)));
                }
                case detail::VecBinaryOp::Min:
                    throw std::runtime_error("forward derivative of nonsmooth min is not supported");
                case detail::VecBinaryOp::Max:
                    throw std::runtime_error("forward derivative of nonsmooth max is not supported");
            }
            break;
        }
        case GraphNodeKind::VectorScale: {
            const Expr scale_expr = as_expr(node->scale);
            const Vec vec_arg = as_vec(node->lhs);
            return add(scale(forward_scalar_node(node->scale, context), vec_arg), scale(scale_expr, forward_vec_node(node->lhs, context)));
        }
        case GraphNodeKind::DenseMatVec:
            return matvec(node->dense, forward_vec_node(node->lhs, context));
        case GraphNodeKind::SparseMatVec:
            return matvec(node->sparse, forward_vec_node(node->lhs, context));
        case GraphNodeKind::Slice:
            return forward_vec_node(node->lhs, context).slice(node->start, node->size);
        case GraphNodeKind::ScatterSlice:
            return scatter_slice(forward_vec_node(node->lhs, context), node->start, node->size);
        case GraphNodeKind::Gather:
            return gather(forward_vec_node(node->lhs, context), node->indices);
        case GraphNodeKind::ScatterAdd:
            return scatter_add(forward_vec_node(node->lhs, context), node->indices, node->size);
        case GraphNodeKind::Concat:
            return concat(forward_vec_node(node->lhs, context), forward_vec_node(node->rhs, context));
        case GraphNodeKind::FunctionCall:
            return forward_through_function_call(node, context);
        case GraphNodeKind::MappedFunctionCall:
            return forward_mapped_function_call(node, context);
        default:
            throw std::runtime_error("unsupported vector graph node while building forward derivative");
    }

    throw std::runtime_error("unhandled vector graph node while building forward derivative");
}

Vec reverse_scalar_variable(const Var &var, const Expr &adjoint, const ReverseContext &context) {
    const auto found = context.adjoint_index.find(var.id());
    if (found == context.adjoint_index.end()) {
        return zero_vec(context.wrt.size());
    }
    return scatter_slice(vec({adjoint}), found->second, context.wrt.size());
}

Vec reverse_vector_variable(const detail::VecNode &node, const Vec &adjoint, const ReverseContext &context) {
    if (node.size == 0 || context.wrt.size() == 0) {
        return zero_vec(context.wrt.size());
    }

    std::vector<int> indices(static_cast<std::size_t>(node.size), -1);
    int matched = 0;
    for (int i = 0; i < node.size; ++i) {
        const auto found = context.adjoint_index.find(node.vars[static_cast<std::size_t>(i)].id());
        if (found != context.adjoint_index.end()) {
            indices[static_cast<std::size_t>(i)] = found->second;
            ++matched;
        }
    }

    if (matched == 0) {
        return zero_vec(context.wrt.size());
    }

    if (matched == node.size) {
        const int start = indices.front();
        bool contiguous = true;
        for (int i = 0; i < node.size; ++i) {
            if (indices[static_cast<std::size_t>(i)] != start + i) {
                contiguous = false;
                break;
            }
        }
        if (contiguous) {
            return scatter_slice(adjoint, start, context.wrt.size());
        }
    }

    Vec total = zero_vec(context.wrt.size());
    for (int i = 0; i < node.size; ++i) {
        const int index = indices[static_cast<std::size_t>(i)];
        if (index >= 0) {
            total = add(total, scatter_slice(vec({adjoint[i]}), index, context.wrt.size()));
        }
    }
    return total;
}

Vec call_reverse_transformed_function(const std::shared_ptr<const detail::VecNode> &node, const Vec &adjoint) {
    if (!node->function) {
        throw std::runtime_error("function call node is missing callee while building reverse derivative");
    }
    if (!adjoint.valid() || adjoint.size() != node->function->info.output_size) {
        throw std::runtime_error("function reverse-call seed size must match callee output size");
    }

    std::vector<Vec> arguments = function_call_arguments(node);
    if (node->function->info.output_size > 0) {
        arguments.push_back(adjoint);
    }
    return detail::function_from_core(node->function).reverse_function().call(std::move(arguments));
}

Vec reverse_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node, const Expr &adjoint, const ReverseContext &context) {
    if (!node) {
        throw std::runtime_error("invalid scalar graph node while building reverse derivative");
    }
    if (!adjoint.valid()) {
        throw std::runtime_error("scalar reverse adjoint must be a valid expression");
    }

    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
        case GraphNodeKind::ScalarParameter:
            return zero_vec(context.wrt.size());
        case GraphNodeKind::ScalarVariable:
            return reverse_scalar_variable(node->var, adjoint, context);
        case GraphNodeKind::ScalarUnary: {
            const Expr arg = as_expr(node->lhs);
            switch (node->unary) {
                case detail::ScalarUnaryOp::Neg:
                    return reverse_scalar_node(node->lhs, -adjoint, context);
                case detail::ScalarUnaryOp::Sin:
                    return reverse_scalar_node(node->lhs, mul(adjoint, cos(arg)), context);
                case detail::ScalarUnaryOp::Cos:
                    return reverse_scalar_node(node->lhs, mul(adjoint, -sin(arg)), context);
                case detail::ScalarUnaryOp::Tan:
                    return reverse_scalar_node(node->lhs, mul(adjoint, add(one_scalar(), tan(arg) * tan(arg))), context);
                case detail::ScalarUnaryOp::Exp:
                    return reverse_scalar_node(node->lhs, mul(adjoint, exp(arg)), context);
                case detail::ScalarUnaryOp::Log:
                    return reverse_scalar_node(node->lhs, div(adjoint, arg), context);
                case detail::ScalarUnaryOp::PowConst:
                    return node->value == 0.0 ? zero_vec(context.wrt.size()) : reverse_scalar_node(node->lhs, mul(adjoint, mul(constant(node->value), pow(arg, node->value - 1.0))), context);
                case detail::ScalarUnaryOp::Abs:
                    throw std::runtime_error("reverse derivative of nonsmooth abs is not supported");
                case detail::ScalarUnaryOp::Sqrt:
                    return reverse_scalar_node(node->lhs, div(adjoint, mul(constant(2.0), sqrt(arg))), context);
                case detail::ScalarUnaryOp::Asin:
                    return reverse_scalar_node(node->lhs, div(adjoint, sqrt(sub(one_scalar(), arg * arg))), context);
                case detail::ScalarUnaryOp::Acos:
                    return reverse_scalar_node(node->lhs, -div(adjoint, sqrt(sub(one_scalar(), arg * arg))), context);
                case detail::ScalarUnaryOp::Atan:
                    return reverse_scalar_node(node->lhs, div(adjoint, add(one_scalar(), arg * arg)), context);
                case detail::ScalarUnaryOp::Sinh:
                    return reverse_scalar_node(node->lhs, mul(adjoint, cosh(arg)), context);
                case detail::ScalarUnaryOp::Cosh:
                    return reverse_scalar_node(node->lhs, mul(adjoint, sinh(arg)), context);
                case detail::ScalarUnaryOp::Tanh: {
                    const Expr t = tanh(arg);
                    return reverse_scalar_node(node->lhs, mul(adjoint, sub(one_scalar(), t * t)), context);
                }
                case detail::ScalarUnaryOp::Log10:
                    return reverse_scalar_node(node->lhs, div(adjoint, mul(arg, constant(std::log(10.0)))), context);
                case detail::ScalarUnaryOp::Sigmoid: {
                    const Expr sig = sigmoid(arg);
                    return reverse_scalar_node(node->lhs, mul(adjoint, mul(sig, sub(one_scalar(), sig))), context);
                }
            }
            break;
        }
        case GraphNodeKind::ScalarBinary: {
            const Expr lhs = as_expr(node->lhs);
            const Expr rhs = as_expr(node->rhs);
            switch (node->binary) {
                case detail::ScalarBinaryOp::Add:
                    return add(reverse_scalar_node(node->lhs, adjoint, context), reverse_scalar_node(node->rhs, adjoint, context));
                case detail::ScalarBinaryOp::Sub:
                    return add(reverse_scalar_node(node->lhs, adjoint, context), reverse_scalar_node(node->rhs, -adjoint, context));
                case detail::ScalarBinaryOp::Mul:
                    return add(reverse_scalar_node(node->lhs, mul(adjoint, rhs), context), reverse_scalar_node(node->rhs, mul(adjoint, lhs), context));
                case detail::ScalarBinaryOp::Div:
                    return add(reverse_scalar_node(node->lhs, div(adjoint, rhs), context),
                               reverse_scalar_node(node->rhs, div(mul(-adjoint, lhs), rhs * rhs), context));
                case detail::ScalarBinaryOp::Pow: {
                    const Expr value = pow(lhs, rhs);
                    return add(reverse_scalar_node(node->lhs, mul(adjoint, div(mul(value, rhs), lhs)), context),
                               reverse_scalar_node(node->rhs, mul(adjoint, mul(value, log(lhs))), context));
                }
                case detail::ScalarBinaryOp::Min:
                    throw std::runtime_error("reverse derivative of nonsmooth min is not supported");
                case detail::ScalarBinaryOp::Max:
                    throw std::runtime_error("reverse derivative of nonsmooth max is not supported");
            }
            break;
        }
        case GraphNodeKind::VectorElement:
            return reverse_vec_node(node->vec, scatter_slice(vec({adjoint}), node->index, node->vec->size), context);
        case GraphNodeKind::Sum:
            return reverse_vec_node(node->vec, scale(adjoint, one_vec(node->vec->size)), context);
        case GraphNodeKind::Dot: {
            const Vec lhs = as_vec(node->vec_lhs);
            const Vec rhs = as_vec(node->vec_rhs);
            return add(reverse_vec_node(node->vec_lhs, scale(adjoint, rhs), context),
                       reverse_vec_node(node->vec_rhs, scale(adjoint, lhs), context));
        }
        default:
            throw std::runtime_error("unsupported scalar graph node while building reverse derivative");
    }

    throw std::runtime_error("unhandled scalar graph node while building reverse derivative");
}

Vec reverse_function_call(const std::shared_ptr<const detail::VecNode> &node, const Vec &adjoint, const ReverseContext &context) {
    const Vec local_adjoint = call_reverse_transformed_function(node, adjoint);
    Vec total = zero_vec(context.wrt.size());
    int offset = 0;
    for (std::size_t i = 0; i < node->arguments.size(); ++i) {
        const int input_size = node->function->info.inputs[i].size;
        const Vec argument_adjoint = local_adjoint.slice(offset, input_size);
        total = add(total, reverse_vec_node(node->arguments[i], argument_adjoint, context));
        offset += input_size;
    }
    return total;
}

Vec reverse_mapped_function_call(const std::shared_ptr<const detail::VecNode> &node, const Vec &adjoint, const ReverseContext &context) {
    if (!node->function) {
        throw std::runtime_error("mapped function call node is missing callee while building reverse derivative");
    }
    if (adjoint.size() != node->size) {
        throw std::runtime_error("mapped reverse adjoint size must match mapped output size");
    }

    Function transformed = detail::function_from_core(node->function).reverse_function();
    const std::shared_ptr<const detail::FunctionCore> &transformed_core = detail::function_core(transformed);
    const int input_count = node->function->info.input_count;
    if (static_cast<int>(node->mapped_bindings.size()) != input_count) {
        throw std::runtime_error("mapped function call binding count does not match callee input count");
    }
    const int expected_transformed_inputs = input_count + (node->function->info.output_size > 0 ? 1 : 0);
    if (transformed_core->info.input_count != expected_transformed_inputs) {
        throw std::runtime_error("mapped reverse transformed function has unexpected input count");
    }

    std::vector<detail::MappedBindingNode> bindings;
    bindings.reserve(static_cast<std::size_t>(expected_transformed_inputs));
    for (int i = 0; i < input_count; ++i) {
        const detail::MappedBindingNode &binding = node->mapped_bindings[static_cast<std::size_t>(i)];
        bindings.push_back(mapped_binding_like(transformed_core->inputs[static_cast<std::size_t>(i)],
                                               binding.source,
                                               binding));
    }
    if (node->function->info.output_size > 0) {
        const int local_output_size = node->function->info.output_size;
        Vec lambda_source = adjoint;
        int lambda_local_size = local_output_size;
        std::vector<int> lambda_indices;
        MapKind lambda_kind = MapKind::Blocks;
        int lambda_rep_stride = local_output_size;

        if (node->mapped_output.mode() == OutputMode::Concat) {
            lambda_indices = {};
            lambda_kind = MapKind::Blocks;
            lambda_rep_stride = local_output_size;
        } else if (node->mapped_output.mode() == OutputMode::Sum) {
            lambda_indices = {};
            lambda_kind = MapKind::Stride;
            lambda_rep_stride = 0;
        } else if (node->mapped_output.mode() == OutputMode::Scatter) {
            lambda_source = gather(adjoint, node->mapped_output.indices());
            lambda_indices = {};
            lambda_kind = MapKind::Blocks;
            lambda_rep_stride = local_output_size;
        } else if (node->mapped_output.mode() == OutputMode::WeightedSum) {
            lambda_source = weighted_repeated_lambda(adjoint, node->reps, local_output_size, node->mapped_output.weights());
            lambda_indices = {};
            lambda_kind = MapKind::Blocks;
            lambda_rep_stride = local_output_size;
        }

        bindings.push_back(mapped_binding(transformed_core->inputs[static_cast<std::size_t>(input_count)],
                                          lambda_source,
                                          node->reps,
                                          lambda_local_size,
                                          lambda_indices,
                                          lambda_kind,
                                          0,
                                          lambda_rep_stride,
                                          1,
                                          0,
                                          {}));
    }

    const Vec local_input_adjoints = detail::mapped_call_from_bindings(transformed_core, std::move(bindings), node->reps);

    Vec total = zero_vec(context.wrt.size());
    int input_offset = 0;
    for (int i = 0; i < input_count; ++i) {
        const detail::MappedBindingNode &binding = node->mapped_bindings[static_cast<std::size_t>(i)];
        const Vec repeated_local_adjoint = gather(local_input_adjoints,
                                                  repeated_local_adjoint_indices(node->reps,
                                                                                 binding.local_size,
                                                                                 node->function->info.input_size,
                                                                                 input_offset));
        const Vec source_adjoint = scatter_add(repeated_local_adjoint, detail::materialize_indices(binding, node->reps), binding.source->size);
        total = add(total, reverse_vec_node(binding.source, source_adjoint, context));
        input_offset += binding.local_size;
    }
    return total;
}

Vec reverse_vec_node(const std::shared_ptr<const detail::VecNode> &node, const Vec &adjoint, const ReverseContext &context) {
    if (!node) {
        throw std::runtime_error("invalid vector graph node while building reverse derivative");
    }
    if (!adjoint.valid()) {
        throw std::runtime_error("vector reverse adjoint must be a valid vector expression");
    }
    if (adjoint.size() != node->size) {
        throw std::runtime_error("vector reverse adjoint size must match vector output size");
    }

    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
            return reverse_vector_variable(*node, adjoint, context);
        case GraphNodeKind::VectorParameter:
        case GraphNodeKind::VectorConstant:
            return zero_vec(context.wrt.size());
        case GraphNodeKind::VectorFromElements: {
            Vec total = zero_vec(context.wrt.size());
            for (int i = 0; i < node->size; ++i) {
                total = add(total, reverse_scalar_node(node->elements[static_cast<std::size_t>(i)], adjoint[i], context));
            }
            return total;
        }
        case GraphNodeKind::VectorUnary: {
            const Vec arg = as_vec(node->lhs);
            switch (node->unary) {
                case detail::VecUnaryOp::Neg:
                    return reverse_vec_node(node->lhs, -adjoint, context);
                case detail::VecUnaryOp::Sin:
                    return reverse_vec_node(node->lhs, mul(adjoint, cos(arg)), context);
                case detail::VecUnaryOp::Cos:
                    return reverse_vec_node(node->lhs, scale(constant(-1.0), mul(adjoint, sin(arg))), context);
                case detail::VecUnaryOp::Tan:
                    return reverse_vec_node(node->lhs, mul(adjoint, one_vec(node->size) + tan(arg) * tan(arg)), context);
                case detail::VecUnaryOp::Exp:
                    return reverse_vec_node(node->lhs, mul(adjoint, exp(arg)), context);
                case detail::VecUnaryOp::Log:
                    return reverse_vec_node(node->lhs, div(adjoint, arg), context);
                case detail::VecUnaryOp::Sigmoid: {
                    const Vec sig = sigmoid(arg);
                    return reverse_vec_node(node->lhs, mul(adjoint, mul(sig, one_vec(node->size) - sig)), context);
                }
                case detail::VecUnaryOp::Abs:
                    throw std::runtime_error("reverse derivative of nonsmooth abs is not supported");
                case detail::VecUnaryOp::Sqrt:
                    return reverse_vec_node(node->lhs, div(adjoint, scale(constant(2.0), sqrt(arg))), context);
                case detail::VecUnaryOp::Asin:
                    return reverse_vec_node(node->lhs, div(adjoint, sqrt(one_vec(node->size) - arg * arg)), context);
                case detail::VecUnaryOp::Acos:
                    return reverse_vec_node(node->lhs, scale(constant(-1.0), div(adjoint, sqrt(one_vec(node->size) - arg * arg))), context);
                case detail::VecUnaryOp::Atan:
                    return reverse_vec_node(node->lhs, div(adjoint, one_vec(node->size) + arg * arg), context);
                case detail::VecUnaryOp::Sinh:
                    return reverse_vec_node(node->lhs, mul(adjoint, cosh(arg)), context);
                case detail::VecUnaryOp::Cosh:
                    return reverse_vec_node(node->lhs, mul(adjoint, sinh(arg)), context);
                case detail::VecUnaryOp::Tanh: {
                    const Vec t = tanh(arg);
                    return reverse_vec_node(node->lhs, mul(adjoint, one_vec(node->size) - t * t), context);
                }
                case detail::VecUnaryOp::Log10:
                    return reverse_vec_node(node->lhs, div(adjoint, scale(constant(std::log(10.0)), arg)), context);
            }
            break;
        }
        case GraphNodeKind::VectorBinary: {
            const Vec lhs = as_vec(node->lhs);
            const Vec rhs = as_vec(node->rhs);
            switch (node->binary) {
                case detail::VecBinaryOp::Add:
                    return add(reverse_vec_node(node->lhs, adjoint, context), reverse_vec_node(node->rhs, adjoint, context));
                case detail::VecBinaryOp::Sub:
                    return add(reverse_vec_node(node->lhs, adjoint, context), reverse_vec_node(node->rhs, scale(constant(-1.0), adjoint), context));
                case detail::VecBinaryOp::Mul:
                    return add(reverse_vec_node(node->lhs, mul(adjoint, rhs), context), reverse_vec_node(node->rhs, mul(adjoint, lhs), context));
                case detail::VecBinaryOp::Div:
                    return add(reverse_vec_node(node->lhs, div(adjoint, rhs), context),
                               reverse_vec_node(node->rhs, div(scale(constant(-1.0), mul(adjoint, lhs)), rhs * rhs), context));
                case detail::VecBinaryOp::Pow: {
                    const Vec value = pow(lhs, rhs);
                    return add(reverse_vec_node(node->lhs, mul(adjoint, div(mul(value, rhs), lhs)), context),
                               reverse_vec_node(node->rhs, mul(adjoint, mul(value, log(lhs))), context));
                }
                case detail::VecBinaryOp::Min:
                    throw std::runtime_error("reverse derivative of nonsmooth min is not supported");
                case detail::VecBinaryOp::Max:
                    throw std::runtime_error("reverse derivative of nonsmooth max is not supported");
            }
            break;
        }
        case GraphNodeKind::VectorScale: {
            const Expr scale_expr = as_expr(node->scale);
            const Vec vec_arg = as_vec(node->lhs);
            return add(reverse_vec_node(node->lhs, scale(scale_expr, adjoint), context),
                       reverse_scalar_node(node->scale, dot(adjoint, vec_arg), context));
        }
        case GraphNodeKind::DenseMatVec:
            return reverse_vec_node(node->lhs, transpose(node->dense) * adjoint, context);
        case GraphNodeKind::SparseMatVec:
            return reverse_vec_node(node->lhs, transpose(node->sparse) * adjoint, context);
        case GraphNodeKind::Slice:
            return reverse_vec_node(node->lhs, scatter_slice(adjoint, node->start, node->lhs->size), context);
        case GraphNodeKind::ScatterSlice:
            return reverse_vec_node(node->lhs, adjoint.slice(node->start, node->lhs->size), context);
        case GraphNodeKind::Gather:
            return reverse_vec_node(node->lhs, scatter_add(adjoint, node->indices, node->lhs->size), context);
        case GraphNodeKind::ScatterAdd:
            return reverse_vec_node(node->lhs, gather(adjoint, node->indices), context);
        case GraphNodeKind::Concat:
            return add(reverse_vec_node(node->lhs, adjoint.slice(0, node->lhs->size), context),
                       reverse_vec_node(node->rhs, adjoint.slice(node->lhs->size, node->rhs->size), context));
        case GraphNodeKind::FunctionCall:
            return reverse_function_call(node, adjoint, context);
        case GraphNodeKind::MappedFunctionCall:
            return reverse_mapped_function_call(node, adjoint, context);
        default:
            throw std::runtime_error("unsupported vector graph node while building reverse derivative");
    }

    throw std::runtime_error("unhandled vector graph node while building reverse derivative");
}

} // namespace

Expr Expr::forward_diff(const Vars &wrt, const Vec &seed) const {
    if (!valid()) {
        throw std::runtime_error("cannot build forward derivative of invalid scalar expression");
    }
    const ForwardContext context = make_context(wrt, seed);
    return forward_scalar_node(detail::scalar_node(*this), context);
}

Vec Vec::forward_diff(const Vars &wrt, const Vec &seed) const {
    if (!valid()) {
        throw std::runtime_error("cannot build forward derivative of invalid vector expression");
    }
    const ForwardContext context = make_context(wrt, seed);
    return forward_vec_node(detail::vec_node(*this), context);
}

Vec Expr::reverse_diff(const Vars &wrt) const {
    if (!valid()) {
        throw std::runtime_error("cannot build reverse derivative of invalid scalar expression");
    }
    const ReverseContext context = make_reverse_context(wrt);
    return reverse_scalar_node(detail::scalar_node(*this), one_scalar(), context);
}

Vec Vec::reverse_diff(const Vars &wrt, const Vec &lambda) const {
    if (!valid()) {
        throw std::runtime_error("cannot build reverse derivative of invalid vector expression");
    }
    if (!lambda.valid()) {
        throw std::runtime_error("reverse seed must be a valid vector expression");
    }
    if (lambda.size() != size()) {
        throw std::runtime_error("reverse seed size must match vector expression size");
    }
    const ReverseContext context = make_reverse_context(wrt);
    return reverse_vec_node(detail::vec_node(*this), lambda, context);
}

} // namespace ad
