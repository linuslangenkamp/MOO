// SPDX-License-Identifier: LGPL-3.0-or-later
#include "values.h"

#include "expr.h"
#include "function.h"
#include "vec.h"
#include "detail/function_core.h"
#include "detail/mapping.h"
#include "detail/node.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

namespace ad {
namespace detail {

struct EvalWorkspaceAccess {
    struct Mark {
        std::size_t buffers = 0;
        std::size_t pointer_buffers = 0;
    };

    static Mark mark(EvalWorkspace &workspace) {
        return Mark{workspace.used_buffers_, workspace.used_pointer_buffers_};
    }

    static double *allocate(EvalWorkspace &workspace, int size) {
        if (size < 0) {
            throw std::runtime_error("cannot allocate negative evaluation workspace size");
        }
        if (size == 0) {
            return nullptr;
        }
        const std::size_t index = workspace.used_buffers_++;
        if (index == workspace.buffers_.size()) {
            EvalWorkspace::Buffer buffer;
            buffer.data = std::make_unique<double[]>(static_cast<std::size_t>(size));
            buffer.capacity = size;
            workspace.buffers_.push_back(std::move(buffer));
        } else if (workspace.buffers_[index].capacity < size) {
            workspace.buffers_[index].data = std::make_unique<double[]>(static_cast<std::size_t>(size));
            workspace.buffers_[index].capacity = size;
        }
        return workspace.buffers_[index].data.get();
    }

    static double **allocate_pointer_array(EvalWorkspace &workspace, int size) {
        if (size < 0) {
            throw std::runtime_error("cannot allocate negative evaluation pointer workspace size");
        }
        if (size == 0) {
            return nullptr;
        }
        const std::size_t index = workspace.used_pointer_buffers_++;
        if (index == workspace.pointer_buffers_.size()) {
            EvalWorkspace::PointerBuffer buffer;
            buffer.data = std::make_unique<double *[]>(static_cast<std::size_t>(size));
            buffer.capacity = size;
            workspace.pointer_buffers_.push_back(std::move(buffer));
        } else if (workspace.pointer_buffers_[index].capacity < size) {
            workspace.pointer_buffers_[index].data = std::make_unique<double *[]>(static_cast<std::size_t>(size));
            workspace.pointer_buffers_[index].capacity = size;
        }
        return workspace.pointer_buffers_[index].data.get();
    }

    static void release(EvalWorkspace &workspace, Mark mark) {
        workspace.used_buffers_ = mark.buffers;
        workspace.used_pointer_buffers_ = mark.pointer_buffers;
    }

    static void push_var_frame(EvalWorkspace &workspace, const Vars &vars, const double *values, int size) {
        if (size != vars.size()) {
            throw std::runtime_error("evaluation variable frame size does not match variable layout");
        }
        if (values == nullptr && size > 0) {
            throw std::runtime_error("evaluation variable frame has null values");
        }
        workspace.var_frames_.push_back(EvalWorkspace::VarFrame{&vars, values, size});
    }

    static void pop_var_frame(EvalWorkspace &workspace) {
        if (workspace.var_frames_.empty()) {
            throw std::runtime_error("evaluation variable frame stack underflow");
        }
        workspace.var_frames_.pop_back();
    }

    static void push_param_frame(EvalWorkspace &workspace, const Params &params, const double *values, int size) {
        if (size != params.size()) {
            throw std::runtime_error("evaluation parameter frame size does not match parameter layout");
        }
        if (values == nullptr && size > 0) {
            throw std::runtime_error("evaluation parameter frame has null values");
        }
        workspace.param_frames_.push_back(EvalWorkspace::ParamFrame{&params, values, size});
    }

    static void pop_param_frame(EvalWorkspace &workspace) {
        if (workspace.param_frames_.empty()) {
            throw std::runtime_error("evaluation parameter frame stack underflow");
        }
        workspace.param_frames_.pop_back();
    }

    static const std::vector<EvalWorkspace::VarFrame> &var_frames(const EvalWorkspace &workspace) {
        return workspace.var_frames_;
    }

    static const std::vector<EvalWorkspace::ParamFrame> &param_frames(const EvalWorkspace &workspace) {
        return workspace.param_frames_;
    }
};

} // namespace detail
namespace {

struct EvalContext {
    EvalWorkspace &workspace;
    const Values *values = nullptr;
};

void require_output(double *out, int size) {
    if (out == nullptr && size > 0) {
        throw std::runtime_error("evaluation output buffer is null");
    }
}

double lookup_var(const Var &var, const EvalContext &context) {
    const auto &frames = detail::EvalWorkspaceAccess::var_frames(context.workspace);
    for (auto frame_it = frames.rbegin(); frame_it != frames.rend(); ++frame_it) {
        for (int i = 0; i < frame_it->size; ++i) {
            if ((*frame_it->vars)[i].id() == var.id()) {
                return frame_it->values[i];
            }
        }
    }
    if (context.values) {
        return context.values->get(var);
    }
    throw std::runtime_error("missing value for variable: " + var.label());
}

double lookup_param(const Param &param, const EvalContext &context) {
    const auto &frames = detail::EvalWorkspaceAccess::param_frames(context.workspace);
    for (auto frame_it = frames.rbegin(); frame_it != frames.rend(); ++frame_it) {
        for (int i = 0; i < frame_it->size; ++i) {
            if ((*frame_it->params)[i].id() == param.id()) {
                return frame_it->values[i];
            }
        }
    }
    if (context.values) {
        return context.values->get(param);
    }
    throw std::runtime_error("missing value for parameter: " + param.label());
}

void eval_vec_node(const std::shared_ptr<const detail::VecNode> &node, EvalContext &context, double *out, int out_size);

double eval_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node, EvalContext &context) {
    if (!node) {
        throw std::runtime_error("cannot evaluate invalid scalar node");
    }

    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
            return node->value;
        case GraphNodeKind::ScalarVariable:
            return lookup_var(node->var, context);
        case GraphNodeKind::ScalarParameter:
            return lookup_param(node->param, context);
        case GraphNodeKind::ScalarUnary: {
            const double arg = eval_scalar_node(node->lhs, context);
            switch (node->unary) {
                case detail::ScalarUnaryOp::Neg:
                    return -arg;
                case detail::ScalarUnaryOp::Sin:
                    return std::sin(arg);
                case detail::ScalarUnaryOp::Cos:
                    return std::cos(arg);
                case detail::ScalarUnaryOp::Tan:
                    return std::tan(arg);
                case detail::ScalarUnaryOp::Exp:
                    return std::exp(arg);
                case detail::ScalarUnaryOp::Log:
                    return std::log(arg);
                case detail::ScalarUnaryOp::PowConst:
                    return std::pow(arg, node->value);
            }
            break;
        }
        case GraphNodeKind::ScalarBinary: {
            const double lhs = eval_scalar_node(node->lhs, context);
            const double rhs = eval_scalar_node(node->rhs, context);
            switch (node->binary) {
                case detail::ScalarBinaryOp::Add:
                    return lhs + rhs;
                case detail::ScalarBinaryOp::Sub:
                    return lhs - rhs;
                case detail::ScalarBinaryOp::Mul:
                    return lhs * rhs;
                case detail::ScalarBinaryOp::Div:
                    return lhs / rhs;
            }
            break;
        }
        case GraphNodeKind::VectorElement: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *values = detail::EvalWorkspaceAccess::allocate(context.workspace, node->vec->size);
            eval_vec_node(node->vec, context, values, node->vec->size);
            const double out = values[node->index];
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            return out;
        }
        case GraphNodeKind::Sum: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *values = detail::EvalWorkspaceAccess::allocate(context.workspace, node->vec->size);
            eval_vec_node(node->vec, context, values, node->vec->size);
            double sum = 0.0;
            for (int i = 0; i < node->vec->size; ++i) {
                sum += values[i];
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            return sum;
        }
        case GraphNodeKind::Dot: {
            if (node->vec_lhs->size != node->vec_rhs->size) {
                throw std::runtime_error("cannot evaluate dot product with mismatched vector sizes");
            }
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *lhs = detail::EvalWorkspaceAccess::allocate(context.workspace, node->vec_lhs->size);
            double *rhs = detail::EvalWorkspaceAccess::allocate(context.workspace, node->vec_rhs->size);
            eval_vec_node(node->vec_lhs, context, lhs, node->vec_lhs->size);
            eval_vec_node(node->vec_rhs, context, rhs, node->vec_rhs->size);
            double dot = 0.0;
            for (int i = 0; i < node->vec_lhs->size; ++i) {
                dot += lhs[i] * rhs[i];
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            return dot;
        }
        default:
            throw std::runtime_error("unsupported scalar node kind in evaluation");
    }

    throw std::runtime_error("unreachable scalar evaluation branch");
}

double eval_vec_unary(detail::VecUnaryOp op, double value) {
    switch (op) {
        case detail::VecUnaryOp::Sin:
            return std::sin(value);
        case detail::VecUnaryOp::Cos:
            return std::cos(value);
        case detail::VecUnaryOp::Tan:
            return std::tan(value);
        case detail::VecUnaryOp::Exp:
            return std::exp(value);
        case detail::VecUnaryOp::Log:
            return std::log(value);
        case detail::VecUnaryOp::Sigmoid:
            return 1.0 / (1.0 + std::exp(-value));
    }
    throw std::runtime_error("unsupported vector unary op in evaluation");
}

double eval_vec_binary(detail::VecBinaryOp op, double lhs, double rhs) {
    switch (op) {
        case detail::VecBinaryOp::Add:
            return lhs + rhs;
        case detail::VecBinaryOp::Sub:
            return lhs - rhs;
        case detail::VecBinaryOp::Mul:
            return lhs * rhs;
        case detail::VecBinaryOp::Div:
            return lhs / rhs;
    }
    throw std::runtime_error("unsupported vector binary op in evaluation");
}

void eval_function_core(const detail::FunctionCore &core, EvalContext &context, const double *input, int input_size, double *out, int out_size) {
    if (input_size != core.info.input_size) {
        throw std::runtime_error("function evaluation input size does not match function layout");
    }
    if (out_size != core.info.output_size) {
        throw std::runtime_error("function evaluation output size does not match function layout");
    }
    require_output(out, out_size);

    detail::EvalWorkspaceAccess::push_var_frame(context.workspace, core.input_vars, input, input_size);
    eval_vec_node(detail::vec_node(core.output), context, out, out_size);
    detail::EvalWorkspaceAccess::pop_var_frame(context.workspace);
}

void eval_vec_node(const std::shared_ptr<const detail::VecNode> &node, EvalContext &context, double *out, int out_size) {
    if (!node) {
        throw std::runtime_error("cannot evaluate invalid vector node");
    }
    if (out_size != node->size) {
        throw std::runtime_error("vector evaluation output size does not match vector size");
    }
    require_output(out, out_size);

    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
            for (int i = 0; i < node->size; ++i) {
                out[i] = lookup_var(node->vars[static_cast<std::size_t>(i)], context);
            }
            break;
        case GraphNodeKind::VectorParameter:
            for (int i = 0; i < node->size; ++i) {
                out[i] = lookup_param(node->params[static_cast<std::size_t>(i)], context);
            }
            break;
        case GraphNodeKind::VectorConstant:
            std::copy(node->constants.begin(), node->constants.end(), out);
            break;
        case GraphNodeKind::VectorFromElements:
            for (int i = 0; i < node->size; ++i) {
                out[i] = eval_scalar_node(node->elements[static_cast<std::size_t>(i)], context);
            }
            break;
        case GraphNodeKind::VectorUnary: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->size);
            eval_vec_node(node->lhs, context, arg, node->size);
            for (int i = 0; i < node->size; ++i) {
                out[i] = eval_vec_unary(node->unary, arg[i]);
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::VectorBinary: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *lhs = detail::EvalWorkspaceAccess::allocate(context.workspace, node->size);
            double *rhs = detail::EvalWorkspaceAccess::allocate(context.workspace, node->size);
            eval_vec_node(node->lhs, context, lhs, node->size);
            eval_vec_node(node->rhs, context, rhs, node->size);
            for (int i = 0; i < node->size; ++i) {
                out[i] = eval_vec_binary(node->binary, lhs[i], rhs[i]);
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::VectorScale: {
            const double scale = eval_scalar_node(node->scale, context);
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->size);
            eval_vec_node(node->lhs, context, arg, node->size);
            for (int i = 0; i < node->size; ++i) {
                out[i] = scale * arg[i];
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::DenseMatVec: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->dense.cols);
            eval_vec_node(node->lhs, context, arg, node->dense.cols);
            for (int row = 0; row < node->dense.rows; ++row) {
                double value = 0.0;
                for (int col = 0; col < node->dense.cols; ++col) {
                    value += node->dense(row, col) * arg[col];
                }
                out[row] = value;
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::SparseMatVec: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->sparse.cols);
            eval_vec_node(node->lhs, context, arg, node->sparse.cols);
            std::fill(out, out + node->size, 0.0);
            for (int k = 0; k < node->sparse.nnz(); ++k) {
                out[node->sparse.row[static_cast<std::size_t>(k)]] += node->sparse.values[static_cast<std::size_t>(k)] * arg[node->sparse.col[static_cast<std::size_t>(k)]];
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::Slice: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->lhs->size);
            eval_vec_node(node->lhs, context, arg, node->lhs->size);
            std::copy(arg + node->start, arg + node->start + node->size, out);
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::ScatterSlice: {
            std::fill(out, out + node->size, 0.0);
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->lhs->size);
            eval_vec_node(node->lhs, context, arg, node->lhs->size);
            std::copy(arg, arg + node->lhs->size, out + node->start);
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::Gather: {
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->lhs->size);
            eval_vec_node(node->lhs, context, arg, node->lhs->size);
            for (int i = 0; i < node->size; ++i) {
                out[i] = arg[node->indices[static_cast<std::size_t>(i)]];
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::ScatterAdd: {
            std::fill(out, out + node->size, 0.0);
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *arg = detail::EvalWorkspaceAccess::allocate(context.workspace, node->lhs->size);
            eval_vec_node(node->lhs, context, arg, node->lhs->size);
            for (int i = 0; i < node->lhs->size; ++i) {
                out[node->indices[static_cast<std::size_t>(i)]] += arg[i];
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::Concat:
            eval_vec_node(node->lhs, context, out, node->lhs->size);
            eval_vec_node(node->rhs, context, out + node->lhs->size, node->rhs->size);
            break;
        case GraphNodeKind::FunctionCall: {
            if (!node->function) {
                throw std::runtime_error("function call node is missing callee during evaluation");
            }
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double *input = detail::EvalWorkspaceAccess::allocate(context.workspace, node->function->info.input_size);
            int offset = 0;
            for (std::size_t i = 0; i < node->arguments.size(); ++i) {
                const int size = node->function->info.inputs[i].size;
                eval_vec_node(node->arguments[i], context, input + offset, size);
                offset += size;
            }
            eval_function_core(*node->function, context, input, node->function->info.input_size, out, out_size);
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        case GraphNodeKind::MappedFunctionCall: {
            if (!node->function) {
                throw std::runtime_error("mapped function call node is missing callee during evaluation");
            }
            const int local_output_size = node->function->info.output_size;
            const auto mark = detail::EvalWorkspaceAccess::mark(context.workspace);
            double **source_values = detail::EvalWorkspaceAccess::allocate_pointer_array(context.workspace, static_cast<int>(node->mapped_bindings.size()));
            for (std::size_t binding_index = 0; binding_index < node->mapped_bindings.size(); ++binding_index) {
                const detail::MappedBindingNode &binding = node->mapped_bindings[binding_index];
                double *source = detail::EvalWorkspaceAccess::allocate(context.workspace, binding.source->size);
                eval_vec_node(binding.source, context, source, binding.source->size);
                source_values[binding_index] = source;
            }
            double *local_input = detail::EvalWorkspaceAccess::allocate(context.workspace, node->function->info.input_size);
            double *local_output = detail::EvalWorkspaceAccess::allocate(context.workspace, local_output_size);
            if (node->mapped_output.mode() != OutputMode::Concat) {
                for (int i = 0; i < out_size; ++i) {
                    out[i] = 0.0;
                }
            }
            for (int rep = 0; rep < node->reps; ++rep) {
                int offset = 0;
                for (std::size_t binding_index = 0; binding_index < node->mapped_bindings.size(); ++binding_index) {
                    const detail::MappedBindingNode &binding = node->mapped_bindings[binding_index];
                    for (int component = 0; component < binding.local_size; ++component) {
                        const int source_index = detail::mapped_index(binding, node->reps, rep, component);
                        local_input[offset + component] = source_values[binding_index][source_index];
                    }
                    offset += binding.local_size;
                }
                double *target = out + rep * local_output_size;
                if (node->mapped_output.mode() == OutputMode::Concat) {
                    eval_function_core(*node->function, context, local_input, node->function->info.input_size, target, local_output_size);
                } else {
                    eval_function_core(*node->function, context, local_input, node->function->info.input_size, local_output, local_output_size);
                    for (int row = 0; row < local_output_size; ++row) {
                        if (node->mapped_output.mode() == OutputMode::Scatter) {
                            const int output_index = node->mapped_output.indices()[static_cast<std::size_t>(rep * local_output_size + row)];
                            out[output_index] += local_output[row];
                        } else if (node->mapped_output.mode() == OutputMode::Sum) {
                            out[row] += local_output[row];
                        } else if (node->mapped_output.mode() == OutputMode::WeightedSum) {
                            out[row] += node->mapped_output.weights()[static_cast<std::size_t>(rep)] * local_output[row];
                        }
                    }
                }
            }
            detail::EvalWorkspaceAccess::release(context.workspace, mark);
            break;
        }
        default:
            throw std::runtime_error("unsupported vector node kind in evaluation");
    }
}

} // namespace

void Expr::eval(const Values &values, EvalWorkspace &workspace, double *out) const {
    if (!valid()) {
        throw std::runtime_error("cannot evaluate invalid scalar expression");
    }
    require_output(out, 1);
    workspace.clear();
    EvalContext context{workspace, &values};
    *out = eval_scalar_node(detail::scalar_node(*this), context);
}

void Vec::eval(const Values &values, EvalWorkspace &workspace, double *out, int out_size) const {
    if (!valid()) {
        throw std::runtime_error("cannot evaluate invalid vector expression");
    }
    workspace.clear();
    EvalContext context{workspace, &values};
    eval_vec_node(detail::vec_node(*this), context, out, out_size);
}

void Function::eval(const double *input, int input_size,
                    const double *params, int param_size,
                    EvalWorkspace &workspace,
                    double *output, int output_size) const {
    const std::shared_ptr<const detail::FunctionCore> &core = detail::function_core(*this);
    if (!core) {
        throw std::runtime_error("cannot evaluate invalid function");
    }
    if (input == nullptr && input_size > 0) {
        throw std::runtime_error("function evaluation input buffer is null");
    }
    if (params == nullptr && param_size > 0) {
        throw std::runtime_error("function evaluation parameter buffer is null");
    }
    if (param_size != core->parameters.size()) {
        throw std::runtime_error("function evaluation parameter size does not match parameter layout");
    }
    workspace.clear();
    EvalContext context{workspace, nullptr};
    detail::EvalWorkspaceAccess::push_param_frame(workspace, core->parameters, params, param_size);
    eval_function_core(*core, context, input, input_size, output, output_size);
    detail::EvalWorkspaceAccess::pop_param_frame(workspace);
}

void Function::eval(const Values &values,
                    EvalWorkspace &workspace,
                    double *output, int output_size) const {
    const std::shared_ptr<const detail::FunctionCore> &core = detail::function_core(*this);
    if (!core) {
        throw std::runtime_error("cannot evaluate invalid function");
    }
    workspace.clear();
    EvalContext context{workspace, &values};
    eval_vec_node(detail::vec_node(core->output), context, output, output_size);
}

} // namespace ad
