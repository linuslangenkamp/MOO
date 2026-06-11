// SPDX-License-Identifier: LGPL-3.0-or-later
#include "traversal.h"

#include "function_core.h"
#include "mapping.h"
#include "node.h"

#include <set>
#include <stdexcept>

namespace ad::detail {
namespace {

void append_var_once(Vars &out, std::set<VarId> &seen, const Var &var) {
    if (var.valid() && seen.insert(var.id()).second) {
        out.append(var);
    }
}

void append_param_once(Params &out, std::set<ParamId> &seen, const Param &param) {
    if (param.valid() && seen.insert(param.id()).second) {
        out.append(param);
    }
}

void collect_vec_vars(const std::shared_ptr<const VecNode> &node, Vars &out, std::set<VarId> &seen);
void collect_vec_params(const std::shared_ptr<const VecNode> &node, Params &out, std::set<ParamId> &seen);

void collect_mapped_source_vars(const MappedBindingNode &binding, Vars &out, std::set<VarId> &seen) {
    if (!binding.source) {
        throw std::runtime_error("mapped function binding is missing source while collecting variables");
    }
    if (binding.source->kind == GraphNodeKind::VectorVariable) {
        for (int rep = 0; rep < binding.reps; ++rep) {
            for (int component = 0; component < binding.local_size; ++component) {
                const int index = mapped_index(binding, binding.reps, rep, component);
                append_var_once(out, seen, binding.source->vars[static_cast<std::size_t>(index)]);
            }
        }
        return;
    }
    if (binding.source->kind == GraphNodeKind::VectorParameter) {
        return;
    }
    collect_vec_vars(binding.source, out, seen);
}

void collect_mapped_source_params(const MappedBindingNode &binding, Params &out, std::set<ParamId> &seen) {
    if (!binding.source) {
        throw std::runtime_error("mapped function binding is missing source while collecting parameters");
    }
    if (binding.source->kind == GraphNodeKind::VectorParameter) {
        for (int rep = 0; rep < binding.reps; ++rep) {
            for (int component = 0; component < binding.local_size; ++component) {
                const int index = mapped_index(binding, binding.reps, rep, component);
                append_param_once(out, seen, binding.source->params[static_cast<std::size_t>(index)]);
            }
        }
        return;
    }
    if (binding.source->kind == GraphNodeKind::VectorVariable) {
        return;
    }
    collect_vec_params(binding.source, out, seen);
}

void collect_scalar_vars(const std::shared_ptr<const ScalarNode> &node, Vars &out, std::set<VarId> &seen) {
    if (!node) {
        throw std::runtime_error("invalid scalar graph node while collecting variables");
    }

    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
        case GraphNodeKind::ScalarParameter:
            break;
        case GraphNodeKind::ScalarVariable:
            append_var_once(out, seen, node->var);
            break;
        case GraphNodeKind::ScalarUnary:
            collect_scalar_vars(node->lhs, out, seen);
            break;
        case GraphNodeKind::ScalarBinary:
            collect_scalar_vars(node->lhs, out, seen);
            collect_scalar_vars(node->rhs, out, seen);
            break;
        case GraphNodeKind::VectorElement:
        case GraphNodeKind::Sum:
            collect_vec_vars(node->vec, out, seen);
            break;
        case GraphNodeKind::Dot:
            collect_vec_vars(node->vec_lhs, out, seen);
            collect_vec_vars(node->vec_rhs, out, seen);
            break;
        default:
            throw std::runtime_error("unsupported scalar graph node while collecting variables");
    }
}

void collect_scalar_params(const std::shared_ptr<const ScalarNode> &node, Params &out, std::set<ParamId> &seen) {
    if (!node) {
        throw std::runtime_error("invalid scalar graph node while collecting parameters");
    }

    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
        case GraphNodeKind::ScalarVariable:
            break;
        case GraphNodeKind::ScalarParameter:
            append_param_once(out, seen, node->param);
            break;
        case GraphNodeKind::ScalarUnary:
            collect_scalar_params(node->lhs, out, seen);
            break;
        case GraphNodeKind::ScalarBinary:
            collect_scalar_params(node->lhs, out, seen);
            collect_scalar_params(node->rhs, out, seen);
            break;
        case GraphNodeKind::VectorElement:
        case GraphNodeKind::Sum:
            collect_vec_params(node->vec, out, seen);
            break;
        case GraphNodeKind::Dot:
            collect_vec_params(node->vec_lhs, out, seen);
            collect_vec_params(node->vec_rhs, out, seen);
            break;
        default:
            throw std::runtime_error("unsupported scalar graph node while collecting parameters");
    }
}

void collect_vec_vars(const std::shared_ptr<const VecNode> &node, Vars &out, std::set<VarId> &seen) {
    if (!node) {
        throw std::runtime_error("invalid vector graph node while collecting variables");
    }

    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
            for (const Var &var : node->vars) {
                append_var_once(out, seen, var);
            }
            break;
        case GraphNodeKind::VectorParameter:
        case GraphNodeKind::VectorConstant:
            break;
        case GraphNodeKind::VectorFromElements:
            for (const auto &element : node->elements) {
                collect_scalar_vars(element, out, seen);
            }
            break;
        case GraphNodeKind::VectorUnary:
        case GraphNodeKind::DenseMatVec:
        case GraphNodeKind::SparseMatVec:
        case GraphNodeKind::Slice:
        case GraphNodeKind::ScatterSlice:
        case GraphNodeKind::Gather:
        case GraphNodeKind::ScatterAdd:
            collect_vec_vars(node->lhs, out, seen);
            break;
        case GraphNodeKind::VectorBinary:
        case GraphNodeKind::Concat:
        case GraphNodeKind::SymbolicMatVec:
        case GraphNodeKind::SymbolicMatMul:
        case GraphNodeKind::OuterProduct:
        case GraphNodeKind::LinearSolve:
            collect_vec_vars(node->lhs, out, seen);
            collect_vec_vars(node->rhs, out, seen);
            break;
        case GraphNodeKind::VectorScale:
            collect_scalar_vars(node->scale, out, seen);
            collect_vec_vars(node->lhs, out, seen);
            break;
        case GraphNodeKind::FunctionCall:
            for (const auto &argument : node->arguments) {
                collect_vec_vars(argument, out, seen);
            }
            break;
        case GraphNodeKind::MappedFunctionCall:
            for (const auto &binding : node->mapped_bindings) {
                collect_mapped_source_vars(binding, out, seen);
            }
            break;
        case GraphNodeKind::MapAccumCall:
            collect_vec_vars(node->lhs, out, seen);
            for (const auto &binding : node->mapped_bindings) {
                collect_mapped_source_vars(binding, out, seen);
            }
            break;
        default:
            throw std::runtime_error("unsupported vector graph node while collecting variables");
    }
}

void collect_vec_params(const std::shared_ptr<const VecNode> &node, Params &out, std::set<ParamId> &seen) {
    if (!node) {
        throw std::runtime_error("invalid vector graph node while collecting parameters");
    }

    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
        case GraphNodeKind::VectorConstant:
            break;
        case GraphNodeKind::VectorParameter:
            for (const Param &param : node->params) {
                append_param_once(out, seen, param);
            }
            break;
        case GraphNodeKind::VectorFromElements:
            for (const auto &element : node->elements) {
                collect_scalar_params(element, out, seen);
            }
            break;
        case GraphNodeKind::VectorUnary:
        case GraphNodeKind::DenseMatVec:
        case GraphNodeKind::SparseMatVec:
        case GraphNodeKind::Slice:
        case GraphNodeKind::ScatterSlice:
        case GraphNodeKind::Gather:
        case GraphNodeKind::ScatterAdd:
            collect_vec_params(node->lhs, out, seen);
            break;
        case GraphNodeKind::VectorBinary:
        case GraphNodeKind::Concat:
        case GraphNodeKind::SymbolicMatVec:
        case GraphNodeKind::SymbolicMatMul:
        case GraphNodeKind::OuterProduct:
        case GraphNodeKind::LinearSolve:
            collect_vec_params(node->lhs, out, seen);
            collect_vec_params(node->rhs, out, seen);
            break;
        case GraphNodeKind::VectorScale:
            collect_scalar_params(node->scale, out, seen);
            collect_vec_params(node->lhs, out, seen);
            break;
        case GraphNodeKind::FunctionCall:
            if (!node->function) {
                throw std::runtime_error("function call node is missing callee while collecting parameters");
            }
            for (const Param &param : node->function->parameters.values()) {
                append_param_once(out, seen, param);
            }
            for (const auto &argument : node->arguments) {
                collect_vec_params(argument, out, seen);
            }
            break;
        case GraphNodeKind::MappedFunctionCall:
            if (!node->function) {
                throw std::runtime_error("mapped function call node is missing callee while collecting parameters");
            }
            for (const Param &param : node->function->parameters.values()) {
                append_param_once(out, seen, param);
            }
            for (const auto &binding : node->mapped_bindings) {
                collect_mapped_source_params(binding, out, seen);
            }
            break;
        case GraphNodeKind::MapAccumCall:
            if (!node->function) {
                throw std::runtime_error("map-accum call node is missing callee while collecting parameters");
            }
            collect_vec_params(node->lhs, out, seen);
            for (const Param &param : node->function->parameters.values()) {
                append_param_once(out, seen, param);
            }
            for (const auto &binding : node->mapped_bindings) {
                collect_mapped_source_params(binding, out, seen);
            }
            break;
        default:
            throw std::runtime_error("unsupported vector graph node while collecting parameters");
    }
}

} // namespace

Vars collect_variables(const Expr &expr) {
    if (!expr.valid()) {
        throw std::runtime_error("cannot collect variables from invalid scalar expression");
    }
    Vars out;
    std::set<VarId> seen;
    collect_scalar_vars(scalar_node(expr), out, seen);
    return out;
}

Vars collect_variables(const Vec &vec) {
    if (!vec.valid()) {
        throw std::runtime_error("cannot collect variables from invalid vector expression");
    }
    Vars out;
    std::set<VarId> seen;
    collect_vec_vars(vec_node(vec), out, seen);
    return out;
}

Params collect_parameters(const Expr &expr) {
    if (!expr.valid()) {
        throw std::runtime_error("cannot collect parameters from invalid scalar expression");
    }
    Params out;
    std::set<ParamId> seen;
    collect_scalar_params(scalar_node(expr), out, seen);
    return out;
}

Params collect_parameters(const Vec &vec) {
    if (!vec.valid()) {
        throw std::runtime_error("cannot collect parameters from invalid vector expression");
    }
    Params out;
    std::set<ParamId> seen;
    collect_vec_params(vec_node(vec), out, seen);
    return out;
}

bool is_vector_variable_group(const Vec &vec) {
    return vec.valid() && vec_node(vec)->kind == GraphNodeKind::VectorVariable;
}

bool is_vector_parameter_group(const Vec &vec) {
    return vec.valid() && vec_node(vec)->kind == GraphNodeKind::VectorParameter;
}

VecVariableGroupInfo vector_variable_group_info(const Vec &vec) {
    if (!is_vector_variable_group(vec)) {
        throw std::runtime_error("vector expression is not a variable group");
    }
    const auto &node = vec_node(vec);
    VecVariableGroupInfo info;
    info.label = node->label;
    info.node_id = node->id;
    info.size = node->size;
    info.vars = node->vars;
    info.group_id = info.vars.empty() ? 0 : info.vars.front().group_id();
    return info;
}

VecParameterGroupInfo vector_parameter_group_info(const Vec &vec) {
    if (!is_vector_parameter_group(vec)) {
        throw std::runtime_error("vector expression is not a parameter group");
    }
    const auto &node = vec_node(vec);
    VecParameterGroupInfo info;
    info.label = node->label;
    info.node_id = node->id;
    info.size = node->size;
    info.params = node->params;
    info.group_id = info.params.empty() ? 0 : info.params.front().group_id();
    return info;
}

} // namespace ad::detail
