// SPDX-License-Identifier: LGPL-3.0-or-later
#include "function.h"

#include "detail/function_core.h"
#include "detail/node.h"
#include "detail/simplify.h"
#include "detail/traversal.h"

#include <mutex>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ad {
namespace {

const detail::FunctionCore &require_core(const std::shared_ptr<const detail::FunctionCore> &core) {
    if (!core) {
        throw std::runtime_error("invalid AD function");
    }
    return *core;
}

void append_input_var(Vars &out, std::set<VarId> &seen, const Var &var) {
    if (!seen.insert(var.id()).second) {
        throw std::runtime_error("duplicate variable ID in function inputs");
    }
    out.append(var);
}

void append_parameter(Params &out, std::set<ParamId> &seen, const Param &param) {
    if (!param.valid()) {
        throw std::runtime_error("function parameter layout contains an invalid parameter");
    }
    if (!seen.insert(param.id()).second) {
        throw std::runtime_error("duplicate parameter ID in function parameter layout");
    }
    out.append(param);
}

bool contains_var_id(const std::set<VarId> &ids, VarId id) {
    return ids.find(id) != ids.end();
}

bool contains_param_id(const std::set<ParamId> &ids, ParamId id) {
    return ids.find(id) != ids.end();
}

std::shared_ptr<const detail::FunctionCore> make_function_core(std::vector<Vec> inputs,
                                                               Params parameters,
                                                               Vec output,
                                                               bool explicit_parameters) {
    if (!output.valid()) {
        throw std::runtime_error("function output must be a valid vector expression");
    }

    auto core = std::make_shared<detail::FunctionCore>();
    core->inputs = std::move(inputs);
    core->output = detail::simplify_vec(output);

    std::set<NodeId> input_nodes;
    std::set<VarId> input_var_ids;
    std::vector<FunctionInputInfo> input_infos;

    for (const Vec &input : core->inputs) {
        if (!input.valid()) {
            throw std::runtime_error("function input must be a valid vector variable group");
        }
        if (!detail::is_vector_variable_group(input)) {
            throw std::runtime_error("function input must be a vector variable group");
        }

        detail::VecVariableGroupInfo group = detail::vector_variable_group_info(input);
        if (group.size <= 0) {
            throw std::runtime_error("function input variable group must be non-empty");
        }
        if (!input_nodes.insert(group.node_id).second) {
            throw std::runtime_error("duplicate function input variable group");
        }

        for (const Var &var : group.vars) {
            append_input_var(core->input_vars, input_var_ids, var);
        }

        FunctionInputInfo input_info;
        input_info.label = group.label;
        input_info.size = group.size;
        input_info.node_id = group.node_id;
        input_infos.push_back(std::move(input_info));
    }

    Vars output_vars = detail::collect_variables(core->output);
    for (const Var &var : output_vars.values()) {
        if (!contains_var_id(input_var_ids, var.id())) {
            throw std::runtime_error("function output depends on undeclared variable: " + var.label());
        }
    }

    Params output_parameters = detail::collect_parameters(core->output);
    if (explicit_parameters) {
        std::set<ParamId> declared_parameter_ids;
        for (const Param &param : parameters.values()) {
            append_parameter(core->parameters, declared_parameter_ids, param);
        }
        for (const Param &param : output_parameters.values()) {
            if (!contains_param_id(declared_parameter_ids, param.id())) {
                throw std::runtime_error("function output depends on undeclared parameter: " + param.label());
            }
        }
    } else {
        core->parameters = std::move(output_parameters);
    }

    core->info.input_count = static_cast<int>(core->inputs.size());
    core->info.input_size = core->input_vars.size();
    core->info.output_size = core->output.size();
    core->info.parameter_count = core->parameters.size();
    core->info.inputs = std::move(input_infos);
    core->info.output_graph = inspect(core->output);

    return core;
}

std::shared_ptr<const detail::FunctionCore> make_function_core(std::vector<Vec> inputs, Vec output) {
    return make_function_core(std::move(inputs), Params{}, std::move(output), false);
}

std::shared_ptr<const detail::FunctionCore> make_function_core(std::vector<Vec> inputs, Params parameters, Vec output) {
    return make_function_core(std::move(inputs), std::move(parameters), std::move(output), true);
}

Vec zero_vec(int size) {
    return vec_constant(std::vector<double>(static_cast<std::size_t>(size), 0.0));
}

Vec concat_all(const std::vector<Vec> &parts) {
    if (parts.empty()) {
        return zero_vec(0);
    }

    Vec out = parts.front();
    for (std::size_t i = 1; i < parts.size(); ++i) {
        out = concat(out, parts[i]);
    }
    return out;
}

std::vector<Vec> make_input_tangent_groups(const detail::FunctionCore &core) {
    std::vector<Vec> groups;
    groups.reserve(core.inputs.size());
    for (std::size_t i = 0; i < core.inputs.size(); ++i) {
        const int size = core.info.inputs[i].size;
        groups.push_back(vec_variable("__d" + std::to_string(i), size));
    }
    return groups;
}

std::vector<Vec> forward_call_arguments(const detail::FunctionCore &core, const Vec &seed) {
    std::vector<Vec> arguments = core.inputs;
    arguments.reserve(core.inputs.size() * 2);

    int offset = 0;
    for (const FunctionInputInfo &input : core.info.inputs) {
        arguments.push_back(seed.slice(offset, input.size));
        offset += input.size;
    }
    return arguments;
}

std::vector<Vec> reverse_call_arguments(const detail::FunctionCore &core, const Vec &lambda) {
    std::vector<Vec> arguments = core.inputs;
    arguments.reserve(core.inputs.size() + (core.info.output_size > 0 ? 1 : 0));
    if (core.info.output_size > 0) {
        arguments.push_back(lambda);
    }
    return arguments;
}

} // namespace

Function::Function(std::vector<Vec> inputs, Vec output)
    : core_(make_function_core(std::move(inputs), std::move(output))) {}

Function::Function(std::vector<Vec> inputs, Params parameters, Vec output)
    : core_(make_function_core(std::move(inputs), std::move(parameters), std::move(output))) {}

Function::Function(std::shared_ptr<const detail::FunctionCore> core)
    : core_(std::move(core)) {}

bool Function::valid() const {
    return static_cast<bool>(core_);
}

int Function::input_count() const {
    return core_ ? core_->info.input_count : 0;
}

int Function::input_size() const {
    return core_ ? core_->info.input_size : 0;
}

int Function::output_size() const {
    return core_ ? core_->info.output_size : 0;
}

const std::vector<Vec> &Function::inputs() const {
    return require_core(core_).inputs;
}

const Vec &Function::output() const {
    return require_core(core_).output;
}

Vars Function::input_vars() const {
    return require_core(core_).input_vars;
}

Params Function::parameters() const {
    return require_core(core_).parameters;
}

FunctionInfo Function::info() const {
    return require_core(core_).info;
}

Vec Function::call(std::vector<Vec> arguments) const {
    const detail::FunctionCore &core = require_core(core_);
    if (static_cast<int>(arguments.size()) != core.info.input_count) {
        throw std::runtime_error("function call argument count does not match input count");
    }

    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::FunctionCall;
    node->size = core.info.output_size;
    node->function = core_;
    node->arguments.reserve(arguments.size());

    for (std::size_t i = 0; i < arguments.size(); ++i) {
        Vec argument = detail::simplify_vec(arguments[i]);
        if (!argument.valid()) {
            throw std::runtime_error("function call argument must be a valid vector expression");
        }
        const int expected_size = core.info.inputs[i].size;
        if (argument.size() != expected_size) {
            throw std::runtime_error("function call argument size does not match input size");
        }
        node->arguments.push_back(detail::vec_node(argument));
    }

    return detail::vec_from_node(node);
}

Function Function::forward_function() const {
    const std::shared_ptr<const detail::FunctionCore> core = core_;
    require_core(core);

    {
        std::lock_guard<std::mutex> lock(core->transform_mutex);
        if (core->forward_cache) {
            return Function(core->forward_cache);
        }
    }

    std::vector<Vec> tangent_groups = make_input_tangent_groups(*core);
    Vec seed = concat_all(tangent_groups);
    Vec output = core->output.forward_diff(core->input_vars, seed);

    std::vector<Vec> inputs = core->inputs;
    inputs.insert(inputs.end(), tangent_groups.begin(), tangent_groups.end());
    Function transformed(std::move(inputs), core->parameters, output);

    {
        std::lock_guard<std::mutex> lock(core->transform_mutex);
        if (!core->forward_cache) {
            core->forward_cache = transformed.core_;
        }
        return Function(core->forward_cache);
    }
}

Function Function::reverse_function() const {
    const std::shared_ptr<const detail::FunctionCore> core = core_;
    require_core(core);

    {
        std::lock_guard<std::mutex> lock(core->transform_mutex);
        if (core->reverse_cache) {
            return Function(core->reverse_cache);
        }
    }

    std::vector<Vec> inputs = core->inputs;
    Vec lambda = zero_vec(0);
    if (core->info.output_size > 0) {
        lambda = vec_variable("__lambda", core->info.output_size);
        inputs.push_back(lambda);
    }

    Vec output = core->output.reverse_diff(core->input_vars, lambda);
    Function transformed(std::move(inputs), core->parameters, output);

    {
        std::lock_guard<std::mutex> lock(core->transform_mutex);
        if (!core->reverse_cache) {
            core->reverse_cache = transformed.core_;
        }
        return Function(core->reverse_cache);
    }
}

Vec Function::forward(const Vec &seed) const {
    const detail::FunctionCore &core = require_core(core_);
    if (!seed.valid()) {
        throw std::runtime_error("function forward seed must be a valid vector expression");
    }
    if (seed.size() != core.info.input_size) {
        throw std::runtime_error("function forward seed size must match flattened input size");
    }
    return forward_function().call(forward_call_arguments(core, seed));
}

Vec Function::reverse(const Vec &lambda) const {
    const detail::FunctionCore &core = require_core(core_);
    if (!lambda.valid()) {
        throw std::runtime_error("function reverse seed must be a valid vector expression");
    }
    if (lambda.size() != core.info.output_size) {
        throw std::runtime_error("function reverse seed size must match output size");
    }
    return reverse_function().call(reverse_call_arguments(core, lambda));
}

SparsityPattern Function::jacobian_sparsity() const {
    const detail::FunctionCore &core = require_core(core_);
    return sparsity(core.output, core.input_vars);
}

namespace detail {

const std::shared_ptr<const FunctionCore> &function_core(const Function &function) {
    return function.core_;
}

Function function_from_core(std::shared_ptr<const FunctionCore> core) {
    return Function(std::move(core));
}

} // namespace detail

} // namespace ad
