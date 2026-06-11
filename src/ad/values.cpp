// SPDX-License-Identifier: LGPL-3.0-or-later
#include "values.h"

#include "expr.h"
#include "function.h"
#include "vec.h"
#include "detail/traversal.h"

#include <stdexcept>

namespace ad {

void Values::clear() {
    vars_.clear();
    params_.clear();
}

void Values::set(const Var &var, double value) {
    if (!var.valid()) {
        throw std::runtime_error("cannot set value for invalid variable");
    }
    vars_[var.id()] = value;
}

void Values::set(const Param &param, double value) {
    if (!param.valid()) {
        throw std::runtime_error("cannot set value for invalid parameter");
    }
    params_[param.id()] = value;
}

void Values::set(const Expr &symbol, double value) {
    if (symbol.is_variable()) {
        set(symbol.var(), value);
        return;
    }
    if (symbol.is_parameter()) {
        set(symbol.param(), value);
        return;
    }
    throw std::runtime_error("Values::set(Expr, double) requires a scalar variable or parameter");
}

void Values::set(const Vec &symbol_group, const double *values, int size) {
    if (!symbol_group.valid()) {
        throw std::runtime_error("cannot set values for invalid vector expression");
    }
    if (values == nullptr && size > 0) {
        throw std::runtime_error("cannot set vector values from null pointer");
    }
    if (detail::is_vector_variable_group(symbol_group)) {
        detail::VecVariableGroupInfo group = detail::vector_variable_group_info(symbol_group);
        if (size != group.size) {
            throw std::runtime_error("variable vector value size does not match group size");
        }
        for (int i = 0; i < size; ++i) {
            set(group.vars[static_cast<std::size_t>(i)], values[i]);
        }
        return;
    }
    if (detail::is_vector_parameter_group(symbol_group)) {
        detail::VecParameterGroupInfo group = detail::vector_parameter_group_info(symbol_group);
        if (size != group.size) {
            throw std::runtime_error("parameter vector value size does not match group size");
        }
        for (int i = 0; i < size; ++i) {
            set(group.params[static_cast<std::size_t>(i)], values[i]);
        }
        return;
    }
    throw std::runtime_error("Values::set(Vec, ...) requires a vector variable or parameter group");
}

void Values::set(const Vec &symbol_group, const std::vector<double> &values) {
    set(symbol_group, values.data(), static_cast<int>(values.size()));
}

bool Values::has(const Var &var) const {
    return var.valid() && vars_.find(var.id()) != vars_.end();
}

bool Values::has(const Param &param) const {
    return param.valid() && params_.find(param.id()) != params_.end();
}

double Values::get(const Var &var) const {
    if (!var.valid()) {
        throw std::runtime_error("cannot read invalid variable value");
    }
    const auto it = vars_.find(var.id());
    if (it == vars_.end()) {
        throw std::runtime_error("missing value for variable: " + var.label());
    }
    return it->second;
}

double Values::get(const Param &param) const {
    if (!param.valid()) {
        throw std::runtime_error("cannot read invalid parameter value");
    }
    const auto it = params_.find(param.id());
    if (it == params_.end()) {
        throw std::runtime_error("missing value for parameter: " + param.label());
    }
    return it->second;
}

void EvalWorkspace::clear() {
    used_buffers_ = 0;
    used_pointer_buffers_ = 0;
    var_frames_.clear();
    param_frames_.clear();
}

void EvalWorkspace::reserve(const Function &function) {
    if (!function.valid()) {
        throw std::runtime_error("cannot reserve evaluation workspace for invalid function");
    }
    const std::size_t baseline = static_cast<std::size_t>(function.input_size() + function.output_size() + function.parameters().size() + 1024);
    var_frames_.reserve(32);
    param_frames_.reserve(8);
    buffers_.reserve(128);
    pointer_buffers_.reserve(16);
    reserve_doubles(baseline);
}

void EvalWorkspace::reserve_doubles(std::size_t count) {
    if (count == 0) {
        return;
    }
    if (buffers_.empty()) {
        Buffer buffer;
        buffer.data = std::make_unique<double[]>(count);
        buffer.capacity = static_cast<int>(count);
        buffers_.push_back(std::move(buffer));
        return;
    }
    if (buffers_.front().capacity < static_cast<int>(count)) {
        buffers_.front().data = std::make_unique<double[]>(count);
        buffers_.front().capacity = static_cast<int>(count);
    }
}

} // namespace ad
