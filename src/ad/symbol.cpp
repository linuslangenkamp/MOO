// SPDX-License-Identifier: LGPL-3.0-or-later
#include "symbol.h"

#include <atomic>
#include <stdexcept>
#include <utility>

namespace ad {
namespace {

VarId next_var_id() {
    static std::atomic<VarId> next{1};
    return next.fetch_add(1);
}

ParamId next_param_id() {
    static std::atomic<ParamId> next{1};
    return next.fetch_add(1);
}

} // namespace

Var::Var(VarId id, std::string label, int component, SymbolGroupId group_id)
    : id_(id),
      label_(std::move(label)),
      component_(component),
      group_id_(group_id) {}

VarId Var::id() const {
    return id_;
}

SymbolGroupId Var::group_id() const {
    return group_id_;
}

const std::string &Var::label() const {
    return label_;
}

int Var::component() const {
    return component_;
}

bool Var::valid() const {
    return id_ != 0;
}

bool Var::operator==(const Var &other) const {
    return id_ == other.id_;
}

bool Var::operator!=(const Var &other) const {
    return !(*this == other);
}

Param::Param(ParamId id, std::string label, int component, SymbolGroupId group_id)
    : id_(id),
      label_(std::move(label)),
      component_(component),
      group_id_(group_id) {}

ParamId Param::id() const {
    return id_;
}

SymbolGroupId Param::group_id() const {
    return group_id_;
}

const std::string &Param::label() const {
    return label_;
}

int Param::component() const {
    return component_;
}

bool Param::valid() const {
    return id_ != 0;
}

bool Param::operator==(const Param &other) const {
    return id_ == other.id_;
}

bool Param::operator!=(const Param &other) const {
    return !(*this == other);
}

Vars::Vars(std::vector<Var> values)
    : values_(std::move(values)) {}

int Vars::size() const {
    return static_cast<int>(values_.size());
}

bool Vars::empty() const {
    return values_.empty();
}

const Var &Vars::operator[](int index) const {
    if (index < 0 || index >= size()) {
        throw std::out_of_range("variable index out of range");
    }
    return values_[static_cast<std::size_t>(index)];
}

const std::vector<Var> &Vars::values() const {
    return values_;
}

void Vars::append(const Var &var) {
    if (!var.valid()) {
        throw std::runtime_error("invalid variable");
    }
    values_.push_back(var);
}

Params::Params(std::vector<Param> values)
    : values_(std::move(values)) {}

int Params::size() const {
    return static_cast<int>(values_.size());
}

bool Params::empty() const {
    return values_.empty();
}

const Param &Params::operator[](int index) const {
    if (index < 0 || index >= size()) {
        throw std::out_of_range("parameter index out of range");
    }
    return values_[static_cast<std::size_t>(index)];
}

const std::vector<Param> &Params::values() const {
    return values_;
}

void Params::append(const Param &param) {
    if (!param.valid()) {
        throw std::runtime_error("invalid parameter");
    }
    values_.push_back(param);
}

namespace detail {

SymbolGroupId next_symbol_group_id() {
    static std::atomic<SymbolGroupId> next{1};
    return next.fetch_add(1);
}

Var make_var(std::string label, int component, SymbolGroupId group_id) {
    return Var(next_var_id(), std::move(label), component, group_id);
}

Param make_param(std::string label, int component, SymbolGroupId group_id) {
    return Param(next_param_id(), std::move(label), component, group_id);
}

} // namespace detail
} // namespace ad
