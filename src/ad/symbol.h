// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_SYMBOL_H
#define MOO_AD_SYMBOL_H

#include <cstdint>
#include <string>
#include <vector>

namespace ad {

using VarId = std::uint64_t;
using ParamId = std::uint64_t;
using SymbolGroupId = std::uint64_t;

class Expr;
class Vec;
class Var;
class Param;
class Vars;
class Params;

namespace detail {
Var make_var(std::string label, int component, SymbolGroupId group_id);
Param make_param(std::string label, int component, SymbolGroupId group_id);
SymbolGroupId next_symbol_group_id();
void append_vars(Vars &out, const Expr &expr);
void append_vars(Vars &out, const Vec &vec);
void append_params(Params &out, const Expr &expr);
void append_params(Params &out, const Vec &vec);
} // namespace detail

class Var {
public:
    Var() = default;

    VarId id() const;
    SymbolGroupId group_id() const;
    const std::string &label() const;
    int component() const;
    bool valid() const;

    bool operator==(const Var &other) const;
    bool operator!=(const Var &other) const;

private:
    Var(VarId id, std::string label, int component, SymbolGroupId group_id);

    VarId id_ = 0;
    std::string label_;
    int component_ = -1;
    SymbolGroupId group_id_ = 0;

    friend Var detail::make_var(std::string label, int component, SymbolGroupId group_id);
};

class Param {
public:
    Param() = default;

    ParamId id() const;
    SymbolGroupId group_id() const;
    const std::string &label() const;
    int component() const;
    bool valid() const;

    bool operator==(const Param &other) const;
    bool operator!=(const Param &other) const;

private:
    Param(ParamId id, std::string label, int component, SymbolGroupId group_id);

    ParamId id_ = 0;
    std::string label_;
    int component_ = -1;
    SymbolGroupId group_id_ = 0;

    friend Param detail::make_param(std::string label, int component, SymbolGroupId group_id);
};

class Vars {
public:
    Vars() = default;
    explicit Vars(std::vector<Var> values);

    int size() const;
    bool empty() const;
    const Var &operator[](int index) const;
    const std::vector<Var> &values() const;
    void append(const Var &var);

private:
    std::vector<Var> values_;
};

class Params {
public:
    Params() = default;
    explicit Params(std::vector<Param> values);

    int size() const;
    bool empty() const;
    const Param &operator[](int index) const;
    const std::vector<Param> &values() const;
    void append(const Param &param);

private:
    std::vector<Param> values_;
};

template <typename... Args>
Vars vars(const Args &...args) {
    Vars out;
    (detail::append_vars(out, args), ...);
    return out;
}

template <typename... Args>
Params params(const Args &...args) {
    Params out;
    (detail::append_params(out, args), ...);
    return out;
}

} // namespace ad

#endif // MOO_AD_SYMBOL_H
