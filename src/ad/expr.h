// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_EXPR_H
#define MOO_AD_EXPR_H

#include "graph_info.h"
#include "symbol.h"

#include <memory>
#include <string>

namespace ad {

class Vec;
class Values;
class EvalWorkspace;

namespace detail {
struct ScalarNode;
const std::shared_ptr<const ScalarNode> &scalar_node(const Expr &expr);
Expr expr_from_node(std::shared_ptr<const ScalarNode> node);
} // namespace detail

class Expr {
public:
    Expr() = default;

    bool valid() const;
    NodeId node_id() const;
    GraphNodeKind node_kind() const;
    bool is_constant(double *value = nullptr) const;
    bool is_variable() const;
    bool is_parameter() const;
    Var var() const;
    Param param() const;
    void eval(const Values &values, EvalWorkspace &workspace, double *out) const;
    Expr forward_diff(const Vars &wrt, const Vec &seed) const;
    Vec reverse_diff(const Vars &wrt) const;

private:
    explicit Expr(std::shared_ptr<const detail::ScalarNode> node);

    std::shared_ptr<const detail::ScalarNode> node_;

    friend const std::shared_ptr<const detail::ScalarNode> &detail::scalar_node(const Expr &expr);
    friend Expr detail::expr_from_node(std::shared_ptr<const detail::ScalarNode> node);
};

Expr constant(double value);
Expr variable(const std::string &label);
Expr parameter(const std::string &label);

Expr operator+(const Expr &lhs, const Expr &rhs);
Expr operator-(const Expr &lhs, const Expr &rhs);
Expr operator*(const Expr &lhs, const Expr &rhs);
Expr operator/(const Expr &lhs, const Expr &rhs);
Expr operator-(const Expr &expr);

Expr operator+(const Expr &lhs, double rhs);
Expr operator-(const Expr &lhs, double rhs);
Expr operator*(const Expr &lhs, double rhs);
Expr operator/(const Expr &lhs, double rhs);
Expr operator+(double lhs, const Expr &rhs);
Expr operator-(double lhs, const Expr &rhs);
Expr operator*(double lhs, const Expr &rhs);
Expr operator/(double lhs, const Expr &rhs);

Expr sin(const Expr &expr);
Expr cos(const Expr &expr);
Expr tan(const Expr &expr);
Expr exp(const Expr &expr);
Expr log(const Expr &expr);
Expr pow(const Expr &expr, double exponent);

} // namespace ad

#endif // MOO_AD_EXPR_H
