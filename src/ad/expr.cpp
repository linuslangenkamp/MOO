// SPDX-License-Identifier: LGPL-3.0-or-later
#include "expr.h"

#include "detail/node.h"

#include <stdexcept>
#include <utility>

namespace ad {
namespace {

void require_valid(const Expr &expr) {
    if (!expr.valid()) {
        throw std::runtime_error("invalid scalar expression");
    }
}

Expr make_unary(detail::ScalarUnaryOp op, const Expr &arg) {
    require_valid(arg);
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::ScalarUnary;
    node->unary = op;
    node->lhs = detail::scalar_node(arg);
    return detail::expr_from_node(node);
}

Expr make_binary(detail::ScalarBinaryOp op, const Expr &lhs, const Expr &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::ScalarBinary;
    node->binary = op;
    node->lhs = detail::scalar_node(lhs);
    node->rhs = detail::scalar_node(rhs);
    return detail::expr_from_node(node);
}

} // namespace

Expr::Expr(std::shared_ptr<const detail::ScalarNode> node)
    : node_(std::move(node)) {}

bool Expr::valid() const {
    return static_cast<bool>(node_);
}

NodeId Expr::node_id() const {
    return node_ ? node_->id : 0;
}

GraphNodeKind Expr::node_kind() const {
    return node_ ? node_->kind : GraphNodeKind::Invalid;
}

bool Expr::is_constant(double *value) const {
    if (!node_ || node_->kind != GraphNodeKind::ScalarConstant) {
        return false;
    }
    if (value) {
        *value = node_->value;
    }
    return true;
}

bool Expr::is_variable() const {
    return node_ && node_->kind == GraphNodeKind::ScalarVariable;
}

bool Expr::is_parameter() const {
    return node_ && node_->kind == GraphNodeKind::ScalarParameter;
}

Var Expr::var() const {
    if (!is_variable()) {
        throw std::runtime_error("scalar expression is not a variable");
    }
    return node_->var;
}

Param Expr::param() const {
    if (!is_parameter()) {
        throw std::runtime_error("scalar expression is not a parameter");
    }
    return node_->param;
}

Expr constant(double value) {
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::ScalarConstant;
    node->value = value;
    return detail::expr_from_node(node);
}

Expr variable(const std::string &label) {
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::ScalarVariable;
    node->var = detail::make_var(label, -1, detail::next_symbol_group_id());
    return detail::expr_from_node(node);
}

Expr parameter(const std::string &label) {
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::ScalarParameter;
    node->param = detail::make_param(label, -1, detail::next_symbol_group_id());
    return detail::expr_from_node(node);
}

Expr operator+(const Expr &lhs, const Expr &rhs) {
    return make_binary(detail::ScalarBinaryOp::Add, lhs, rhs);
}

Expr operator-(const Expr &lhs, const Expr &rhs) {
    return make_binary(detail::ScalarBinaryOp::Sub, lhs, rhs);
}

Expr operator*(const Expr &lhs, const Expr &rhs) {
    return make_binary(detail::ScalarBinaryOp::Mul, lhs, rhs);
}

Expr operator/(const Expr &lhs, const Expr &rhs) {
    return make_binary(detail::ScalarBinaryOp::Div, lhs, rhs);
}

Expr operator-(const Expr &expr) {
    return make_unary(detail::ScalarUnaryOp::Neg, expr);
}

Expr operator+(const Expr &lhs, double rhs) {
    return lhs + constant(rhs);
}

Expr operator-(const Expr &lhs, double rhs) {
    return lhs - constant(rhs);
}

Expr operator*(const Expr &lhs, double rhs) {
    return lhs * constant(rhs);
}

Expr operator/(const Expr &lhs, double rhs) {
    return lhs / constant(rhs);
}

Expr operator+(double lhs, const Expr &rhs) {
    return constant(lhs) + rhs;
}

Expr operator-(double lhs, const Expr &rhs) {
    return constant(lhs) - rhs;
}

Expr operator*(double lhs, const Expr &rhs) {
    return constant(lhs) * rhs;
}

Expr operator/(double lhs, const Expr &rhs) {
    return constant(lhs) / rhs;
}

Expr sin(const Expr &expr) {
    return make_unary(detail::ScalarUnaryOp::Sin, expr);
}

Expr cos(const Expr &expr) {
    return make_unary(detail::ScalarUnaryOp::Cos, expr);
}

Expr tan(const Expr &expr) {
    return make_unary(detail::ScalarUnaryOp::Tan, expr);
}

Expr exp(const Expr &expr) {
    return make_unary(detail::ScalarUnaryOp::Exp, expr);
}

Expr log(const Expr &expr) {
    return make_unary(detail::ScalarUnaryOp::Log, expr);
}

Expr pow(const Expr &expr, double exponent) {
    require_valid(expr);
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::ScalarUnary;
    node->unary = detail::ScalarUnaryOp::PowConst;
    node->lhs = detail::scalar_node(expr);
    node->value = exponent;
    return detail::expr_from_node(node);
}

namespace detail {

const std::shared_ptr<const ScalarNode> &scalar_node(const Expr &expr) {
    return expr.node_;
}

Expr expr_from_node(std::shared_ptr<const ScalarNode> node) {
    return Expr(std::move(node));
}

void append_vars(Vars &out, const Expr &expr) {
    if (expr.is_variable()) {
        out.append(expr.var());
    }
}

void append_params(Params &out, const Expr &expr) {
    if (expr.is_parameter()) {
        out.append(expr.param());
    }
}

const char *to_string(ScalarUnaryOp op) {
    switch (op) {
        case ScalarUnaryOp::Neg:
            return "neg";
        case ScalarUnaryOp::Sin:
            return "sin";
        case ScalarUnaryOp::Cos:
            return "cos";
        case ScalarUnaryOp::Tan:
            return "tan";
        case ScalarUnaryOp::Exp:
            return "exp";
        case ScalarUnaryOp::Log:
            return "log";
        case ScalarUnaryOp::PowConst:
            return "pow_const";
    }
    return "unknown";
}

const char *to_string(ScalarBinaryOp op) {
    switch (op) {
        case ScalarBinaryOp::Add:
            return "add";
        case ScalarBinaryOp::Sub:
            return "sub";
        case ScalarBinaryOp::Mul:
            return "mul";
        case ScalarBinaryOp::Div:
            return "div";
    }
    return "unknown";
}

} // namespace detail
} // namespace ad
