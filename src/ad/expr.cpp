// SPDX-License-Identifier: LGPL-3.0-or-later
#include "expr.h"

#include "detail/node.h"
#include "detail/simplify.h"

#include <stdexcept>
#include <utility>

namespace ad {

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
    return detail::make_scalar_constant(value);
}

Expr variable(const std::string &label) {
    return detail::make_scalar_variable(label);
}

Expr parameter(const std::string &label) {
    return detail::make_scalar_parameter(label);
}

Expr operator+(const Expr &lhs, const Expr &rhs) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Add, lhs, rhs);
}

Expr operator-(const Expr &lhs, const Expr &rhs) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Sub, lhs, rhs);
}

Expr operator*(const Expr &lhs, const Expr &rhs) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Mul, lhs, rhs);
}

Expr operator/(const Expr &lhs, const Expr &rhs) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Div, lhs, rhs);
}

Expr operator-(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Neg, expr);
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
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Sin, expr);
}

Expr cos(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Cos, expr);
}

Expr tan(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Tan, expr);
}

Expr exp(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Exp, expr);
}

Expr log(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Log, expr);
}

Expr abs(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Abs, expr);
}

Expr sqrt(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Sqrt, expr);
}

Expr asin(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Asin, expr);
}

Expr acos(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Acos, expr);
}

Expr atan(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Atan, expr);
}

Expr sinh(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Sinh, expr);
}

Expr cosh(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Cosh, expr);
}

Expr tanh(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Tanh, expr);
}

Expr log10(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Log10, expr);
}

Expr sigmoid(const Expr &expr) {
    return detail::make_scalar_unary(detail::ScalarUnaryOp::Sigmoid, expr);
}

Expr pow(const Expr &expr, double exponent) {
    return detail::make_scalar_pow_const(expr, exponent);
}

Expr pow(const Expr &base, const Expr &exponent) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Pow, base, exponent);
}

Expr min(const Expr &lhs, const Expr &rhs) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Min, lhs, rhs);
}

Expr max(const Expr &lhs, const Expr &rhs) {
    return detail::make_scalar_binary(detail::ScalarBinaryOp::Max, lhs, rhs);
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
        case ScalarUnaryOp::Abs:
            return "abs";
        case ScalarUnaryOp::Sqrt:
            return "sqrt";
        case ScalarUnaryOp::Asin:
            return "asin";
        case ScalarUnaryOp::Acos:
            return "acos";
        case ScalarUnaryOp::Atan:
            return "atan";
        case ScalarUnaryOp::Sinh:
            return "sinh";
        case ScalarUnaryOp::Cosh:
            return "cosh";
        case ScalarUnaryOp::Tanh:
            return "tanh";
        case ScalarUnaryOp::Log10:
            return "log10";
        case ScalarUnaryOp::Sigmoid:
            return "sigmoid";
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
        case ScalarBinaryOp::Pow:
            return "pow";
        case ScalarBinaryOp::Min:
            return "min";
        case ScalarBinaryOp::Max:
            return "max";
    }
    return "unknown";
}

} // namespace detail
} // namespace ad
