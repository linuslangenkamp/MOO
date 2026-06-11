// SPDX-License-Identifier: LGPL-3.0-or-later
#include "vec.h"

#include "detail/node.h"
#include "detail/simplify.h"
#include "detail/traversal.h"
#include <utility>

namespace ad {

Vec::Vec(std::vector<Expr> elements) {
    node_ = detail::vec_node(detail::make_vec_from_elements(std::move(elements)));
}

Vec::Vec(std::initializer_list<Expr> elements)
    : Vec(std::vector<Expr>(elements)) {}

Vec::Vec(std::shared_ptr<const detail::VecNode> node)
    : node_(std::move(node)) {}

bool Vec::valid() const {
    return static_cast<bool>(node_);
}

bool Vec::empty() const {
    return size() == 0;
}

int Vec::size() const {
    return node_ ? node_->size : 0;
}

NodeId Vec::node_id() const {
    return node_ ? node_->id : 0;
}

GraphNodeKind Vec::node_kind() const {
    return node_ ? node_->kind : GraphNodeKind::Invalid;
}

Expr Vec::operator[](int index) const {
    return detail::make_vector_element(*this, index);
}

Vec Vec::slice(int start, int length) const {
    return detail::make_slice(*this, start, length);
}

Vars Vec::variables() const {
    return detail::collect_variables(*this);
}

Params Vec::parameters() const {
    return detail::collect_parameters(*this);
}

Vec vec(std::vector<Expr> elements) {
    return detail::make_vec_from_elements(std::move(elements));
}

Vec vec_variable(const std::string &label, int size) {
    return detail::make_vec_variable(label, size);
}

Vec vec_parameter(const std::string &label, int size) {
    return detail::make_vec_parameter(label, size);
}

Vec vec_constant(const std::vector<double> &values) {
    return detail::make_vec_constant(values);
}

Vec concat(const Vec &lhs, const Vec &rhs) {
    return detail::make_concat(lhs, rhs);
}

Vec gather(const Vec &source, std::vector<int> indices) {
    return detail::make_gather(source, std::move(indices));
}

Vec scatter_add(const Vec &values, std::vector<int> indices, int output_size) {
    return detail::make_scatter_add(values, std::move(indices), output_size);
}

Expr sum(const Vec &values) {
    return detail::make_sum(values);
}

Expr dot(const Vec &lhs, const Vec &rhs) {
    return detail::make_dot(lhs, rhs);
}

Vec operator+(const Vec &lhs, const Vec &rhs) {
    return detail::make_vec_binary(detail::VecBinaryOp::Add, lhs, rhs);
}

Vec operator-(const Vec &lhs, const Vec &rhs) {
    return detail::make_vec_binary(detail::VecBinaryOp::Sub, lhs, rhs);
}

Vec operator*(const Vec &lhs, const Vec &rhs) {
    return detail::make_vec_binary(detail::VecBinaryOp::Mul, lhs, rhs);
}

Vec operator/(const Vec &lhs, const Vec &rhs) {
    return detail::make_vec_binary(detail::VecBinaryOp::Div, lhs, rhs);
}

Vec operator+(const Vec &lhs, double rhs) {
    return lhs + vec_constant(std::vector<double>(static_cast<std::size_t>(lhs.size()), rhs));
}

Vec operator-(const Vec &lhs, double rhs) {
    return lhs - vec_constant(std::vector<double>(static_cast<std::size_t>(lhs.size()), rhs));
}

Vec operator*(const Vec &lhs, double rhs) {
    return constant(rhs) * lhs;
}

Vec operator/(const Vec &lhs, double rhs) {
    return constant(1.0 / rhs) * lhs;
}

Vec operator+(double lhs, const Vec &rhs) {
    return vec_constant(std::vector<double>(static_cast<std::size_t>(rhs.size()), lhs)) + rhs;
}

Vec operator-(double lhs, const Vec &rhs) {
    return vec_constant(std::vector<double>(static_cast<std::size_t>(rhs.size()), lhs)) - rhs;
}

Vec operator*(double lhs, const Vec &rhs) {
    return constant(lhs) * rhs;
}

Vec operator/(double lhs, const Vec &rhs) {
    return vec_constant(std::vector<double>(static_cast<std::size_t>(rhs.size()), lhs)) / rhs;
}

Vec operator*(const Expr &lhs, const Vec &rhs) {
    return detail::make_vec_scale(lhs, rhs);
}

Vec operator*(const Vec &lhs, const Expr &rhs) {
    return rhs * lhs;
}

Vec operator-(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Neg, expr);
}

Vec operator*(const DenseMatrix &matrix, const Vec &rhs) {
    return detail::make_dense_matvec(matrix, rhs);
}

Vec operator*(const SparseMatrix &matrix, const Vec &rhs) {
    return detail::make_sparse_matvec(matrix, rhs);
}

Vec sin(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Sin, expr);
}

Vec cos(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Cos, expr);
}

Vec tan(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Tan, expr);
}

Vec exp(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Exp, expr);
}

Vec log(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Log, expr);
}

Vec sigmoid(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Sigmoid, expr);
}

Vec abs(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Abs, expr);
}

Vec sqrt(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Sqrt, expr);
}

Vec asin(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Asin, expr);
}

Vec acos(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Acos, expr);
}

Vec atan(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Atan, expr);
}

Vec sinh(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Sinh, expr);
}

Vec cosh(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Cosh, expr);
}

Vec tanh(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Tanh, expr);
}

Vec log10(const Vec &expr) {
    return detail::make_vec_unary(detail::VecUnaryOp::Log10, expr);
}

Vec pow(const Vec &base, const Vec &exponent) {
    return detail::make_vec_binary(detail::VecBinaryOp::Pow, base, exponent);
}

Vec pow(const Vec &base, double exponent) {
    return pow(base, vec_constant(std::vector<double>(static_cast<std::size_t>(base.size()), exponent)));
}

Vec min(const Vec &lhs, const Vec &rhs) {
    return detail::make_vec_binary(detail::VecBinaryOp::Min, lhs, rhs);
}

Vec max(const Vec &lhs, const Vec &rhs) {
    return detail::make_vec_binary(detail::VecBinaryOp::Max, lhs, rhs);
}

namespace detail {

const std::shared_ptr<const VecNode> &vec_node(const Vec &vec) {
    return vec.node_;
}

Vec vec_from_node(std::shared_ptr<const VecNode> node) {
    return Vec(std::move(node));
}

void append_vars(Vars &out, const Vec &vec) {
    Vars collected = vec.variables();
    for (const Var &var : collected.values()) {
        out.append(var);
    }
}

void append_params(Params &out, const Vec &vec) {
    Params collected = vec.parameters();
    for (const Param &param : collected.values()) {
        out.append(param);
    }
}

const char *to_string(VecUnaryOp op) {
    switch (op) {
        case VecUnaryOp::Neg:
            return "neg";
        case VecUnaryOp::Sin:
            return "sin";
        case VecUnaryOp::Cos:
            return "cos";
        case VecUnaryOp::Tan:
            return "tan";
        case VecUnaryOp::Exp:
            return "exp";
        case VecUnaryOp::Log:
            return "log";
        case VecUnaryOp::Sigmoid:
            return "sigmoid";
        case VecUnaryOp::Abs:
            return "abs";
        case VecUnaryOp::Sqrt:
            return "sqrt";
        case VecUnaryOp::Asin:
            return "asin";
        case VecUnaryOp::Acos:
            return "acos";
        case VecUnaryOp::Atan:
            return "atan";
        case VecUnaryOp::Sinh:
            return "sinh";
        case VecUnaryOp::Cosh:
            return "cosh";
        case VecUnaryOp::Tanh:
            return "tanh";
        case VecUnaryOp::Log10:
            return "log10";
    }
    return "unknown";
}

const char *to_string(VecBinaryOp op) {
    switch (op) {
        case VecBinaryOp::Add:
            return "add";
        case VecBinaryOp::Sub:
            return "sub";
        case VecBinaryOp::Mul:
            return "mul";
        case VecBinaryOp::Div:
            return "div";
        case VecBinaryOp::Pow:
            return "pow";
        case VecBinaryOp::Min:
            return "min";
        case VecBinaryOp::Max:
            return "max";
    }
    return "unknown";
}

} // namespace detail
} // namespace ad
