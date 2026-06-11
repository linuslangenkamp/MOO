// SPDX-License-Identifier: LGPL-3.0-or-later
#include "vec.h"

#include "detail/node.h"
#include "detail/traversal.h"

#include <stdexcept>
#include <utility>

namespace ad {
namespace {

void require_valid(const Vec &vec) {
    if (!vec.valid()) {
        throw std::runtime_error("invalid vector expression");
    }
}

void require_same_size(const Vec &lhs, const Vec &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("vector sizes differ");
    }
}

void validate_indices_in_range(const std::vector<int> &indices, int upper_bound, const char *message) {
    for (int index : indices) {
        if (index < 0 || index >= upper_bound) {
            throw std::runtime_error(message);
        }
    }
}

Vec make_binary(detail::VecBinaryOp op, const Vec &lhs, const Vec &rhs) {
    require_same_size(lhs, rhs);
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorBinary;
    node->size = lhs.size();
    node->binary = op;
    node->lhs = detail::vec_node(lhs);
    node->rhs = detail::vec_node(rhs);
    return detail::vec_from_node(node);
}

Vec make_unary(detail::VecUnaryOp op, const Vec &arg) {
    require_valid(arg);
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorUnary;
    node->size = arg.size();
    node->unary = op;
    node->lhs = detail::vec_node(arg);
    return detail::vec_from_node(node);
}

} // namespace

Vec::Vec(std::vector<Expr> elements) {
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorFromElements;
    node->size = static_cast<int>(elements.size());
    node->elements.reserve(elements.size());
    for (const Expr &element : elements) {
        if (!element.valid()) {
            throw std::runtime_error("invalid scalar element in vector expression");
        }
        node->elements.push_back(detail::scalar_node(element));
    }
    node_ = node;
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
    require_valid(*this);
    if (index < 0 || index >= size()) {
        throw std::out_of_range("vector element index out of range");
    }
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::VectorElement;
    node->vec = node_;
    node->index = index;
    return detail::expr_from_node(node);
}

Vec Vec::slice(int start, int length) const {
    require_valid(*this);
    if (start < 0 || length < 0 || start + length > size()) {
        throw std::runtime_error("invalid vector slice");
    }
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::Slice;
    node->size = length;
    node->lhs = node_;
    node->start = start;
    return detail::vec_from_node(node);
}

Vars Vec::variables() const {
    return detail::collect_variables(*this);
}

Params Vec::parameters() const {
    return detail::collect_parameters(*this);
}

Vec vec(std::vector<Expr> elements) {
    return Vec(std::move(elements));
}

Vec vec_variable(const std::string &label, int size) {
    if (size < 0) {
        throw std::runtime_error("vector variable size must be non-negative");
    }
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorVariable;
    node->size = size;
    node->label = label;
    node->vars.reserve(static_cast<std::size_t>(size));
    const SymbolGroupId group_id = detail::next_symbol_group_id();
    for (int i = 0; i < size; ++i) {
        node->vars.push_back(detail::make_var(label, i, group_id));
    }
    return detail::vec_from_node(node);
}

Vec vec_parameter(const std::string &label, int size) {
    if (size < 0) {
        throw std::runtime_error("vector parameter size must be non-negative");
    }
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorParameter;
    node->size = size;
    node->label = label;
    node->params.reserve(static_cast<std::size_t>(size));
    const SymbolGroupId group_id = detail::next_symbol_group_id();
    for (int i = 0; i < size; ++i) {
        node->params.push_back(detail::make_param(label, i, group_id));
    }
    return detail::vec_from_node(node);
}

Vec vec_constant(const std::vector<double> &values) {
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorConstant;
    node->size = static_cast<int>(values.size());
    node->constants = values;
    return detail::vec_from_node(node);
}

Vec concat(const Vec &lhs, const Vec &rhs) {
    require_valid(lhs);
    require_valid(rhs);
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::Concat;
    node->size = lhs.size() + rhs.size();
    node->lhs = detail::vec_node(lhs);
    node->rhs = detail::vec_node(rhs);
    return detail::vec_from_node(node);
}

Vec gather(const Vec &source, std::vector<int> indices) {
    require_valid(source);
    validate_indices_in_range(indices, source.size(), "gather index is out of source vector range");
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::Gather;
    node->size = static_cast<int>(indices.size());
    node->lhs = detail::vec_node(source);
    node->indices = std::move(indices);
    return detail::vec_from_node(node);
}

Vec scatter_add(const Vec &values, std::vector<int> indices, int output_size) {
    require_valid(values);
    if (output_size < 0) {
        throw std::runtime_error("scatter_add output size must be non-negative");
    }
    if (static_cast<int>(indices.size()) != values.size()) {
        throw std::runtime_error("scatter_add index count must match value vector size");
    }
    validate_indices_in_range(indices, output_size, "scatter_add index is out of output vector range");
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::ScatterAdd;
    node->size = output_size;
    node->lhs = detail::vec_node(values);
    node->indices = std::move(indices);
    return detail::vec_from_node(node);
}

Expr sum(const Vec &values) {
    require_valid(values);
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::Sum;
    node->vec = detail::vec_node(values);
    return detail::expr_from_node(node);
}

Expr dot(const Vec &lhs, const Vec &rhs) {
    require_same_size(lhs, rhs);
    auto node = std::make_shared<detail::ScalarNode>();
    node->kind = GraphNodeKind::Dot;
    node->vec_lhs = detail::vec_node(lhs);
    node->vec_rhs = detail::vec_node(rhs);
    return detail::expr_from_node(node);
}

Vec operator+(const Vec &lhs, const Vec &rhs) {
    return make_binary(detail::VecBinaryOp::Add, lhs, rhs);
}

Vec operator-(const Vec &lhs, const Vec &rhs) {
    return make_binary(detail::VecBinaryOp::Sub, lhs, rhs);
}

Vec operator*(const Vec &lhs, const Vec &rhs) {
    return make_binary(detail::VecBinaryOp::Mul, lhs, rhs);
}

Vec operator/(const Vec &lhs, const Vec &rhs) {
    return make_binary(detail::VecBinaryOp::Div, lhs, rhs);
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
    if (!lhs.valid()) {
        throw std::runtime_error("invalid scalar scale expression");
    }
    require_valid(rhs);
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::VectorScale;
    node->size = rhs.size();
    node->scale = detail::scalar_node(lhs);
    node->lhs = detail::vec_node(rhs);
    return detail::vec_from_node(node);
}

Vec operator*(const Vec &lhs, const Expr &rhs) {
    return rhs * lhs;
}

Vec operator*(const DenseMatrix &matrix, const Vec &rhs) {
    require_valid(rhs);
    if (matrix.cols != rhs.size()) {
        throw std::runtime_error("dense matrix/vector dimensions differ");
    }
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::DenseMatVec;
    node->size = matrix.rows;
    node->dense = matrix;
    node->lhs = detail::vec_node(rhs);
    return detail::vec_from_node(node);
}

Vec operator*(const SparseMatrix &matrix, const Vec &rhs) {
    require_valid(rhs);
    if (matrix.cols != rhs.size()) {
        throw std::runtime_error("sparse matrix/vector dimensions differ");
    }
    auto node = std::make_shared<detail::VecNode>();
    node->kind = GraphNodeKind::SparseMatVec;
    node->size = matrix.rows;
    node->sparse = matrix;
    node->lhs = detail::vec_node(rhs);
    return detail::vec_from_node(node);
}

Vec sin(const Vec &expr) {
    return make_unary(detail::VecUnaryOp::Sin, expr);
}

Vec cos(const Vec &expr) {
    return make_unary(detail::VecUnaryOp::Cos, expr);
}

Vec tan(const Vec &expr) {
    return make_unary(detail::VecUnaryOp::Tan, expr);
}

Vec exp(const Vec &expr) {
    return make_unary(detail::VecUnaryOp::Exp, expr);
}

Vec log(const Vec &expr) {
    return make_unary(detail::VecUnaryOp::Log, expr);
}

Vec sigmoid(const Vec &expr) {
    return make_unary(detail::VecUnaryOp::Sigmoid, expr);
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
    }
    return "unknown";
}

} // namespace detail
} // namespace ad
