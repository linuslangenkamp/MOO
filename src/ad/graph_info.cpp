// SPDX-License-Identifier: LGPL-3.0-or-later
#include "graph_info.h"

#include "expr.h"
#include "vec.h"
#include "detail/node.h"

#include <stdexcept>

namespace ad {
namespace {

GraphInfo inspect_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node);
GraphInfo inspect_vec_node(const std::shared_ptr<const detail::VecNode> &node);

GraphInfo inspect_scalar_node(const std::shared_ptr<const detail::ScalarNode> &node) {
    if (!node) {
        return {};
    }

    GraphInfo info;
    info.kind = node->kind;
    info.id = node->id;
    info.size = 1;

    switch (node->kind) {
        case GraphNodeKind::ScalarUnary:
            info.op = detail::to_string(node->unary);
            info.children.push_back(inspect_scalar_node(node->lhs));
            break;
        case GraphNodeKind::ScalarBinary:
            info.op = detail::to_string(node->binary);
            info.children.push_back(inspect_scalar_node(node->lhs));
            info.children.push_back(inspect_scalar_node(node->rhs));
            break;
        case GraphNodeKind::VectorElement:
            info.op = std::to_string(node->index);
            info.children.push_back(inspect_vec_node(node->vec));
            break;
        case GraphNodeKind::Sum:
            info.op = "sum";
            info.children.push_back(inspect_vec_node(node->vec));
            break;
        case GraphNodeKind::Dot:
            info.op = "dot";
            info.children.push_back(inspect_vec_node(node->vec_lhs));
            info.children.push_back(inspect_vec_node(node->vec_rhs));
            break;
        default:
            break;
    }

    return info;
}

GraphInfo inspect_vec_node(const std::shared_ptr<const detail::VecNode> &node) {
    if (!node) {
        return {};
    }

    GraphInfo info;
    info.kind = node->kind;
    info.id = node->id;
    info.size = node->size;

    switch (node->kind) {
        case GraphNodeKind::VectorFromElements:
            info.op = "elements";
            for (const auto &element : node->elements) {
                info.children.push_back(inspect_scalar_node(element));
            }
            break;
        case GraphNodeKind::VectorUnary:
            info.op = detail::to_string(node->unary);
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::VectorBinary:
            info.op = detail::to_string(node->binary);
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::VectorScale:
            info.op = "scale";
            info.children.push_back(inspect_scalar_node(node->scale));
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::DenseMatVec:
            info.op = "dense_matvec";
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::SparseMatVec:
            info.op = "sparse_matvec";
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::SymbolicMatVec:
            info.op = "symbolic_matvec";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::SymbolicMatMul:
            info.op = "symbolic_matmul";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::SymbolicSparseMatVec:
            info.op = "symbolic_sparse_matvec";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::SymbolicSparseMatMul:
            info.op = node->symbolic_sparse_lhs ? "symbolic_sparse_dense_matmul" : "symbolic_dense_sparse_matmul";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::OuterProduct:
            info.op = "outer_product";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::LinearSolve:
            info.op = node->linear_solve_transpose ? "solve_transpose" : "solve";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::Slice:
            info.op = std::to_string(node->start);
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::ScatterSlice:
            info.op = std::to_string(node->start);
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::Gather:
            info.op = "gather";
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::ScatterAdd:
            info.op = "scatter_add";
            info.children.push_back(inspect_vec_node(node->lhs));
            break;
        case GraphNodeKind::Concat:
            info.op = "concat";
            info.children.push_back(inspect_vec_node(node->lhs));
            info.children.push_back(inspect_vec_node(node->rhs));
            break;
        case GraphNodeKind::FunctionCall:
            info.op = "function_call";
            for (const auto &argument : node->arguments) {
                info.children.push_back(inspect_vec_node(argument));
            }
            break;
        case GraphNodeKind::MappedFunctionCall:
            info.op = "mapped_function_call";
            for (const auto &binding : node->mapped_bindings) {
                info.children.push_back(inspect_vec_node(binding.source));
            }
            break;
        case GraphNodeKind::MapAccumCall:
            info.op = "map_accum_call";
            info.children.push_back(inspect_vec_node(node->lhs));
            for (const auto &binding : node->mapped_bindings) {
                info.children.push_back(inspect_vec_node(binding.source));
            }
            break;
        default:
            break;
    }

    return info;
}

} // namespace

GraphInfo inspect(const Expr &expr) {
    if (!expr.valid()) {
        throw std::runtime_error("cannot inspect invalid scalar expression");
    }
    return inspect_scalar_node(detail::scalar_node(expr));
}

GraphInfo inspect(const Vec &vec) {
    if (!vec.valid()) {
        throw std::runtime_error("cannot inspect invalid vector expression");
    }
    return inspect_vec_node(detail::vec_node(vec));
}

const char *to_string(GraphNodeKind kind) {
    switch (kind) {
        case GraphNodeKind::Invalid:
            return "Invalid";
        case GraphNodeKind::ScalarConstant:
            return "ScalarConstant";
        case GraphNodeKind::ScalarVariable:
            return "ScalarVariable";
        case GraphNodeKind::ScalarParameter:
            return "ScalarParameter";
        case GraphNodeKind::ScalarUnary:
            return "ScalarUnary";
        case GraphNodeKind::ScalarBinary:
            return "ScalarBinary";
        case GraphNodeKind::VectorElement:
            return "VectorElement";
        case GraphNodeKind::Sum:
            return "Sum";
        case GraphNodeKind::Dot:
            return "Dot";
        case GraphNodeKind::VectorVariable:
            return "VectorVariable";
        case GraphNodeKind::VectorParameter:
            return "VectorParameter";
        case GraphNodeKind::VectorConstant:
            return "VectorConstant";
        case GraphNodeKind::VectorFromElements:
            return "VectorFromElements";
        case GraphNodeKind::VectorUnary:
            return "VectorUnary";
        case GraphNodeKind::VectorBinary:
            return "VectorBinary";
        case GraphNodeKind::VectorScale:
            return "VectorScale";
        case GraphNodeKind::DenseMatVec:
            return "DenseMatVec";
        case GraphNodeKind::SparseMatVec:
            return "SparseMatVec";
        case GraphNodeKind::SymbolicMatVec:
            return "SymbolicMatVec";
        case GraphNodeKind::SymbolicMatMul:
            return "SymbolicMatMul";
        case GraphNodeKind::SymbolicSparseMatVec:
            return "SymbolicSparseMatVec";
        case GraphNodeKind::SymbolicSparseMatMul:
            return "SymbolicSparseMatMul";
        case GraphNodeKind::OuterProduct:
            return "OuterProduct";
        case GraphNodeKind::LinearSolve:
            return "LinearSolve";
        case GraphNodeKind::Slice:
            return "Slice";
        case GraphNodeKind::ScatterSlice:
            return "ScatterSlice";
        case GraphNodeKind::Gather:
            return "Gather";
        case GraphNodeKind::ScatterAdd:
            return "ScatterAdd";
        case GraphNodeKind::Concat:
            return "Concat";
        case GraphNodeKind::FunctionCall:
            return "FunctionCall";
        case GraphNodeKind::MappedFunctionCall:
            return "MappedFunctionCall";
        case GraphNodeKind::MapAccumCall:
            return "MapAccumCall";
    }
    return "Unknown";
}

int count_nodes(const GraphInfo &info, GraphNodeKind kind) {
    int count = info.kind == kind ? 1 : 0;
    for (const GraphInfo &child : info.children) {
        count += count_nodes(child, kind);
    }
    return count;
}

} // namespace ad
