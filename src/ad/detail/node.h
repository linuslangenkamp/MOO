// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_DETAIL_NODE_H
#define MOO_AD_DETAIL_NODE_H

#include "../graph_info.h"
#include "../map_kind.h"
#include "../matrix.h"
#include "../symbol.h"

#include <atomic>
#include <memory>
#include <string>
#include <vector>

namespace ad::detail {

inline NodeId next_node_id() {
    static std::atomic<NodeId> next{1};
    return next.fetch_add(1);
}

enum class ScalarUnaryOp {
    Neg,
    Sin,
    Cos,
    Tan,
    Exp,
    Log,
    PowConst
};

enum class ScalarBinaryOp {
    Add,
    Sub,
    Mul,
    Div
};

enum class VecUnaryOp {
    Sin,
    Cos,
    Tan,
    Exp,
    Log,
    Sigmoid
};

enum class VecBinaryOp {
    Add,
    Sub,
    Mul,
    Div
};

struct VecNode;
struct FunctionCore;

struct MappedBindingNode {
    std::shared_ptr<const VecNode> local_input;
    std::shared_ptr<const VecNode> source;
    int reps = 0;
    int local_size = 0;
    std::vector<int> indices;
    MapKind map_kind = MapKind::Invalid;
    int base = 0;
    int rep_stride = 0;
    int component_stride = 1;
    int shift = 0;
    std::vector<int> offsets;
};

struct ScalarNode {
    GraphNodeKind kind = GraphNodeKind::Invalid;
    NodeId id = next_node_id();
    double value = 0.0;
    Var var;
    Param param;
    ScalarUnaryOp unary = ScalarUnaryOp::Neg;
    ScalarBinaryOp binary = ScalarBinaryOp::Add;
    std::shared_ptr<const ScalarNode> lhs;
    std::shared_ptr<const ScalarNode> rhs;
    std::shared_ptr<const VecNode> vec;
    std::shared_ptr<const VecNode> vec_lhs;
    std::shared_ptr<const VecNode> vec_rhs;
    int index = -1;
};

struct VecNode {
    GraphNodeKind kind = GraphNodeKind::Invalid;
    NodeId id = next_node_id();
    int size = 0;
    std::string label;
    std::vector<Var> vars;
    std::vector<Param> params;
    std::vector<double> constants;
    std::vector<std::shared_ptr<const ScalarNode>> elements;
    std::shared_ptr<const VecNode> lhs;
    std::shared_ptr<const VecNode> rhs;
    std::shared_ptr<const ScalarNode> scale;
    VecUnaryOp unary = VecUnaryOp::Sin;
    VecBinaryOp binary = VecBinaryOp::Add;
    DenseMatrix dense;
    SparseMatrix sparse;
    int start = 0;
    std::vector<int> indices;
    std::shared_ptr<const FunctionCore> function;
    std::vector<std::shared_ptr<const VecNode>> arguments;
    int reps = 0;
    std::vector<MappedBindingNode> mapped_bindings;
    MappedOutput mapped_output;
};

const char *to_string(ScalarUnaryOp op);
const char *to_string(ScalarBinaryOp op);
const char *to_string(VecUnaryOp op);
const char *to_string(VecBinaryOp op);

} // namespace ad::detail

#endif // MOO_AD_DETAIL_NODE_H
