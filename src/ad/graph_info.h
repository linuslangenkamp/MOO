// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_GRAPH_INFO_H
#define MOO_AD_GRAPH_INFO_H

#include <cstdint>
#include <string>
#include <vector>

namespace ad {

using NodeId = std::uint64_t;

class Expr;
class Vec;

enum class GraphNodeKind {
    Invalid,
    ScalarConstant,
    ScalarVariable,
    ScalarParameter,
    ScalarUnary,
    ScalarBinary,
    VectorElement,
    Sum,
    Dot,
    VectorVariable,
    VectorParameter,
    VectorConstant,
    VectorFromElements,
    VectorUnary,
    VectorBinary,
    VectorScale,
    DenseMatVec,
    SparseMatVec,
    SymbolicMatVec,
    SymbolicMatMul,
    SymbolicSparseMatVec,
    SymbolicSparseMatMul,
    OuterProduct,
    LinearSolve,
    Slice,
    ScatterSlice,
    Gather,
    ScatterAdd,
    Concat,
    FunctionCall,
    MappedFunctionCall,
    MapAccumCall
};

struct GraphInfo {
    GraphNodeKind kind = GraphNodeKind::Invalid;
    NodeId id = 0;
    int size = 0;
    std::string op;
    std::vector<GraphInfo> children;
};

GraphInfo inspect(const Expr &expr);
GraphInfo inspect(const Vec &vec);

const char *to_string(GraphNodeKind kind);
int count_nodes(const GraphInfo &info, GraphNodeKind kind);

} // namespace ad

#endif // MOO_AD_GRAPH_INFO_H
