// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
#include <string>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool test_scalar_local_simplification() {
    ad::Expr x = ad::variable("x");
    ad::Expr c = ad::sin(ad::constant(0.0)) + (x * 1.0) + (0.0 * x) - 0.0;
    ad::GraphInfo info = ad::inspect(c);

    bool ok = true;
    ok &= check(c.node_id() == x.node_id(), "scalar identities should simplify back to x");
    ok &= check(info.kind == ad::GraphNodeKind::ScalarVariable, "scalar simplification should remove identity operations");
    ok &= check(ad::pow(x, 1.0).node_id() == x.node_id(), "pow(x, 1) should simplify to x");
    ok &= check(ad::pow(x, 0.0).node_kind() == ad::GraphNodeKind::ScalarConstant, "pow(x, 0) should simplify to constant one");
    ok &= check((-(-x)).node_id() == x.node_id(), "double negation should simplify");
    return ok;
}

bool test_vector_local_simplification() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec zero = ad::vec_constant({0.0, 0.0, 0.0});
    ad::Vec one = ad::vec_constant({1.0, 1.0, 1.0});

    ad::Vec add_zero = X + zero;
    ad::Vec mul_one = X * one;
    ad::Vec scale_zero = ad::constant(0.0) * X;
    ad::Vec scale_one = ad::constant(1.0) * X;

    bool ok = true;
    ok &= check(add_zero.node_id() == X.node_id(), "X + 0 should simplify to X");
    ok &= check(mul_one.node_id() == X.node_id(), "X * 1 should simplify to X");
    ok &= check(scale_zero.node_kind() == ad::GraphNodeKind::VectorConstant, "0 * X should simplify to zero vector");
    ok &= check(scale_one.node_id() == X.node_id(), "1 * X should simplify to X");
    ok &= check(ad::concat(ad::vec_constant({}), X).node_id() == X.node_id(), "concat(empty, X) should simplify to X");
    ok &= check(ad::gather(X, {0, 1, 2}).node_id() == X.node_id(), "identity gather should simplify to X");
    return ok;
}

bool test_matvec_and_reductions_simplify_known_zero() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec zero = ad::vec_constant({0.0, 0.0, 0.0});
    ad::DenseMatrix A(2, 3, {1.0, 2.0, 0.0,
                             0.0, 0.0, 3.0});
    ad::SparseMatrix S(2, 3, {0, 1}, {1, 2}, {4.0, 5.0});

    bool ok = true;
    ok &= check((A * zero).node_kind() == ad::GraphNodeKind::VectorConstant, "dense matvec of zero should simplify to zero vector");
    ok &= check((S * zero).node_kind() == ad::GraphNodeKind::VectorConstant, "sparse matvec of zero should simplify to zero vector");
    ok &= check(ad::sum(zero).node_kind() == ad::GraphNodeKind::ScalarConstant, "sum(zero) should simplify to zero scalar");
    ok &= check(ad::dot(X, zero).node_kind() == ad::GraphNodeKind::ScalarConstant, "dot(X, zero) should simplify to zero scalar");
    return ok;
}

bool test_function_construction_simplifies_output() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec output = ad::constant(1.0) * (X + ad::vec_constant({0.0, 0.0}));
    ad::Function f({X}, output);
    ad::GraphInfo info = f.info().output_graph;

    bool ok = true;
    ok &= check(info.kind == ad::GraphNodeKind::VectorVariable, "Function should store simplified output graph");
    ok &= check(f.output().node_id() == X.node_id(), "Function output should be the simplified vector handle");
    ok &= check(f.jacobian_sparsity().nnz() == 2, "simplified identity function sparsity should remain correct");
    return ok;
}

bool test_derivative_simplification_removes_zero_branches() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec dX = ad::vec_variable("dX", 2);
    ad::Expr h = ad::parameter("h");

    ad::Vec y = h * X;
    ad::Vec dy = y.forward_diff(ad::vars(X), dX);
    ad::Vec lambda = ad::vec_variable("lambda", 2);
    ad::Vec adj = y.reverse_diff(ad::vars(X), lambda);

    bool ok = true;
    ok &= check(ad::inspect(dy).kind == ad::GraphNodeKind::VectorScale, "forward h * X should simplify to h * dX");
    ok &= check(ad::count_nodes(ad::inspect(dy), ad::GraphNodeKind::VectorBinary) == 0, "forward zero branch should be removed");
    ok &= check(ad::inspect(adj).kind == ad::GraphNodeKind::VectorScale, "reverse h * X should simplify to h * lambda");
    ok &= check(ad::count_nodes(ad::inspect(adj), ad::GraphNodeKind::VectorBinary) == 0, "reverse zero branch should be removed");
    return ok;
}

bool test_linear_hvp_simplifies_to_zero() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});
    ad::Vec y = A * X;
    ad::Vec lambda = ad::vec_variable("lambda", y.size());
    ad::Vec seed = ad::vec_variable("seed", X.size());

    ad::Vec grad = y.reverse_diff(ad::vars(X), lambda);
    ad::Vec hvp = grad.forward_diff(ad::vars(X), seed);
    ad::GraphInfo info = ad::inspect(hvp);

    bool ok = true;
    ok &= check(hvp.size() == X.size(), "linear HVP size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorConstant, "direct linear HVP should simplify to zero");
    return ok;
}

bool test_call_and_mapped_boundaries_survive_simplification() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x + ad::vec_constant({0.0, 0.0})));

    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec call = local.call({X.slice(0, 2) + ad::vec_constant({0.0, 0.0})});
    ad::Vec mapped = ad::map(local, {ad::bind(x, ad::Map::blocks(X + ad::vec_constant({0.0, 0.0, 0.0, 0.0}), 2, 2))});

    bool ok = true;
    ok &= check(ad::inspect(call).kind == ad::GraphNodeKind::FunctionCall, "FunctionCall boundary should survive simplification");
    ok &= check(contains_kind(ad::inspect(call), ad::GraphNodeKind::Slice), "simplified call argument should remain structured");
    ok &= check(ad::inspect(mapped).kind == ad::GraphNodeKind::MappedFunctionCall, "MappedFunctionCall boundary should survive simplification");
    ok &= check(!contains_kind(ad::inspect(mapped), ad::GraphNodeKind::VectorBinary), "mapped zero source branch should be removed");
    return ok;
}

bool test_explicit_optimize_api() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Expr e = ad::optimize((X[0] + ad::constant(0.0)) * ad::constant(1.0));
    ad::Vec v = ad::optimize(X + ad::vec_constant({0.0, 0.0}));
    ad::Function f({X}, ad::constant(1.0) * X);
    ad::Function g = ad::optimize(f);

    bool ok = true;
    ok &= check(e.node_kind() == ad::GraphNodeKind::VectorElement, "optimize(Expr) should simplify scalar identities");
    ok &= check(v.node_id() == X.node_id(), "optimize(Vec) should simplify vector identities");
    ok &= check(g.output().node_id() == X.node_id(), "optimize(Function) should simplify function output");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_local_simplification();
    ok &= test_vector_local_simplification();
    ok &= test_matvec_and_reductions_simplify_known_zero();
    ok &= test_function_construction_simplifies_output();
    ok &= test_derivative_simplification_removes_zero_branches();
    ok &= test_linear_hvp_simplifies_to_zero();
    ok &= test_call_and_mapped_boundaries_survive_simplification();
    ok &= test_explicit_optimize_api();

    if (!ok) {
        return 1;
    }

    std::cout << "AD simplification tests passed\n";
    return 0;
}
