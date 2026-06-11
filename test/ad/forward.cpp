// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

template <typename Fn>
bool throws_cleanly(Fn &&fn) {
    try {
        fn();
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool contains_op(const ad::GraphInfo &info, ad::GraphNodeKind kind, const std::string &op) {
    if (info.kind == kind && info.op == op) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_op(child, kind, op);
    });
}

bool test_scalar_forward_basics() {
    ad::Expr x = ad::variable("x");
    ad::Expr y = ad::variable("y");
    ad::Expr f = ad::sin(x) + ad::tan(y) + x * y;
    ad::Vec seed = ad::vec_variable("seed", 2);

    ad::Expr df = f.forward_diff(ad::vars(x, y), seed);
    ad::GraphInfo info = ad::inspect(df);

    bool ok = true;
    ok &= check(df.valid(), "scalar forward result should be valid");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorElement), "scalar forward should reference seed components symbolically");
    ok &= check(contains_op(info, ad::GraphNodeKind::ScalarUnary, "cos"), "sin derivative should introduce cos");
    ok &= check(contains_op(info, ad::GraphNodeKind::ScalarUnary, "tan"), "tan derivative should retain tan structurally");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::ScalarBinary) > 0, "scalar forward should use symbolic scalar binary rules");
    ok &= check(throws_cleanly([&] {
        ad::Vec bad_seed = ad::vec_variable("bad_seed", 1);
        (void)f.forward_diff(ad::vars(x, y), bad_seed);
    }), "scalar forward should reject seed size mismatch");
    return ok;
}

bool test_vector_variable_forward() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec dX = ad::vec_variable("dX", 3);
    ad::Vec tangent = X.forward_diff(ad::vars(X), dX);
    ad::GraphInfo info = ad::inspect(tangent);

    bool ok = true;
    ok &= check(tangent.size() == 3, "vector variable forward size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorVariable, "full vector variable forward should return the seed vector directly");
    ok &= check(tangent.variables().size() == 3 && tangent.variables()[0] == dX.variables()[0], "vector variable forward should depend on seed variables");
    return ok;
}

bool test_parameter_zero_tangent() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec dX = ad::vec_variable("dX", 2);
    ad::Expr h = ad::parameter("h");

    ad::Vec y = h * X;
    ad::Vec dy = y.forward_diff(ad::vars(X), dX);
    ad::GraphInfo info = ad::inspect(dy);

    bool ok = true;
    ok &= check(dy.size() == 2, "scaled parameter forward size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorScale, "h * X forward should be h * dX after zero tangent removal");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScalarParameter), "h * X forward should retain parameter h");
    ok &= check(dy.variables().size() == 2 && dy.variables()[0] == dX.variables()[0], "parameter tangent should not introduce primal variable dependency");
    ok &= check(dy.parameters().size() == 1 && dy.parameters()[0] == h.param(), "parameter should remain a runtime dependency");
    return ok;
}

bool test_scalar_dependent_scale_forward() {
    ad::Expr a = ad::variable("a");
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec seed = ad::vec_variable("seed", 3);

    ad::Vec dy = (a * X).forward_diff(ad::vars(a, X), seed);
    ad::GraphInfo info = ad::inspect(dy);

    bool ok = true;
    ok &= check(dy.size() == 2, "variable scale forward size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorBinary, "variable scale forward should add da * X and a * dX");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorScale) == 2, "variable scale forward should keep both scalar-vector scale branches");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorVariable), "variable scale forward should retain vector variable/seed dependencies");
    return ok;
}

bool test_matvec_forward() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec dX = ad::vec_variable("dX", 3);
    ad::DenseMatrix D(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});
    ad::SparseMatrix S(2, 3, {0, 1}, {0, 2}, {2.0, -1.0});

    ad::Vec dense_dy = (D * X).forward_diff(ad::vars(X), dX);
    ad::GraphInfo dense_info = ad::inspect(dense_dy);
    ad::Vec sparse_dy = (S * X).forward_diff(ad::vars(X), dX);
    ad::GraphInfo sparse_info = ad::inspect(sparse_dy);

    bool ok = true;
    ok &= check(dense_dy.size() == 2, "DenseMatVec forward size is wrong");
    ok &= check(dense_info.kind == ad::GraphNodeKind::DenseMatVec, "DenseMatVec forward should preserve DenseMatVec");
    ok &= check(dense_info.children.size() == 1 && dense_info.children[0].kind == ad::GraphNodeKind::VectorVariable, "DenseMatVec forward should use seed vector as child");
    ok &= check(ad::count_nodes(dense_info, ad::GraphNodeKind::VectorElement) == 0, "DenseMatVec forward should not scalar-lower rows");
    ok &= check(sparse_dy.size() == 2, "SparseMatVec forward size is wrong");
    ok &= check(sparse_info.kind == ad::GraphNodeKind::SparseMatVec, "SparseMatVec forward should preserve SparseMatVec");
    ok &= check(ad::count_nodes(sparse_info, ad::GraphNodeKind::VectorElement) == 0, "SparseMatVec forward should not scalar-lower rows");
    return ok;
}

bool test_elementwise_unary_forward() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec B = ad::vec_variable("B", 2);
    ad::Vec seed = ad::vec_variable("seed", 5);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});

    ad::Vec y = ad::sigmoid(A * X + B);
    ad::Vec dy = y.forward_diff(ad::vars(X, B), seed);
    ad::GraphInfo info = ad::inspect(dy);
    ad::Vec t = ad::tan(A * X + B).forward_diff(ad::vars(X, B), seed);
    ad::GraphInfo tan_info = ad::inspect(t);

    bool ok = true;
    ok &= check(dy.size() == 2, "sigmoid forward size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorUnary), "sigmoid forward should retain structured unary nodes");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorBinary), "sigmoid forward should use structured vector binary rules");
    ok &= check(contains_kind(info, ad::GraphNodeKind::DenseMatVec), "sigmoid(A * X + B) forward should retain DenseMatVec");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorElement) == 0, "sigmoid vector forward should not scalar-lower");
    ok &= check(contains_op(tan_info, ad::GraphNodeKind::VectorUnary, "tan"), "tan vector forward should retain tan structurally");
    ok &= check(ad::count_nodes(tan_info, ad::GraphNodeKind::VectorElement) == 0, "tan vector forward should not scalar-lower");
    return ok;
}

bool test_slice_and_concat_forward() {
    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec dX = ad::vec_variable("dX", 6);
    ad::Vec slice = X.slice(2, 3).forward_diff(ad::vars(X), dX);
    ad::Vec gathered = ad::gather(X, {0, 4, 2}).forward_diff(ad::vars(X), dX);

    ad::Vec A = ad::vec_variable("A", 2);
    ad::Vec B = ad::vec_variable("B", 3);
    ad::Vec seed = ad::vec_variable("seed", 5);
    ad::Vec cat = ad::concat(A, B).forward_diff(ad::vars(A, B), seed);
    ad::Vec values = ad::vec_variable("values", 3);
    ad::Vec dvalues = ad::vec_variable("dvalues", 3);
    ad::Vec scattered = ad::scatter_add(values, {0, 4, 4}, 5).forward_diff(ad::vars(values), dvalues);

    bool ok = true;
    ok &= check(slice.size() == 3 && ad::inspect(slice).kind == ad::GraphNodeKind::Slice, "slice forward should be a slice of the seed");
    ok &= check(gathered.size() == 3 && ad::inspect(gathered).kind == ad::GraphNodeKind::Gather, "gather forward should be a gather of the seed");
    ok &= check(cat.size() == 5 && ad::inspect(cat).kind == ad::GraphNodeKind::Concat, "concat forward should be concat of child tangents");
    ok &= check(scattered.size() == 5 && ad::inspect(scattered).kind == ad::GraphNodeKind::ScatterAdd, "scatter_add forward should scatter-add the seed");
    ok &= check(ad::count_nodes(ad::inspect(cat), ad::GraphNodeKind::VectorElement) == 0, "concat forward should not scalar-lower");
    ok &= check(ad::count_nodes(ad::inspect(gathered), ad::GraphNodeKind::VectorElement) == 0, "gather forward should not scalar-lower");
    ok &= check(ad::count_nodes(ad::inspect(scattered), ad::GraphNodeKind::VectorElement) == 0, "scatter_add forward should not scalar-lower");
    return ok;
}

bool test_sum_and_dot_forward() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec Y = ad::vec_variable("Y", 2);
    ad::Vec seed = ad::vec_variable("seed", 4);

    ad::Expr s = ad::sum(X).forward_diff(ad::vars(X), seed.slice(0, 2));
    ad::Expr d = ad::dot(X, Y).forward_diff(ad::vars(X, Y), seed);

    ad::GraphInfo sum_info = ad::inspect(s);
    ad::GraphInfo dot_info = ad::inspect(d);

    bool ok = true;
    ok &= check(sum_info.kind == ad::GraphNodeKind::Sum, "sum forward should remain a structured Sum reduction");
    ok &= check(ad::count_nodes(sum_info, ad::GraphNodeKind::VectorElement) == 0, "sum forward should not scalar-lower");
    ok &= check(dot_info.kind == ad::GraphNodeKind::ScalarBinary, "dot forward should combine two structured dot products");
    ok &= check(ad::count_nodes(dot_info, ad::GraphNodeKind::Dot) == 2, "dot forward should use the structured dot product rule");
    ok &= check(ad::count_nodes(dot_info, ad::GraphNodeKind::VectorElement) == 0, "dot forward should not scalar-lower");
    return ok;
}

bool test_function_forward() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));
    ad::Vec seed = ad::vec_variable("seed", rhs.input_size());

    ad::Vec jvp = rhs.forward(seed);
    ad::GraphInfo info = ad::inspect(jvp);

    bool ok = true;
    ok &= check(jvp.size() == rhs.output_size(), "function forward output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::FunctionCall, "function forward convenience should call the transformed function");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "function forward convenience should not inline transformed function body");
    ok &= check(throws_cleanly([&] {
        ad::Vec bad_seed = ad::vec_variable("bad_seed", rhs.input_size() + 1);
        (void)rhs.forward(bad_seed);
    }), "function forward should reject seed size mismatch");
    return ok;
}

bool test_function_call_forward() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Vec seed = ad::vec_variable("seed", 4);

    ad::Vec call = rhs.call({xc, uc});
    ad::Vec dcall = call.forward_diff(ad::vars(xc, uc), seed);
    ad::GraphInfo info = ad::inspect(dcall);
    ad::Vars deps = dcall.variables();

    bool ok = true;
    ok &= check(dcall.size() == 2, "function call forward output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::FunctionCall, "function call forward should create a generic FunctionCall boundary node");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "function call forward inspection should not inline callee derivative body");
    ok &= check(deps.size() == 8, "function call forward should depend on primal arguments and seed variables");
    ok &= check(deps[0] == xc.variables()[0] && deps[2] == uc.variables()[0], "function call forward should depend on outer argument variables");
    ok &= check(std::none_of(deps.values().begin(), deps.values().end(), [&](const ad::Var &var) {
        return var == x.variables()[0] || var == u.variables()[0];
    }), "callee local variables must not leak through function call forward");
    return ok;
}

bool test_local_defect_forward() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec z = ad::vec_variable("z", 4);
    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(2, 4, {1.0, 0.0, -1.0, 0.0,
                                0.0, 1.0, 0.0, -1.0});

    ad::Vec defect_expr = Dloc * z - h * rhs.call({xc, uc});
    ad::Function defect({z, xc, uc}, defect_expr);
    ad::Vec seed = ad::vec_variable("seed", defect.input_size());
    ad::Vec ddefect = defect.forward(seed);
    ad::GraphInfo info = ad::inspect(ddefect);

    bool ok = true;
    ok &= check(ddefect.size() == defect.output_size(), "local defect forward output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::FunctionCall, "local defect forward convenience should call the transformed function");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::DenseMatVec) == 0, "local defect forward convenience should not inline transformed function body");
    ok &= check(ddefect.parameters().size() == 1 && ddefect.parameters()[0] == h.param(), "local defect forward should not create parameter tangent dependencies");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_forward_basics();
    ok &= test_vector_variable_forward();
    ok &= test_parameter_zero_tangent();
    ok &= test_scalar_dependent_scale_forward();
    ok &= test_matvec_forward();
    ok &= test_elementwise_unary_forward();
    ok &= test_slice_and_concat_forward();
    ok &= test_sum_and_dot_forward();
    ok &= test_function_forward();
    ok &= test_function_call_forward();
    ok &= test_local_defect_forward();
    return ok ? 0 : 1;
}
