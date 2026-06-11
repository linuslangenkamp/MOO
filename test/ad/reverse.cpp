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

bool test_scalar_reverse_basics() {
    ad::Expr x = ad::variable("x");
    ad::Expr y = ad::variable("y");
    ad::Expr f = ad::sin(x) + x * y;

    ad::Vec grad = f.reverse_diff(ad::vars(x, y));
    ad::GraphInfo info = ad::inspect(grad);

    bool ok = true;
    ok &= check(grad.size() == 2, "scalar reverse gradient size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScatterSlice), "scalar reverse should scatter scalar adjoints into wrt layout");
    ok &= check(contains_op(info, ad::GraphNodeKind::ScalarUnary, "cos"), "sin reverse should introduce cos");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScalarBinary), "scalar reverse should use symbolic scalar product rules");
    return ok;
}

bool test_vector_variable_reverse() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec lambda = ad::vec_variable("lambda", 3);

    ad::Vec adj = X.reverse_diff(ad::vars(X), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 3, "vector variable reverse size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorVariable, "reverse of vector variable should return lambda directly when layouts match");
    ok &= check(adj.variables().size() == 3 && adj.variables()[0] == lambda.variables()[0], "vector variable reverse should depend on lambda variables");
    return ok;
}

bool test_parameter_zero_adjoint() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Expr h = ad::parameter("h");
    ad::Vec lambda = ad::vec_variable("lambda", 2);

    ad::Vec adj = (h * X).reverse_diff(ad::vars(X), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 2, "h * X reverse size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorScale, "h * X reverse should be h * lambda");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScalarParameter), "h * X reverse should retain parameter h");
    ok &= check(adj.parameters().size() == 1 && adj.parameters()[0] == h.param(), "parameter should remain a runtime dependency");
    ok &= check(adj.variables().size() == 2 && adj.variables()[0] == lambda.variables()[0], "parameter should not create parameter adjoint output");
    return ok;
}

bool test_scalar_dependent_scale_reverse() {
    ad::Expr a = ad::variable("a");
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec lambda = ad::vec_variable("lambda", 2);

    ad::Vec adj = (a * X).reverse_diff(ad::vars(a, X), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 3, "variable scale reverse size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::Dot), "variable scale reverse should use dot(lambda, X) for scalar adjoint");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorScale), "variable scale reverse should use a * lambda for vector adjoint");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScatterSlice), "variable scale reverse should scatter scalar/vector adjoints into flattened layout");
    return ok;
}

bool test_matvec_reverse() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix D(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});
    ad::Vec dense_lambda = ad::vec_variable("dense_lambda", 2);
    ad::Vec dense_adj = (D * X).reverse_diff(ad::vars(X), dense_lambda);
    ad::GraphInfo dense_info = ad::inspect(dense_adj);

    ad::SparseMatrix S(2, 3, {0, 1}, {0, 2}, {2.0, -1.0});
    ad::Vec sparse_lambda = ad::vec_variable("sparse_lambda", 2);
    ad::Vec sparse_adj = (S * X).reverse_diff(ad::vars(X), sparse_lambda);
    ad::GraphInfo sparse_info = ad::inspect(sparse_adj);

    bool ok = true;
    ok &= check(dense_adj.size() == 3, "DenseMatVec reverse size is wrong");
    ok &= check(dense_info.kind == ad::GraphNodeKind::DenseMatVec, "DenseMatVec reverse should preserve transpose DenseMatVec");
    ok &= check(dense_info.size == 3, "transpose DenseMatVec output size is wrong");
    ok &= check(dense_info.children.size() == 1 && dense_info.children[0].kind == ad::GraphNodeKind::VectorVariable, "transpose DenseMatVec should use lambda as child");
    ok &= check(ad::count_nodes(dense_info, ad::GraphNodeKind::VectorElement) == 0, "DenseMatVec reverse should not scalar-lower rows");

    ok &= check(sparse_adj.size() == 3, "SparseMatVec reverse size is wrong");
    ok &= check(sparse_info.kind == ad::GraphNodeKind::SparseMatVec, "SparseMatVec reverse should preserve transpose SparseMatVec");
    ok &= check(sparse_info.size == 3, "transpose SparseMatVec output size is wrong");
    ok &= check(ad::count_nodes(sparse_info, ad::GraphNodeKind::VectorElement) == 0, "SparseMatVec reverse should not scalar-lower rows");
    return ok;
}

bool test_elementwise_unary_reverse() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec B = ad::vec_variable("B", 2);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});
    ad::Vec lambda = ad::vec_variable("lambda", 2);

    ad::Vec adj = ad::sigmoid(A * X + B).reverse_diff(ad::vars(X, B), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 5, "elementwise unary reverse size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorUnary), "sigmoid reverse should retain structured unary nodes");
    ok &= check(contains_kind(info, ad::GraphNodeKind::VectorBinary), "sigmoid reverse should use structured vector binary rules");
    ok &= check(contains_kind(info, ad::GraphNodeKind::DenseMatVec), "sigmoid(A * X + B) reverse should retain transpose DenseMatVec");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScatterSlice), "multi-input reverse should scatter into flattened wrt layout");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorElement) == 0, "elementwise unary reverse should not scalar-lower");
    return ok;
}

bool test_slice_and_concat_reverse() {
    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec slice_lambda = ad::vec_variable("slice_lambda", 3);
    ad::Vec slice_adj = X.slice(2, 3).reverse_diff(ad::vars(X), slice_lambda);
    ad::Vec gather_lambda = ad::vec_variable("gather_lambda", 3);
    ad::Vec gather_adj = ad::gather(X, {0, 4, 4}).reverse_diff(ad::vars(X), gather_lambda);

    ad::Vec A = ad::vec_variable("A", 2);
    ad::Vec B = ad::vec_variable("B", 3);
    ad::Vec concat_lambda = ad::vec_variable("concat_lambda", 5);
    ad::Vec concat_adj = ad::concat(A, B).reverse_diff(ad::vars(A, B), concat_lambda);
    ad::GraphInfo concat_info = ad::inspect(concat_adj);
    ad::Vec values = ad::vec_variable("values", 3);
    ad::Vec scatter_lambda = ad::vec_variable("scatter_lambda", 6);
    ad::Vec scatter_adj = ad::scatter_add(values, {0, 4, 4}, 6).reverse_diff(ad::vars(values), scatter_lambda);

    bool ok = true;
    ok &= check(slice_adj.size() == 6, "slice reverse size is wrong");
    ok &= check(ad::inspect(slice_adj).kind == ad::GraphNodeKind::ScatterSlice, "slice reverse should create ScatterSlice");
    ok &= check(gather_adj.size() == 6, "gather reverse size is wrong");
    ok &= check(ad::inspect(gather_adj).kind == ad::GraphNodeKind::ScatterAdd, "gather reverse should create ScatterAdd");
    ok &= check(ad::count_nodes(ad::inspect(gather_adj), ad::GraphNodeKind::VectorElement) == 0, "gather reverse should not scalar-lower");
    ok &= check(concat_adj.size() == 5, "concat reverse size is wrong");
    ok &= check(contains_kind(concat_info, ad::GraphNodeKind::Slice), "concat reverse should split lambda with slices");
    ok &= check(contains_kind(concat_info, ad::GraphNodeKind::ScatterSlice), "concat reverse should scatter child adjoints into wrt layout");
    ok &= check(scatter_adj.size() == 3, "scatter_add reverse size is wrong");
    ok &= check(ad::inspect(scatter_adj).kind == ad::GraphNodeKind::Gather, "scatter_add reverse should gather output adjoints");
    ok &= check(ad::count_nodes(ad::inspect(scatter_adj), ad::GraphNodeKind::VectorElement) == 0, "scatter_add reverse should not scalar-lower");
    return ok;
}

bool test_sum_and_dot_reverse() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec Y = ad::vec_variable("Y", 2);

    ad::Vec sum_adj = ad::sum(X).reverse_diff(ad::vars(X));
    ad::Vec dot_adj = ad::dot(X, Y).reverse_diff(ad::vars(X, Y));
    ad::GraphInfo sum_info = ad::inspect(sum_adj);
    ad::GraphInfo dot_info = ad::inspect(dot_adj);

    bool ok = true;
    ok &= check(sum_adj.size() == 2, "sum reverse size is wrong");
    ok &= check(sum_info.kind == ad::GraphNodeKind::VectorScale, "sum reverse should broadcast scalar adjoint structurally");
    ok &= check(ad::count_nodes(sum_info, ad::GraphNodeKind::VectorElement) == 0, "sum reverse should not scalar-lower");
    ok &= check(dot_adj.size() == 4, "dot reverse size is wrong");
    ok &= check(contains_kind(dot_info, ad::GraphNodeKind::VectorScale), "dot reverse should use structured scalar-vector product rules");
    ok &= check(contains_kind(dot_info, ad::GraphNodeKind::ScatterSlice), "dot reverse should scatter into flattened wrt layout");
    ok &= check(ad::count_nodes(dot_info, ad::GraphNodeKind::VectorElement) == 0, "dot reverse should not scalar-lower");
    return ok;
}

bool test_dot_shared_variable_reverse_accumulates() {
    ad::Vec X = ad::vec_variable("X", 2);

    ad::Vec adj = ad::dot(X, X).reverse_diff(ad::vars(X));
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 2, "shared-variable dot reverse size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorBinary, "shared-variable dot reverse should add both operand contributions");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorScale) == 2, "shared-variable dot reverse should keep two structured scale contributions");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorElement) == 0, "shared-variable dot reverse should not scalar-lower");
    return ok;
}

bool test_function_reverse() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));
    ad::Vec lambda = ad::vec_variable("lambda", rhs.output_size());

    ad::Vec vjp = rhs.reverse(lambda);
    ad::GraphInfo info = ad::inspect(vjp);

    bool ok = true;
    ok &= check(vjp.size() == rhs.input_size(), "function reverse output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::FunctionCall, "function reverse convenience should call the transformed function");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "function reverse convenience should not inline transformed function body");
    ok &= check(throws_cleanly([&] {
        ad::Vec bad_lambda = ad::vec_variable("bad_lambda", rhs.output_size() + 1);
        (void)rhs.reverse(bad_lambda);
    }), "function reverse should reject lambda size mismatch");
    return ok;
}

bool test_function_call_reverse() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Vec lambda = ad::vec_variable("lambda", 2);

    ad::Vec call = rhs.call({xc, uc});
    ad::Vec adj = call.reverse_diff(ad::vars(xc, uc), lambda);
    ad::GraphInfo info = ad::inspect(adj);
    ad::Vars deps = adj.variables();

    bool ok = true;
    ok &= check(adj.size() == 4, "function call reverse output size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "function call reverse should create a generic FunctionCall boundary node");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "function call reverse inspection should not inline callee derivative body");
    ok &= check(deps.size() == 6, "function call reverse should depend on primal arguments and lambda variables");
    ok &= check(deps[0] == lambda.variables()[0] || deps[0] == xc.variables()[0], "function call reverse should expose symbolic dependencies");
    ok &= check(std::none_of(deps.values().begin(), deps.values().end(), [&](const ad::Var &var) {
        return var == x.variables()[0] || var == u.variables()[0];
    }), "callee local variables must not leak through function call reverse");
    return ok;
}

bool test_local_defect_reverse() {
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
    ad::Vec lambda = ad::vec_variable("lambda", defect.output_size());
    ad::Vec adj = defect.reverse(lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == defect.input_size(), "local defect reverse output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::FunctionCall, "local defect reverse convenience should call the transformed function");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::DenseMatVec) == 0, "local defect reverse convenience should not inline transformed function body");
    ok &= check(adj.parameters().size() == 1 && adj.parameters()[0] == h.param(), "local defect reverse should not create parameter adjoint dependencies");
    return ok;
}

bool test_duplicate_adjoint_accumulation() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec lambda = ad::vec_variable("lambda", 2);

    ad::Vec adj = (X + X).reverse_diff(ad::vars(X), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 2, "duplicate adjoint size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::VectorBinary, "duplicate contributions should accumulate with vector addition");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorVariable) == 2, "duplicate contributions should keep both lambda branches");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_reverse_basics();
    ok &= test_vector_variable_reverse();
    ok &= test_parameter_zero_adjoint();
    ok &= test_scalar_dependent_scale_reverse();
    ok &= test_matvec_reverse();
    ok &= test_elementwise_unary_reverse();
    ok &= test_slice_and_concat_reverse();
    ok &= test_sum_and_dot_reverse();
    ok &= test_dot_shared_variable_reverse_accumulates();
    ok &= test_function_reverse();
    ok &= test_function_call_reverse();
    ok &= test_local_defect_reverse();
    ok &= test_duplicate_adjoint_accumulation();
    return ok ? 0 : 1;
}
