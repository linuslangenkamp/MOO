// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

bool equals_entries(const ad::SparsityPattern &pattern, std::vector<std::pair<int, int>> expected) {
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    return pattern.entries() == expected;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool test_scalar_sparsity_basics() {
    ad::Expr x = ad::variable("x");
    ad::Expr y = ad::variable("y");
    ad::Expr p = ad::parameter("p");

    ad::SparsityPattern sp = ad::sparsity(ad::sin(x) + p * y, ad::vars(x, y));

    bool ok = true;
    ok &= check(sp.rows() == 1 && sp.cols() == 2, "scalar sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {{0, 0}, {0, 1}}), "scalar sparsity entries are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 1), "scalar sparsity contains() failed");
    ok &= check(!sp.contains(0, 2), "scalar sparsity should not contain a parameter column");
    return ok;
}

bool test_same_name_variables_use_ids() {
    ad::Expr x1 = ad::variable("x");
    ad::Expr x2 = ad::variable("x");

    ad::SparsityPattern sp1 = ad::sparsity(x1, ad::vars(x2));
    ad::SparsityPattern sp2 = ad::sparsity(x1 + x2, ad::vars(x1, x2));

    bool ok = true;
    ok &= check(sp1.empty(), "same-name different-ID variable should not match wrt by label");
    ok &= check(equals_entries(sp2, {{0, 0}, {0, 1}}), "same-name variables should produce separate ID columns");
    return ok;
}

bool test_vector_variable_and_parameter() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec P = ad::vec_parameter("P", 3);

    ad::SparsityPattern x_sp = ad::sparsity(X, ad::vars(X));
    ad::SparsityPattern p_sp = ad::sparsity(P, ad::vars(X));

    bool ok = true;
    ok &= check(x_sp.rows() == 3 && x_sp.cols() == 3, "vector variable sparsity dimensions are wrong");
    ok &= check(equals_entries(x_sp, {{0, 0}, {1, 1}, {2, 2}}), "vector variable should be identity sparse pattern");
    ok &= check(p_sp.rows() == 3 && p_sp.cols() == 3 && p_sp.empty(), "vector parameter should be empty w.r.t. variables");
    return ok;
}

bool test_matvec_sparsity() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, 2.0,
                             0.0, 0.0, 3.0});
    ad::SparseMatrix S(2, 3, {0, 0, 1, 1}, {0, 0, 1, 2}, {4.0, 4.0, 0.0, 5.0});

    ad::SparsityPattern dense = ad::sparsity(A * X, ad::vars(X));
    ad::SparsityPattern sparse = ad::sparsity(S * X, ad::vars(X));

    bool ok = true;
    ok &= check(equals_entries(dense, {{0, 0}, {0, 2}, {1, 2}}), "dense matvec sparsity should respect nonzero entries only");
    ok &= check(equals_entries(sparse, {{0, 0}, {1, 2}}), "sparse matvec sparsity should ignore zero triplets and deduplicate entries");
    return ok;
}

bool test_slice_scatter_concat_sparsity() {
    ad::Vec X = ad::vec_variable("X", 5);
    ad::Vec lambda = ad::vec_variable("lambda", 2);
    ad::Vec A = ad::vec_variable("A", 2);
    ad::Vec B = ad::vec_variable("B", 3);
    ad::Vec values = ad::vec_variable("values", 3);

    ad::SparsityPattern slice = ad::sparsity(X.slice(1, 3), ad::vars(X));
    ad::SparsityPattern scatter = ad::sparsity(X.slice(2, 2).reverse_diff(ad::vars(X), lambda), ad::vars(lambda));
    ad::SparsityPattern concat = ad::sparsity(ad::concat(A, B), ad::vars(A, B));
    ad::SparsityPattern gather = ad::sparsity(ad::gather(X, {0, 4, 2}), ad::vars(X));
    ad::SparsityPattern scatter_add = ad::sparsity(ad::scatter_add(values, {0, 4, 4}, 5), ad::vars(values));

    bool ok = true;
    ok &= check(equals_entries(slice, {{0, 1}, {1, 2}, {2, 3}}), "slice sparsity row mapping is wrong");
    ok &= check(equals_entries(scatter, {{2, 0}, {3, 1}}), "scatter-slice sparsity row mapping is wrong");
    ok &= check(equals_entries(concat, {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}}), "concat sparsity row mapping is wrong");
    ok &= check(equals_entries(gather, {{0, 0}, {1, 4}, {2, 2}}), "gather sparsity row mapping is wrong");
    ok &= check(equals_entries(scatter_add, {{0, 0}, {4, 1}, {4, 2}}), "scatter_add sparsity row mapping is wrong");
    return ok;
}

bool test_sum_dot_and_scale_sparsity() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec Y = ad::vec_variable("Y", 3);
    ad::Expr a = ad::variable("a");

    ad::SparsityPattern sum_sp = ad::sparsity(ad::sum(X), ad::vars(X));
    ad::SparsityPattern dot_sp = ad::sparsity(ad::dot(X, Y), ad::vars(X, Y));
    ad::SparsityPattern dot_shared = ad::sparsity(ad::dot(X, X), ad::vars(X));
    ad::SparsityPattern scale_sp = ad::sparsity(a * X, ad::vars(a, X));

    bool ok = true;
    ok &= check(equals_entries(sum_sp, {{0, 0}, {0, 1}, {0, 2}}), "sum sparsity should depend on all vector components");
    ok &= check(equals_entries(dot_sp, {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}}), "dot sparsity should depend on both operands");
    ok &= check(equals_entries(dot_shared, {{0, 0}, {0, 1}, {0, 2}}), "dot(X, X) sparsity should deduplicate columns");
    ok &= check(equals_entries(scale_sp, {{0, 0}, {0, 1}, {1, 0}, {1, 2}, {2, 0}, {2, 3}}), "scalar-vector scale sparsity is wrong");
    return ok;
}

bool test_function_local_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function f({x, u}, ad::sigmoid(x + u));

    ad::SparsityPattern sp = f.jacobian_sparsity();

    bool ok = true;
    ok &= check(sp.rows() == 2 && sp.cols() == 4, "function sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {{0, 0}, {0, 2}, {1, 1}, {1, 3}}), "function local sparsity entries are wrong");
    return ok;
}

bool test_function_call_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Vec y = rhs.call({xc, uc});

    ad::SparsityPattern sp = ad::sparsity(y, ad::vars(xc, uc));

    bool ok = true;
    ok &= check(sp.rows() == 2 && sp.cols() == 4, "function call sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {{0, 0}, {0, 2}, {1, 1}, {1, 3}}), "function call sparsity should compose callee sparsity through arguments");
    return ok;
}

bool test_function_call_argument_expression_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec U = ad::vec_variable("U", 2);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, 2.0,
                             0.0, 3.0, 0.0});

    ad::SparsityPattern sp = ad::sparsity(rhs.call({A * X, U}), ad::vars(X, U));

    bool ok = true;
    ok &= check(sp.rows() == 2 && sp.cols() == 5, "function call expression sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {{0, 0}, {0, 2}, {0, 3}, {1, 1}, {1, 4}}), "function call should compose through matvec argument sparsity");
    return ok;
}

bool test_function_forward_transformed_call_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Vec seed = ad::vec_variable("seed", 4);
    ad::Vec jvp = rhs.call({xc, uc}).forward_diff(ad::vars(xc, uc), seed);

    ad::SparsityPattern sp = ad::sparsity(jvp, ad::vars(xc, uc, seed));

    bool ok = true;
    ok &= check(ad::inspect(jvp).kind == ad::GraphNodeKind::FunctionCall, "forward through call should build generic FunctionCall");
    ok &= check(sp.rows() == 2 && sp.cols() == 8, "transformed forward FunctionCall sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2) && sp.contains(0, 4) && sp.contains(0, 6), "transformed forward call row 0 should include primal and tangent dependencies");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 3) && sp.contains(1, 5) && sp.contains(1, 7), "transformed forward call row 1 should include primal and tangent dependencies");
    return ok;
}

bool test_function_reverse_transformed_call_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Vec lambda = ad::vec_variable("lambda", 2);
    ad::Vec vjp = rhs.call({xc, uc}).reverse_diff(ad::vars(xc, uc), lambda);

    ad::SparsityPattern sp = ad::sparsity(vjp, ad::vars(xc, uc, lambda));

    bool ok = true;
    ok &= check(contains_kind(ad::inspect(vjp), ad::GraphNodeKind::FunctionCall), "reverse through call should build generic FunctionCall");
    ok &= check(sp.rows() == 4 && sp.cols() == 6, "transformed reverse FunctionCall sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2) && sp.contains(0, 4), "transformed reverse call x0 adjoint dependencies are wrong");
    ok &= check(sp.contains(2, 0) && sp.contains(2, 2) && sp.contains(2, 4), "transformed reverse call u0 adjoint dependencies are wrong");
    return ok;
}

bool test_hvp_transformed_call_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function rhs({x}, ad::sigmoid(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Function outer({y}, rhs.call({y}));
    ad::Vec lambda = ad::vec_variable("lambda", outer.output_size());
    ad::Vec seed = ad::vec_variable("seed", outer.input_size());
    ad::Vec hvp = outer.reverse(lambda).forward_diff(outer.input_vars(), seed);

    ad::SparsityPattern sp = ad::sparsity(hvp, ad::vars(y, lambda, seed));

    bool ok = true;
    ok &= check(contains_kind(ad::inspect(hvp), ad::GraphNodeKind::FunctionCall), "HVP composition should build generic transformed FunctionCall");
    ok &= check(sp.rows() == 2 && sp.cols() == 6, "HVP transformed FunctionCall sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2) && sp.contains(0, 4), "HVP transformed call row 0 dependencies are wrong");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 3) && sp.contains(1, 5), "HVP transformed call row 1 dependencies are wrong");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_sparsity_basics();
    ok &= test_same_name_variables_use_ids();
    ok &= test_vector_variable_and_parameter();
    ok &= test_matvec_sparsity();
    ok &= test_slice_scatter_concat_sparsity();
    ok &= test_sum_dot_and_scale_sparsity();
    ok &= test_function_local_sparsity();
    ok &= test_function_call_sparsity();
    ok &= test_function_call_argument_expression_sparsity();
    ok &= test_function_forward_transformed_call_sparsity();
    ok &= test_function_reverse_transformed_call_sparsity();
    ok &= test_hvp_transformed_call_sparsity();
    return ok ? 0 : 1;
}
