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

bool contains_var(const ad::Vars &vars, const ad::Var &var) {
    return std::any_of(vars.values().begin(), vars.values().end(), [&](const ad::Var &candidate) {
        return candidate == var;
    });
}

bool contains_param(const ad::Params &params, const ad::Param &param) {
    return std::any_of(params.values().begin(), params.values().end(), [&](const ad::Param &candidate) {
        return candidate == param;
    });
}

bool test_mapped_forward_block_rhs() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });
    ad::Vec seed = ad::vec_variable("seed", 12);
    ad::Vec dF = F.forward_diff(ad::vars(X, U), seed);
    ad::GraphInfo info = ad::inspect(dF);

    bool ok = true;
    ok &= check(dF.size() == F.size(), "mapped forward output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::MappedFunctionCall, "mapped forward should remain a mapped call");
    ok &= check(!contains_kind(info, ad::GraphNodeKind::VectorUnary), "mapped forward inspection should not inline local RHS body");
    ok &= check(contains_var(dF.variables(), X.variables()[0]) && contains_var(dF.variables(), U.variables()[5]), "mapped forward should depend on primal sources");
    ok &= check(contains_var(dF.variables(), seed.variables()[0]) && contains_var(dF.variables(), seed.variables()[11]), "mapped forward should depend on tangent source");
    ok &= check(!contains_var(dF.variables(), x.variables()[0]) && !contains_var(dF.variables(), u.variables()[0]), "mapped forward should not leak local variables");
    return ok;
}

bool test_mapped_forward_expression_source() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));

    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix A(4, 3, {1.0, 0.0, 0.0,
                             0.0, 2.0, 0.0,
                             0.0, 0.0, 3.0,
                             1.0, 1.0, 0.0});
    ad::Vec mapped = ad::map(local, {
        ad::bind(x, ad::Map::blocks(A * X, 2, 2)),
    });
    ad::Vec dX = ad::vec_variable("dX", 3);
    ad::Vec dmapped = mapped.forward_diff(ad::vars(X), dX);
    ad::GraphInfo info = ad::inspect(dmapped);

    bool ok = true;
    ok &= check(dmapped.size() == 4, "mapped expression-source forward size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::MappedFunctionCall, "mapped expression-source forward should remain mapped");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::DenseMatVec) >= 2, "mapped forward should preserve primal and tangent matvec sources");
    return ok;
}

bool test_mapped_forward_stencil_defect_parameter() {
    ad::Vec x = ad::vec_variable("x", 1);
    ad::Vec u = ad::vec_variable("u", 1);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec z = ad::vec_variable("z", 3);
    ad::Vec xc = ad::vec_variable("xc", 1);
    ad::Vec uc = ad::vec_variable("uc", 1);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(1, 3, {1.0, 2.0, 3.0});
    ad::Function defect({z, xc, uc}, Dloc * z - h * rhs.call({xc, uc}));

    ad::Vec X = ad::vec_variable("X", 20);
    ad::Vec U = ad::vec_variable("U", 4);
    ad::Vec defects = ad::map(defect, {
        ad::bind(z, ad::Map::stencil(X, 4, 0, 4, {0, 4, 7})),
        ad::bind(xc, ad::Map::stride(X, 4, 1, 4, 4, 1)),
        ad::bind(uc, ad::Map::blocks(U, 4, 1)),
    });
    ad::Vec seed = ad::vec_variable("seed", 24);
    ad::Vec ddefects = defects.forward_diff(ad::vars(X, U), seed);

    bool ok = true;
    ok &= check(ddefects.size() == defects.size(), "mapped stencil defect forward size is wrong");
    ok &= check(ad::inspect(ddefects).kind == ad::GraphNodeKind::MappedFunctionCall, "mapped stencil defect forward should remain mapped");
    ok &= check(contains_param(ddefects.parameters(), h.param()), "mapped stencil defect forward should retain scalar parameter dependency");
    ok &= check(!contains_var(ddefects.variables(), z.variables()[0]), "mapped stencil defect forward should not leak local stencil variable");
    return ok;
}

bool test_mapped_reverse_block_rhs() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });
    ad::Vec lambda = ad::vec_variable("lambda", F.size());
    ad::Vec adj = F.reverse_diff(ad::vars(X, U), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == 12, "mapped reverse adjoint size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::MappedFunctionCall), "mapped reverse should call a mapped local reverse kernel");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScatterAdd), "mapped reverse should scatter-add local input adjoints");
    ok &= check(contains_var(adj.variables(), lambda.variables()[0]), "mapped reverse should depend on output adjoint lambda");
    ok &= check(!contains_var(adj.variables(), x.variables()[0]) && !contains_var(adj.variables(), u.variables()[0]), "mapped reverse should not leak local variables");
    return ok;
}

bool test_mapped_reverse_duplicate_accumulation() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function identity({x}, x);

    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec F = ad::map(identity, {
        ad::bind(x, ad::Map::explicit_indices(X, 2, 2, {0, 1, 1, 2})),
    });
    ad::Vec lambda = ad::vec_variable("lambda", F.size());
    ad::Vec adj = F.reverse_diff(ad::vars(X), lambda);

    bool ok = true;
    ok &= check(adj.size() == X.size(), "duplicate mapped reverse adjoint size is wrong");
    ok &= check(contains_kind(ad::inspect(adj), ad::GraphNodeKind::ScatterAdd), "duplicate mapped reverse should use scatter_add accumulation");
    return ok;
}

bool test_mapped_reverse_expression_source() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));

    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix A(4, 3, {1.0, 0.0, 0.0,
                             0.0, 2.0, 0.0,
                             0.0, 0.0, 3.0,
                             1.0, 1.0, 0.0});
    ad::Vec mapped = ad::map(local, {
        ad::bind(x, ad::Map::blocks(A * X, 2, 2)),
    });
    ad::Vec lambda = ad::vec_variable("lambda", mapped.size());
    ad::Vec adj = mapped.reverse_diff(ad::vars(X), lambda);
    ad::GraphInfo info = ad::inspect(adj);

    bool ok = true;
    ok &= check(adj.size() == X.size(), "mapped expression-source reverse size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::DenseMatVec), "mapped expression-source reverse should propagate through transpose matvec");
    ok &= check(contains_kind(info, ad::GraphNodeKind::ScatterAdd), "mapped expression-source reverse should scatter-add mapped local adjoints");
    return ok;
}

bool test_mapped_reverse_parameter_behavior() {
    ad::Vec x = ad::vec_variable("x", 1);
    ad::Vec u = ad::vec_variable("u", 1);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec z = ad::vec_variable("z", 3);
    ad::Vec xc = ad::vec_variable("xc", 1);
    ad::Vec uc = ad::vec_variable("uc", 1);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(1, 3, {1.0, 2.0, 3.0});
    ad::Function defect({z, xc, uc}, Dloc * z - h * rhs.call({xc, uc}));

    ad::Vec X = ad::vec_variable("X", 20);
    ad::Vec U = ad::vec_variable("U", 4);
    ad::Vec defects = ad::map(defect, {
        ad::bind(z, ad::Map::stencil(X, 4, 0, 4, {0, 4, 7})),
        ad::bind(xc, ad::Map::stride(X, 4, 1, 4, 4, 1)),
        ad::bind(uc, ad::Map::blocks(U, 4, 1)),
    });
    ad::Vec lambda = ad::vec_variable("lambda", defects.size());
    ad::Vec adj = defects.reverse_diff(ad::vars(X, U), lambda);

    bool ok = true;
    ok &= check(adj.size() == 24, "mapped parameterized defect reverse size is wrong");
    ok &= check(contains_param(adj.parameters(), h.param()), "mapped parameterized defect reverse should retain scalar parameter dependency");
    ok &= check(!contains_var(adj.variables(), z.variables()[0]), "mapped parameterized defect reverse should not leak local variables");
    return ok;
}

bool test_mapped_output_mode_derivatives_construct() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));
    ad::Vec X = ad::vec_variable("X", 6);

    ad::Vec sum = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::OutputMode::Sum);
    ad::Vec weighted = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::weighted_sum({1.0, 2.0, 3.0}));
    ad::Vec scatter = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::scatter({0, 2, 2, 3, 1, 3}, 4));

    ad::Vec seed = ad::vec_variable("seed", X.size());
    ad::Vec lambda_sum = ad::vec_variable("lambda_sum", sum.size());
    ad::Vec lambda_weighted = ad::vec_variable("lambda_weighted", weighted.size());
    ad::Vec lambda_scatter = ad::vec_variable("lambda_scatter", scatter.size());

    ad::Vec dsum = sum.forward_diff(ad::vars(X), seed);
    ad::Vec dweighted = weighted.forward_diff(ad::vars(X), seed);
    ad::Vec dscatter = scatter.forward_diff(ad::vars(X), seed);
    ad::Vec adj_sum = sum.reverse_diff(ad::vars(X), lambda_sum);
    ad::Vec adj_weighted = weighted.reverse_diff(ad::vars(X), lambda_weighted);
    ad::Vec adj_scatter = scatter.reverse_diff(ad::vars(X), lambda_scatter);

    bool ok = true;
    ok &= check(dsum.size() == 2 && dweighted.size() == 2 && dscatter.size() == 4, "mapped output mode forward sizes are wrong");
    ok &= check(adj_sum.size() == X.size() && adj_weighted.size() == X.size() && adj_scatter.size() == X.size(), "mapped output mode reverse sizes are wrong");
    ok &= check(ad::inspect(dsum).kind == ad::GraphNodeKind::MappedFunctionCall, "mapped sum forward should remain mapped");
    ok &= check(contains_kind(ad::inspect(adj_sum), ad::GraphNodeKind::MappedFunctionCall), "mapped sum reverse should call mapped reverse kernel");
    ok &= check(contains_kind(ad::inspect(adj_scatter), ad::GraphNodeKind::Gather), "mapped scatter reverse should gather output adjoints");
    ok &= check(contains_kind(ad::inspect(adj_weighted), ad::GraphNodeKind::SparseMatVec), "mapped weighted reverse should preserve structured weight operator");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_mapped_forward_block_rhs();
    ok &= test_mapped_forward_expression_source();
    ok &= test_mapped_forward_stencil_defect_parameter();
    ok &= test_mapped_reverse_block_rhs();
    ok &= test_mapped_reverse_duplicate_accumulation();
    ok &= test_mapped_reverse_expression_source();
    ok &= test_mapped_reverse_parameter_behavior();
    ok &= test_mapped_output_mode_derivatives_construct();
    return ok ? 0 : 1;
}
