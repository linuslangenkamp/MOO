// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
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

bool row_empty(const ad::SparsityPattern &pattern, int row) {
    return std::none_of(pattern.entries().begin(), pattern.entries().end(), [&](const std::pair<int, int> &entry) {
        return entry.first == row;
    });
}

bool test_basic_mapped_rhs_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });

    ad::SparsityPattern sp = ad::sparsity(F, ad::vars(X, U));

    bool ok = true;
    ok &= check(sp.rows() == 6 && sp.cols() == 12, "mapped RHS sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {
        {0, 0}, {0, 6},
        {1, 1}, {1, 7},
        {2, 2}, {2, 8},
        {3, 3}, {3, 9},
        {4, 4}, {4, 10},
        {5, 5}, {5, 11},
    }), "mapped RHS sparsity entries are wrong");
    return ok;
}

bool test_shifted_block_defect_sparsity() {
    ad::Vec x0 = ad::vec_variable("x0", 1);
    ad::Vec x1 = ad::vec_variable("x1", 1);
    ad::Vec u = ad::vec_variable("u", 1);
    ad::Function defect({x0, x1, u}, x1 - x0 - u);

    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec U = ad::vec_variable("U", 3);
    ad::Vec defects = ad::map(defect, {
        ad::bind(x0, ad::Map::blocks(X, 3, 1)),
        ad::bind(x1, ad::Map::shifted_blocks(X, 3, 1, 1)),
        ad::bind(u, ad::Map::blocks(U, 3, 1)),
    });

    ad::SparsityPattern sp = ad::sparsity(defects, ad::vars(X, U));

    bool ok = true;
    ok &= check(sp.rows() == 3 && sp.cols() == 7, "shifted defect sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {
        {0, 0}, {0, 1}, {0, 4},
        {1, 1}, {1, 2}, {1, 5},
        {2, 2}, {2, 3}, {2, 6},
    }), "shifted defect sparsity entries are wrong");
    return ok;
}

bool test_stencil_and_duplicate_sparsity() {
    ad::Vec z = ad::vec_variable("z", 3);
    ad::Function stencil_sum({z}, ad::Vec{ad::sum(z)});

    ad::Vec X = ad::vec_variable("X", 10);
    ad::Vec stenciled = ad::map(stencil_sum, {
        ad::bind(z, ad::Map::stencil(X, 2, 0, 1, {0, 4, 7})),
    });
    ad::Vec duplicated = ad::map(stencil_sum, {
        ad::bind(z, ad::Map::explicit_indices(X, 1, 3, {0, 4, 4})),
    });

    ad::SparsityPattern stencil_sp = ad::sparsity(stenciled, ad::vars(X));
    ad::SparsityPattern duplicate_sp = ad::sparsity(duplicated, ad::vars(X));

    bool ok = true;
    ok &= check(equals_entries(stencil_sp, {
        {0, 0}, {0, 4}, {0, 7},
        {1, 1}, {1, 5}, {1, 8},
    }), "stencil mapped sparsity entries are wrong");
    ok &= check(equals_entries(duplicate_sp, {{0, 0}, {0, 4}}), "duplicate mapped indices should be deduplicated");
    return ok;
}

bool test_expression_source_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, x);

    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix A(4, 3, {1.0, 0.0, 0.0,
                             0.0, 2.0, 0.0,
                             0.0, 0.0, 3.0,
                             1.0, 1.0, 0.0});
    ad::Vec mapped = ad::map(local, {
        ad::bind(x, ad::Map::blocks(A * X, 2, 2)),
    });

    ad::SparsityPattern sp = ad::sparsity(mapped, ad::vars(X));

    bool ok = true;
    ok &= check(sp.rows() == 4 && sp.cols() == 3, "mapped expression-source sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {
        {0, 0},
        {1, 1},
        {2, 2},
        {3, 0}, {3, 1},
    }), "mapped expression-source sparsity entries are wrong");
    return ok;
}

bool test_parameters_do_not_create_columns() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec p = ad::vec_parameter("p", 2);
    ad::Function parameterized({x}, x + p);

    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec P = ad::vec_parameter("P", 4);
    ad::Vec mapped_variable_source = ad::map(parameterized, {
        ad::bind(x, ad::Map::blocks(X, 2, 2)),
    });
    ad::Vec mapped_parameter_source = ad::map(parameterized, {
        ad::bind(x, ad::Map::blocks(P, 2, 2)),
    });

    ad::SparsityPattern variable_source = ad::sparsity(mapped_variable_source, ad::vars(X));
    ad::SparsityPattern parameter_source = ad::sparsity(mapped_parameter_source, ad::vars(X));

    bool ok = true;
    ok &= check(equals_entries(variable_source, {{0, 0}, {1, 1}, {2, 2}, {3, 3}}), "callee parameters should not add columns");
    ok &= check(parameter_source.rows() == 4 && parameter_source.cols() == 4 && parameter_source.empty(), "parameter map source should be empty w.r.t. variables");
    return ok;
}

bool test_mapped_forward_sparsity() {
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

    ad::SparsityPattern sp = ad::sparsity(dF, ad::vars(X, U, seed));

    bool ok = true;
    ok &= check(sp.rows() == 6 && sp.cols() == 24, "mapped forward sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 6) && sp.contains(0, 12) && sp.contains(0, 18), "mapped forward row 0 dependencies are wrong");
    ok &= check(sp.contains(5, 5) && sp.contains(5, 11) && sp.contains(5, 17) && sp.contains(5, 23), "mapped forward row 5 dependencies are wrong");
    return ok;
}

bool test_mapped_reverse_sparsity() {
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

    ad::SparsityPattern sp = ad::sparsity(adj, ad::vars(X, U, lambda));

    bool ok = true;
    ok &= check(sp.rows() == 12 && sp.cols() == 18, "mapped reverse sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 6) && sp.contains(0, 12), "mapped reverse X0 adjoint dependencies are wrong");
    ok &= check(sp.contains(6, 0) && sp.contains(6, 6) && sp.contains(6, 12), "mapped reverse U0 adjoint dependencies are wrong");
    ok &= check(sp.contains(5, 5) && sp.contains(5, 11) && sp.contains(5, 17), "mapped reverse X5 adjoint dependencies are wrong");
    ok &= check(sp.contains(11, 5) && sp.contains(11, 11) && sp.contains(11, 17), "mapped reverse U5 adjoint dependencies are wrong");
    return ok;
}

bool test_function_jacobian_sparsity_through_map() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));

    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec F = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 2, 2)),
    });
    ad::Function full({X}, F);

    ad::SparsityPattern sp = full.jacobian_sparsity();

    bool ok = true;
    ok &= check(sp.rows() == 4 && sp.cols() == 4, "function mapped jacobian sparsity dimensions are wrong");
    ok &= check(equals_entries(sp, {{0, 0}, {1, 1}, {2, 2}, {3, 3}}), "function mapped jacobian sparsity entries are wrong");
    return ok;
}

bool test_mapped_output_mode_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, x);
    ad::Vec X = ad::vec_variable("X", 6);

    ad::Vec sum = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::OutputMode::Sum);
    ad::Vec weighted = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::weighted_sum({1.0, 0.0, 3.0}));
    ad::Vec scatter = ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::scatter({0, 2, 2, 3, 1, 3}, 4));

    ad::SparsityPattern sum_sp = ad::sparsity(sum, ad::vars(X));
    ad::SparsityPattern weighted_sp = ad::sparsity(weighted, ad::vars(X));
    ad::SparsityPattern scatter_sp = ad::sparsity(scatter, ad::vars(X));

    bool ok = true;
    ok &= check(equals_entries(sum_sp, {
        {0, 0}, {0, 2}, {0, 4},
        {1, 1}, {1, 3}, {1, 5},
    }), "mapped sum output sparsity entries are wrong");
    ok &= check(equals_entries(weighted_sp, {
        {0, 0}, {0, 4},
        {1, 1}, {1, 5},
    }), "mapped weighted-sum output sparsity should skip zero-weight repetitions");
    ok &= check(equals_entries(scatter_sp, {
        {0, 0},
        {1, 4},
        {2, 1}, {2, 2},
        {3, 3}, {3, 5},
    }), "mapped scatter output sparsity entries are wrong");
    return ok;
}

bool test_mapped_lagrangian_hessian_vector_sparsity() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });
    ad::Function full({X, U}, F);

    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec grad = full.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(full.input_vars(), seed);

    ad::SparsityPattern sp = ad::sparsity(hvp, ad::vars(X, U, lambda, seed));

    bool ok = true;
    ok &= check(grad.size() == full.input_size(), "mapped Lagrangian gradient size is wrong");
    ok &= check(hvp.size() == full.input_size(), "mapped Lagrangian HVP size is wrong");
    ok &= check(sp.rows() == full.input_size() && sp.cols() == 30, "mapped Lagrangian HVP sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 6) && sp.contains(0, 12) && sp.contains(0, 18) && sp.contains(0, 24), "mapped Lagrangian HVP X0 dependencies are wrong");
    ok &= check(sp.contains(6, 0) && sp.contains(6, 6) && sp.contains(6, 12) && sp.contains(6, 18) && sp.contains(6, 24), "mapped Lagrangian HVP U0 dependencies are wrong");
    ok &= check(sp.contains(5, 5) && sp.contains(5, 11) && sp.contains(5, 17) && sp.contains(5, 23) && sp.contains(5, 29), "mapped Lagrangian HVP X5 dependencies are wrong");
    ok &= check(sp.contains(11, 5) && sp.contains(11, 11) && sp.contains(11, 17) && sp.contains(11, 23) && sp.contains(11, 29), "mapped Lagrangian HVP U5 dependencies are wrong");
    return ok;
}

bool test_stencil_defect_lagrangian_hessian_vector_sparsity() {
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
    ad::Function full({X, U}, defects);

    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec hvp = full.reverse(lambda).forward_diff(full.input_vars(), seed);
    ad::SparsityPattern sp = ad::sparsity(hvp, ad::vars(X, U, lambda, seed));

    bool ok = true;
    ok &= check(hvp.size() == 24, "stencil defect HVP size is wrong");
    ok &= check(sp.rows() == 24 && sp.cols() == 52, "stencil defect HVP sparsity dimensions are wrong");
    ok &= check(row_empty(sp, 0), "linear-only stencil variable X0 should not contribute to HVP sparsity");
    ok &= check(sp.contains(4, 4) && sp.contains(4, 20) && sp.contains(4, 24) && sp.contains(4, 32) && sp.contains(4, 48), "stencil defect HVP center-state row dependencies are wrong");
    ok &= check(sp.contains(20, 4) && sp.contains(20, 20) && sp.contains(20, 24) && sp.contains(20, 32) && sp.contains(20, 48), "stencil defect HVP control row dependencies are wrong");
    ok &= check(sp.contains(16, 16) && sp.contains(16, 23) && sp.contains(16, 27) && sp.contains(16, 44) && sp.contains(16, 51), "stencil defect HVP final center-state row dependencies are wrong");
    ok &= check(sp.contains(23, 16) && sp.contains(23, 23) && sp.contains(23, 27) && sp.contains(23, 44) && sp.contains(23, 51), "stencil defect HVP final control row dependencies are wrong");
    return ok;
}

bool test_overlapping_stencil_lagrangian_hessian_vector_sparsity() {
    ad::Vec z = ad::vec_variable("z", 2);
    ad::Function local({z}, ad::sigmoid(ad::Vec{z[0] + z[1]}));

    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec Y = ad::map(local, {
        ad::bind(z, ad::Map::explicit_indices(X, 2, 2, {0, 1, 1, 2})),
    });
    ad::Function full({X}, Y);

    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec hvp = full.reverse(lambda).forward_diff(full.input_vars(), seed);
    ad::SparsityPattern sp = ad::sparsity(hvp, ad::vars(X, lambda, seed));

    bool ok = true;
    ok &= check(hvp.size() == 3, "overlapping stencil HVP size is wrong");
    ok &= check(sp.rows() == 3 && sp.cols() == 8, "overlapping stencil HVP sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 1) && sp.contains(0, 3) && sp.contains(0, 5) && sp.contains(0, 6), "overlapping stencil HVP row 0 dependencies are wrong");
    ok &= check(sp.contains(1, 0) && sp.contains(1, 1) && sp.contains(1, 2) && sp.contains(1, 3) && sp.contains(1, 4) && sp.contains(1, 5) && sp.contains(1, 6) && sp.contains(1, 7), "overlapping stencil HVP shared row dependencies are wrong");
    ok &= check(sp.contains(2, 1) && sp.contains(2, 2) && sp.contains(2, 4) && sp.contains(2, 6) && sp.contains(2, 7), "overlapping stencil HVP row 2 dependencies are wrong");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_basic_mapped_rhs_sparsity();
    ok &= test_shifted_block_defect_sparsity();
    ok &= test_stencil_and_duplicate_sparsity();
    ok &= test_expression_source_sparsity();
    ok &= test_parameters_do_not_create_columns();
    ok &= test_mapped_forward_sparsity();
    ok &= test_mapped_reverse_sparsity();
    ok &= test_function_jacobian_sparsity_through_map();
    ok &= test_mapped_output_mode_sparsity();
    ok &= test_mapped_lagrangian_hessian_vector_sparsity();
    ok &= test_stencil_defect_lagrangian_hessian_vector_sparsity();
    ok &= test_overlapping_stencil_lagrangian_hessian_vector_sparsity();
    return ok ? 0 : 1;
}
