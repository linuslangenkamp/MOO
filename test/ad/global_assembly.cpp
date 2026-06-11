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

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool contains_param(const ad::Params &params, const ad::Param &param) {
    return std::any_of(params.values().begin(), params.values().end(), [&](const ad::Param &candidate) {
        return candidate == param;
    });
}

bool equals_entries(const ad::SparsityPattern &pattern, std::vector<std::pair<int, int>> expected) {
    std::sort(expected.begin(), expected.end());
    expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
    return pattern.entries() == expected;
}

ad::DenseMatrix diagonal_dense(int size) {
    std::vector<double> values(static_cast<std::size_t>(size * size), 0.0);
    for (int i = 0; i < size; ++i) {
        values[static_cast<std::size_t>(i * size + i)] = 1.0;
    }
    return ad::DenseMatrix(size, size, std::move(values));
}

ad::DenseMatrix radau_stage_derivative_kron_identity(int state_size) {
    constexpr int stages = 3;
    constexpr int nodes = stages + 1;
    const double radau[nodes * nodes] = {
        -9.0, 10.048809399827416, -1.382142733160749, 0.3333333333333333,
        -4.139387691339814, 3.224744871391589, 1.1678400846904055, -0.2531972647421808,
        1.7393876913398137, -3.5678400846904055, 0.775255128608411, 1.0531972647421808,
        -3.0, 5.531972647421808, -7.531972647421809, 5.0,
    };
    const int rows = stages * state_size;
    const int cols = nodes * state_size;
    std::vector<double> values(static_cast<std::size_t>(rows * cols), 0.0);
    for (int stage = 0; stage < stages; ++stage) {
        for (int state = 0; state < state_size; ++state) {
            const int row = stage * state_size + state;
            for (int node = 0; node < nodes; ++node) {
                const int col = node * state_size + state;
                values[static_cast<std::size_t>(row * cols + col)] = radau[(stage + 1) * nodes + node];
            }
        }
    }
    return ad::DenseMatrix(rows, cols, std::move(values));
}

bool test_rhs_map_global_collocation_residual() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec X = ad::vec_variable("X", 6);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Expr h = ad::parameter("h");

    ad::Vec F = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });
    ad::Vec defects = diagonal_dense(6) * X - h * F;
    ad::Vec ic = X.slice(0, 2) - ad::vec_constant(std::vector<double>{1.0, 2.0});
    ad::Vec fc = X.slice(4, 2);
    ad::Vec residual = ad::concat(ad::concat(ad::concat(ic, defects), fc), ad::Vec{ad::dot(F, F)});
    ad::Function full({X, U}, residual);

    ad::Function df = full.forward_function();
    ad::Function rf = full.reverse_function();
    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec jvp = full.forward(seed);
    ad::Vec adj = full.reverse(lambda);
    ad::Vec hvp = adj.forward_diff(full.input_vars(), seed);
    ad::SparsityPattern jac = full.jacobian_sparsity();
    ad::SparsityPattern hvp_sp = ad::sparsity(hvp, ad::vars(X, U, lambda, seed));

    bool ok = true;
    ok &= check(full.output_size() == 11 && full.input_size() == 12, "global RHS-map residual layout is wrong");
    ok &= check(contains_param(full.parameters(), h.param()), "global RHS-map residual should collect h parameter");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "global RHS-map residual should keep mapped RHS boundary");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::DenseMatVec), "global RHS-map residual should keep global DenseMatVec");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::Slice), "global RHS-map residual should keep boundary slices");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::Dot), "global RHS-map residual should keep objective-style dot reduction");
    ok &= check(!contains_kind(full.info().output_graph, ad::GraphNodeKind::VectorUnary), "global mapped RHS body should not be inlined");
    ok &= check(contains_kind(df.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "global RHS-map forward function should keep mapped derivative boundary");
    ok &= check(contains_kind(rf.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "global RHS-map reverse function should keep mapped derivative boundary");
    ok &= check(jvp.size() == full.output_size() && adj.size() == full.input_size() && hvp.size() == full.input_size(), "global RHS-map derivative sizes are wrong");
    ok &= check(jac.rows() == 11 && jac.cols() == 12, "global RHS-map jacobian sparsity dimensions are wrong");
    ok &= check(jac.contains(0, 0) && jac.contains(1, 1), "global RHS-map initial condition sparsity is wrong");
    ok &= check(jac.contains(2, 0) && jac.contains(2, 6), "global RHS-map first defect sparsity is wrong");
    ok &= check(jac.contains(7, 5) && jac.contains(7, 11), "global RHS-map last defect sparsity is wrong");
    ok &= check(jac.contains(8, 4) && jac.contains(9, 5), "global RHS-map final condition sparsity is wrong");
    ok &= check(jac.contains(10, 0) && jac.contains(10, 5) && jac.contains(10, 6) && jac.contains(10, 11), "global RHS-map reduction sparsity is wrong");
    ok &= check(hvp_sp.rows() == full.input_size() && hvp_sp.cols() == 35, "global RHS-map HVP sparsity dimensions are wrong");
    ok &= check(hvp_sp.contains(0, 0) && hvp_sp.contains(0, 6) && hvp_sp.contains(0, 14) && hvp_sp.contains(0, 22) && hvp_sp.contains(0, 23), "global RHS-map HVP representative dependencies are wrong");
    return ok;
}

bool test_mapped_local_stencil_defect_global_residual() {
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
    ad::Vec residual = ad::concat(defects, X.slice(0, 1));
    ad::Function full({X, U}, residual);

    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec hvp = full.reverse(lambda).forward_diff(full.input_vars(), seed);
    ad::SparsityPattern jac = full.jacobian_sparsity();
    ad::SparsityPattern hvp_sp = ad::sparsity(hvp, ad::vars(X, U, lambda, seed));

    bool ok = true;
    ok &= check(full.output_size() == 5 && full.input_size() == 24, "mapped local stencil residual layout is wrong");
    ok &= check(contains_param(full.parameters(), h.param()), "mapped local stencil residual should collect h parameter");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "mapped local stencil residual should keep mapped boundary");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::Slice), "mapped local stencil residual should keep boundary slice");
    ok &= check(!contains_kind(full.info().output_graph, ad::GraphNodeKind::FunctionCall), "mapped inspection should not inline nested RHS call from defect kernel");
    ok &= check(jac.rows() == 5 && jac.cols() == 24, "mapped local stencil jacobian dimensions are wrong");
    ok &= check(jac.contains(0, 0) && jac.contains(0, 4) && jac.contains(0, 7) && jac.contains(0, 20), "mapped local stencil first row sparsity is wrong");
    ok &= check(jac.contains(3, 12) && jac.contains(3, 16) && jac.contains(3, 19) && jac.contains(3, 23), "mapped local stencil last defect row sparsity is wrong");
    ok &= check(jac.contains(4, 0), "mapped local stencil boundary sparsity is wrong");
    ok &= check(hvp.size() == full.input_size(), "mapped local stencil HVP size is wrong");
    ok &= check(hvp_sp.rows() == full.input_size() && hvp_sp.cols() == 53, "mapped local stencil HVP sparsity dimensions are wrong");
    ok &= check(hvp_sp.contains(4, 4) && hvp_sp.contains(4, 20) && hvp_sp.contains(4, 24) && hvp_sp.contains(4, 33) && hvp_sp.contains(4, 49), "mapped local stencil HVP representative dependencies are wrong");
    return ok;
}

bool test_current_next_multiple_shooting_residual() {
    ad::Vec x0 = ad::vec_variable("x0", 2);
    ad::Vec x1 = ad::vec_variable("x1", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Expr h = ad::parameter("h");
    ad::Function defect({x0, x1, u}, x1 - x0 - h * ad::sigmoid(x0 + u));

    ad::Vec X = ad::vec_variable("X", 8);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec defects = ad::map(defect, {
        ad::bind(x0, ad::Map::blocks(X, 3, 2)),
        ad::bind(x1, ad::Map::shifted_blocks(X, 3, 2, 1)),
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });
    ad::Vec residual = ad::concat(X.slice(0, 2), defects);
    ad::Function full({X, U}, residual);

    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec adj = full.reverse(lambda);
    ad::Vec hvp = adj.forward_diff(full.input_vars(), seed);
    ad::SparsityPattern jac = full.jacobian_sparsity();
    ad::SparsityPattern adj_sp = ad::sparsity(adj, ad::vars(X, U, lambda));

    bool ok = true;
    ok &= check(full.output_size() == 8 && full.input_size() == 14, "multiple-shooting residual layout is wrong");
    ok &= check(contains_param(full.parameters(), h.param()), "multiple-shooting residual should collect h parameter");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "multiple-shooting residual should keep mapped boundary");
    ok &= check(jac.rows() == 8 && jac.cols() == 14, "multiple-shooting jacobian dimensions are wrong");
    ok &= check(jac.contains(0, 0) && jac.contains(1, 1), "multiple-shooting boundary sparsity is wrong");
    ok &= check(jac.contains(2, 0) && jac.contains(2, 2) && jac.contains(2, 8), "multiple-shooting first defect sparsity is wrong");
    ok &= check(jac.contains(6, 4) && jac.contains(6, 6) && jac.contains(6, 12), "multiple-shooting last defect sparsity is wrong");
    ok &= check(adj.size() == full.input_size() && hvp.size() == full.input_size(), "multiple-shooting derivative sizes are wrong");
    ok &= check(adj_sp.rows() == 14 && adj_sp.cols() == 22, "multiple-shooting reverse sparsity dimensions are wrong");
    ok &= check(adj_sp.contains(2, 2) && adj_sp.contains(2, 10) && adj_sp.contains(2, 16) && adj_sp.contains(2, 18), "multiple-shooting shared-state reverse sparsity is wrong");
    return ok;
}

bool test_expression_source_global_assembly() {
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
    ad::Vec residual = ad::concat(mapped, X.slice(0, 1));
    ad::Function full({X}, residual);

    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Vec jvp = full.forward(seed);
    ad::Vec adj = full.reverse(lambda);
    ad::SparsityPattern jac = full.jacobian_sparsity();

    bool ok = true;
    ok &= check(full.output_size() == 5 && full.input_size() == 3, "expression-source global residual layout is wrong");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "expression-source global residual should keep mapped boundary");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::DenseMatVec), "expression-source global residual should keep source matvec");
    ok &= check(jvp.size() == full.output_size() && adj.size() == full.input_size(), "expression-source global derivative sizes are wrong");
    ok &= check(equals_entries(jac, {
        {0, 0},
        {1, 1},
        {2, 2},
        {3, 0}, {3, 1},
        {4, 0},
    }), "expression-source global jacobian sparsity entries are wrong");
    return ok;
}

bool test_local_rhs_lagrangian_hessian_sparsity_is_sparse_blocks() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Expr rhs0 = x[0] * x[0] + x[0] * u[0];
    ad::Expr rhs1 = x[1] * x[1] + x[1] * u[1];
    ad::Function rhs({x, u}, ad::Vec{rhs0, rhs1});

    ad::Vec lambda = ad::vec_variable("lambda", rhs.output_size());
    ad::Vec seed = ad::vec_variable("seed", rhs.input_size());
    ad::Vec grad = rhs.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(rhs.input_vars(), seed);
    ad::SparsityPattern hvp_sp = ad::sparsity(hvp, ad::vars(x, u, lambda, seed));

    const int seed_col_offset = rhs.input_size() + rhs.output_size();
    std::vector<std::pair<int, int>> actual_seed_block;
    for (const auto &entry : hvp_sp.entries()) {
        if (entry.second >= seed_col_offset && entry.second < seed_col_offset + rhs.input_size()) {
            actual_seed_block.emplace_back(entry.first, entry.second - seed_col_offset);
        }
    }

    bool ok = true;
    ok &= check(hvp.size() == 4, "local RHS Lagrangian HVP size is wrong");
    ok &= check(hvp_sp.rows() == 4 && hvp_sp.cols() == 10, "local RHS Lagrangian HVP sparsity dimensions are wrong");
    ok &= check(actual_seed_block == std::vector<std::pair<int, int>>{
                      {0, 0}, {0, 2},
                      {1, 1}, {1, 3},
                      {2, 0},
                      {3, 1},
                  },
                  "local RHS Lagrangian Hessian seed block is wrong");
    return ok;
}

bool test_stage_interval_collocation_defect_shape() {
    constexpr int intervals = 3;
    constexpr int stages = 3;
    constexpr int nx = 2;
    constexpr int nu = 2;
    constexpr int state_nodes = intervals * stages + 1;
    constexpr int control_nodes = intervals * stages;
    constexpr int x_size = state_nodes * nx;
    constexpr int u_size = control_nodes * nu;
    constexpr int stage_state_size = stages * nx;
    constexpr int stage_control_size = stages * nu;
    constexpr int stencil_state_size = (stages + 1) * nx;

    ad::Vec x = ad::vec_variable("x", nx);
    ad::Vec u = ad::vec_variable("u", nu);
    ad::Expr rhs0 = ad::pow(x[0] + x[1] + u[0] + u[1], 2.0);
    ad::Expr rhs1 = x[0] * u[0] - x[1] * ad::sin(u[1]);
    ad::Function rhs({x, u}, ad::Vec{rhs0, rhs1});

    ad::Vec z = ad::vec_variable("z", stencil_state_size);
    ad::Vec stage_x = ad::vec_variable("stage_x", stage_state_size);
    ad::Vec stage_u = ad::vec_variable("stage_u", stage_control_size);
    ad::Expr h = ad::parameter("h");
    ad::Vec stage_rhs = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(stage_x, stages, nx)),
        ad::bind(u, ad::Map::blocks(stage_u, stages, nu)),
    });
    ad::Function interval_defect({z, stage_x, stage_u},
                                 radau_stage_derivative_kron_identity(nx) * z - h * stage_rhs);

    ad::Vec X = ad::vec_variable("X", x_size);
    ad::Vec U = ad::vec_variable("U", u_size);
    ad::Vec defects = ad::map(interval_defect, {
        ad::bind(z, ad::Map::stride(X, intervals, stencil_state_size, 0, stages * nx, 1)),
        ad::bind(stage_x, ad::Map::stride(X, intervals, stage_state_size, nx, stages * nx, 1)),
        ad::bind(stage_u, ad::Map::blocks(U, intervals, stage_control_size)),
    });
    ad::Function full({X, U}, defects);

    ad::Vec seed = ad::vec_variable("seed", full.input_size());
    ad::Vec lambda = ad::vec_variable("lambda", full.output_size());
    ad::Function df = full.forward_function();
    ad::Function rf = full.reverse_function();
    ad::Vec jvp = full.forward(seed);
    ad::Vec adj = full.reverse(lambda);
    ad::Vec hvp = adj.forward_diff(full.input_vars(), seed);
    ad::SparsityPattern jac = full.jacobian_sparsity();
    ad::SparsityPattern hvp_sp = ad::sparsity(hvp, ad::vars(X, U, lambda, seed));

    bool ok = true;
    ok &= check(full.output_size() == intervals * stage_state_size && full.input_size() == x_size + u_size, "stage-interval collocation layout is wrong");
    ok &= check(contains_param(full.parameters(), h.param()), "stage-interval collocation should collect h parameter");
    ok &= check(contains_kind(full.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "stage-interval collocation should keep interval mapped boundary");
    ok &= check(contains_kind(interval_defect.info().output_graph, ad::GraphNodeKind::DenseMatVec), "stage-interval local defect should keep D kron I matvec");
    ok &= check(contains_kind(interval_defect.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "stage-interval local defect should keep mapped stage RHS boundary");
    ok &= check(!contains_kind(full.info().output_graph, ad::GraphNodeKind::FunctionCall), "stage-interval global inspection should not inline RHS calls");
    ok &= check(contains_kind(df.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "stage-interval forward function should keep mapped boundary");
    ok &= check(contains_kind(rf.info().output_graph, ad::GraphNodeKind::MappedFunctionCall), "stage-interval reverse function should keep mapped boundary");
    ok &= check(jvp.size() == full.output_size() && adj.size() == full.input_size() && hvp.size() == full.input_size(), "stage-interval derivative sizes are wrong");
    ok &= check(jac.rows() == intervals * stage_state_size && jac.cols() == x_size + u_size, "stage-interval jacobian dimensions are wrong");
    ok &= check(jac.contains(0, 0) && jac.contains(0, 2) && jac.contains(0, 4) && jac.contains(0, 6) && jac.contains(0, x_size) && jac.contains(0, x_size + 1), "stage-interval first stage/state sparsity is wrong");
    ok &= check(jac.contains(1, 1) && jac.contains(1, 2) && jac.contains(1, 3) && jac.contains(1, 5) && jac.contains(1, 7) && jac.contains(1, x_size) && jac.contains(1, x_size + 1), "stage-interval first stage second-state sparsity is wrong");
    ok &= check(jac.contains(17, 13) && jac.contains(17, 15) && jac.contains(17, 17) && jac.contains(17, 18) && jac.contains(17, 19) && jac.contains(17, x_size + 16) && jac.contains(17, x_size + 17), "stage-interval last stage sparsity is wrong");
    ok &= check(hvp_sp.rows() == full.input_size() && hvp_sp.cols() == full.input_size() + full.output_size() + full.input_size(), "stage-interval HVP sparsity dimensions are wrong");
    ok &= check(hvp_sp.contains(2, 2) && hvp_sp.contains(2, 20) && hvp_sp.contains(2, 38) && hvp_sp.contains(2, 58) && hvp_sp.contains(2, 76), "stage-interval HVP representative dependencies are wrong");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_rhs_map_global_collocation_residual();
    ok &= test_mapped_local_stencil_defect_global_residual();
    ok &= test_current_next_multiple_shooting_residual();
    ok &= test_expression_source_global_assembly();
    ok &= test_local_rhs_lagrangian_hessian_sparsity_is_sparse_blocks();
    ok &= test_stage_interval_collocation_defect_shape();
    return ok ? 0 : 1;
}
