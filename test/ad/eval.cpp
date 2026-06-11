// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <cmath>
#include <functional>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

bool near(double lhs, double rhs, double tolerance = 1e-12) {
    return std::abs(lhs - rhs) <= tolerance;
}

bool near_array(const double *actual, const std::vector<double> &expected, double tolerance = 1e-12) {
    for (std::size_t i = 0; i < expected.size(); ++i) {
        if (!near(actual[i], expected[i], tolerance)) {
            std::cerr << "  at " << i << ": expected " << expected[i] << ", got " << actual[i] << '\n';
            return false;
        }
    }
    return true;
}

bool throws_runtime(const std::function<void()> &call) {
    try {
        call();
    } catch (const std::runtime_error &) {
        return true;
    }
    return false;
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

std::vector<double> rhs_numeric(double x0, double x1, double u0, double u1) {
    return {
        (x0 + x1 + u0 + u1) * (x0 + x1 + u0 + u1),
        x0 * u0 - x1 * std::sin(u1),
    };
}

bool test_scalar_and_id_based_values() {
    ad::Expr x1 = ad::variable("x");
    ad::Expr x2 = ad::variable("x");
    ad::Expr h = ad::parameter("x");
    ad::Expr expr = x1 - x2 + h;

    ad::Values values;
    values.set(x1, 5.0);
    values.set(x2, 2.0);
    values.set(h, 0.25);

    ad::EvalWorkspace work;
    double out = 0.0;
    expr.eval(values, work, &out);

    bool ok = true;
    ok &= check(near(out, 3.25), "scalar eval should use IDs, not labels");
    ok &= check(throws_runtime([&] {
        ad::Values missing;
        missing.set(x1, 1.0);
        expr.eval(missing, work, &out);
    }), "scalar eval should reject missing values");
    return ok;
}

bool test_vector_nodes_values_eval() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec P = ad::vec_parameter("P", 3);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, 2.0,
                             0.0, -1.0, 3.0});
    ad::SparseMatrix S(3, 3, {0, 1, 2}, {2, 0, 1}, {4.0, -2.0, 0.5});
    ad::Vec y = ad::concat(ad::concat(A * X + P.slice(0, 2),
                                      ad::gather(S * X, {1, 0})),
                           ad::scatter_add(ad::vec_constant({1.0, 2.0, 3.0}), {0, 2, 2}, 3));

    ad::Values values;
    values.set(X, std::vector<double>{1.0, 2.0, 3.0});
    values.set(P, std::vector<double>{10.0, 20.0, 30.0});

    ad::EvalWorkspace work;
    double out[7] = {};
    y.eval(values, work, out, 7);

    bool ok = true;
    ok &= check(near_array(out, {17.0, 27.0, -2.0, 12.0, 1.0, 0.0, 5.0}), "vector eval result is wrong");
    ok &= check(throws_runtime([&] {
        double short_out[6] = {};
        y.eval(values, work, short_out, 6);
    }), "vector eval should reject output size mismatch");
    return ok;
}

bool test_function_buffer_eval_and_values_eval() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Expr h = ad::parameter("h");
    ad::Function f({x, u}, h * ad::Vec{x[0] + u[0], x[1] * u[1]});

    double input[4] = {1.0, 2.0, 3.0, 4.0};
    double params[1] = {0.5};
    double out[2] = {};
    ad::EvalWorkspace work;
    work.reserve(f);
    f.eval(input, 4, params, 1, work, out, 2);

    ad::Values values;
    values.set(x, std::vector<double>{1.0, 2.0});
    values.set(u, std::vector<double>{3.0, 4.0});
    values.set(h, 0.5);
    double out_values[2] = {};
    f.eval(values, work, out_values, 2);

    bool ok = true;
    ok &= check(near_array(out, {2.0, 4.0}), "function buffer eval result is wrong");
    ok &= check(near_array(out_values, {2.0, 4.0}), "function Values eval result is wrong");
    ok &= check(throws_runtime([&] {
        f.eval(input, 3, params, 1, work, out, 2);
    }), "function eval should reject input size mismatch");
    ok &= check(throws_runtime([&] {
        f.eval(input, 4, params, 0, work, out, 2);
    }), "function eval should reject parameter size mismatch");
    return ok;
}

bool test_function_call_and_mapped_eval() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::Vec{x[0] + u[0], x[1] * u[1]});

    ad::Vec z = ad::vec_variable("z", 2);
    ad::Vec v = ad::vec_variable("v", 2);
    ad::Function nested({z, v}, rhs.call({z, v}) + ad::vec_constant({1.0, -1.0}));

    double nested_input[4] = {2.0, 3.0, 4.0, 5.0};
    double nested_out[2] = {};
    ad::EvalWorkspace work;
    nested.eval(nested_input, 4, nullptr, 0, work, nested_out, 2);

    ad::Vec X = ad::vec_variable("X", 4);
    ad::Vec U = ad::vec_variable("U", 4);
    ad::Function mapped_fun({X, U}, ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(X, 2, 2)),
        ad::bind(u, ad::Map::blocks(U, 2, 2)),
    }));
    double mapped_input[8] = {1.0, 2.0, 3.0, 4.0,
                              10.0, 20.0, 30.0, 40.0};
    double mapped_out[4] = {};
    mapped_fun.eval(mapped_input, 8, nullptr, 0, work, mapped_out, 4);

    bool ok = true;
    ok &= check(near_array(nested_out, {7.0, 14.0}), "nested function call eval result is wrong");
    ok &= check(near_array(mapped_out, {11.0, 40.0, 33.0, 160.0}), "mapped function eval result is wrong");
    return ok;
}

bool test_mapped_output_mode_eval() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, x);
    ad::Vec X = ad::vec_variable("X", 6);

    ad::Function sum_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::OutputMode::Sum));
    ad::Function weighted_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::weighted_sum({1.0, 2.0, 3.0})));
    ad::Function scatter_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::scatter({0, 2, 2, 3, 1, 3}, 4)));

    double input[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double sum_out[2] = {};
    double weighted_out[2] = {};
    double scatter_out[4] = {};
    ad::EvalWorkspace work;
    sum_fun.eval(input, 6, nullptr, 0, work, sum_out, 2);
    weighted_fun.eval(input, 6, nullptr, 0, work, weighted_out, 2);
    scatter_fun.eval(input, 6, nullptr, 0, work, scatter_out, 4);

    bool ok = true;
    ok &= check(near_array(sum_out, {9.0, 12.0}), "mapped sum output eval result is wrong");
    ok &= check(near_array(weighted_out, {22.0, 28.0}), "mapped weighted-sum output eval result is wrong");
    ok &= check(near_array(scatter_out, {1.0, 5.0, 5.0, 10.0}), "mapped scatter output eval result is wrong");
    return ok;
}

bool test_radau_collocation_eval() {
    constexpr int stages = 3;
    constexpr int nx = 2;
    constexpr int nu = 2;
    constexpr int x_size = (stages + 1) * nx;
    constexpr int u_size = stages * nu;

    ad::Vec x = ad::vec_variable("x", nx);
    ad::Vec u = ad::vec_variable("u", nu);
    ad::Function rhs({x, u}, ad::Vec{
        ad::pow(x[0] + x[1] + u[0] + u[1], 2.0),
        x[0] * u[0] - x[1] * ad::sin(u[1]),
    });

    ad::Vec z = ad::vec_variable("z", x_size);
    ad::Vec stage_x = ad::vec_variable("stage_x", stages * nx);
    ad::Vec stage_u = ad::vec_variable("stage_u", stages * nu);
    ad::Expr h = ad::parameter("h");
    ad::Vec stage_rhs = ad::map(rhs, {
        ad::bind(x, ad::Map::blocks(stage_x, stages, nx)),
        ad::bind(u, ad::Map::blocks(stage_u, stages, nu)),
    });
    ad::Function defect({z, stage_x, stage_u}, radau_stage_derivative_kron_identity(nx) * z - h * stage_rhs);

    double input[x_size + stages * nx + u_size] = {
        1.0, 2.0, 1.1, 2.1, 1.2, 2.2, 1.3, 2.3,
        1.1, 2.1, 1.2, 2.2, 1.3, 2.3,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
    };
    double params[1] = {0.25};
    double out[stages * nx] = {};
    ad::EvalWorkspace work;
    defect.eval(input, static_cast<int>(std::size(input)), params, 1, work, out, stages * nx);

    const ad::DenseMatrix D = radau_stage_derivative_kron_identity(nx);
    std::vector<double> expected(stages * nx, 0.0);
    for (int row = 0; row < stages * nx; ++row) {
        for (int col = 0; col < x_size; ++col) {
            expected[static_cast<std::size_t>(row)] += D(row, col) * input[col];
        }
    }
    for (int stage = 0; stage < stages; ++stage) {
        const double sx0 = input[x_size + stage * nx];
        const double sx1 = input[x_size + stage * nx + 1];
        const double su0 = input[x_size + stages * nx + stage * nu];
        const double su1 = input[x_size + stages * nx + stage * nu + 1];
        std::vector<double> f = rhs_numeric(sx0, sx1, su0, su1);
        expected[static_cast<std::size_t>(stage * nx)] -= params[0] * f[0];
        expected[static_cast<std::size_t>(stage * nx + 1)] -= params[0] * f[1];
    }

    return check(near_array(out, expected, 1e-10), "Radau collocation defect eval result is wrong");
}

bool test_derivative_graph_eval() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, ad::Vec{x[0] * x[0] + x[1], x[0] * x[1]});

    ad::Vec seed = ad::vec_variable("seed", f.input_size());
    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Vec jvp = f.forward(seed);
    ad::Vec adj = f.reverse(lambda);
    ad::Vec hvp = adj.forward_diff(f.input_vars(), seed);

    ad::Values values;
    values.set(x, std::vector<double>{2.0, 3.0});
    values.set(seed, std::vector<double>{5.0, 7.0});
    values.set(lambda, std::vector<double>{11.0, 13.0});

    ad::EvalWorkspace work;
    double jvp_out[2] = {};
    double adj_out[2] = {};
    double hvp_out[2] = {};
    jvp.eval(values, work, jvp_out, 2);
    adj.eval(values, work, adj_out, 2);
    hvp.eval(values, work, hvp_out, 2);

    bool ok = true;
    ok &= check(near_array(jvp_out, {27.0, 29.0}), "JVP graph eval result is wrong");
    ok &= check(near_array(adj_out, {83.0, 37.0}), "VJP graph eval result is wrong");
    ok &= check(near_array(hvp_out, {201.0, 65.0}), "HVP graph eval result is wrong");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_and_id_based_values();
    ok &= test_vector_nodes_values_eval();
    ok &= test_function_buffer_eval_and_values_eval();
    ok &= test_function_call_and_mapped_eval();
    ok &= test_mapped_output_mode_eval();
    ok &= test_radau_collocation_eval();
    ok &= test_derivative_graph_eval();
    return ok ? 0 : 1;
}
