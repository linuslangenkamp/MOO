// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <sstream>
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

bool throws_runtime(const std::function<void()> &call) {
    try {
        call();
    } catch (const std::runtime_error &) {
        return true;
    }
    return false;
}

std::string shell_quote(const std::string &text) {
    std::string out = "'";
    for (char c : text) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out += c;
        }
    }
    out += "'";
    return out;
}

std::string number_literal(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    std::string text = out.str();
    if (text.find_first_of(".eE") == std::string::npos) {
        text += ".0";
    }
    return text;
}

std::string c_array(const std::vector<double> &values) {
    std::ostringstream out;
    out << "{";
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (i != 0) {
            out << ", ";
        }
        out << number_literal(values[i]);
    }
    out << "}";
    return out.str();
}

std::vector<double> eval_expected(const ad::Function &function,
                                  const std::vector<double> &input,
                                  const std::vector<double> &params) {
    ad::EvalWorkspace work;
    std::vector<double> out(static_cast<std::size_t>(function.output_size()), 0.0);
    function.eval(input.empty() ? nullptr : input.data(), static_cast<int>(input.size()),
                  params.empty() ? nullptr : params.data(), static_cast<int>(params.size()),
                  work,
                  out.empty() ? nullptr : out.data(), static_cast<int>(out.size()));
    return out;
}

bool compile_and_run(const std::string &name,
                     const ad::Function &function,
                     const std::vector<double> &input,
                     const std::vector<double> &params,
                     const std::string *required_substring = nullptr,
                     double tolerance = 1e-9) {
    ad::CodegenOptions options;
    options.function_name = name;
    const ad::CodegenResult generated = ad::generate_c(function, options);
    const std::vector<double> expected = eval_expected(function, input, params);
    const ad::SparsityPattern sp = function.jacobian_sparsity();

    bool ok = true;
    ok &= check(generated.entry_name == name, "generated entry name is wrong");
    ok &= check(generated.input_size == function.input_size(), "generated input size is wrong");
    ok &= check(generated.parameter_size == function.parameters().size(), "generated parameter size is wrong");
    ok &= check(generated.output_size == function.output_size(), "generated output size is wrong");
    if (required_substring) {
        ok &= check(generated.source.find(*required_substring) != std::string::npos,
                    "generated source is missing expected structural code: " + *required_substring);
    }
    if (!ok) {
        return false;
    }

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_codegen_tmp";
    std::filesystem::create_directories(dir);
    const std::filesystem::path source_path = dir / (name + ".c");
    const std::filesystem::path exe_path = dir / name;

    std::ofstream file(source_path);
    file << generated.source << "\n";
    file << "#include <math.h>\n";
    file << "#include <stdio.h>\n";
    file << "int main(void) {\n";
    file << "    if (" << name << "_input_size() != " << input.size() << ") return 10;\n";
    file << "    if (" << name << "_parameter_size() != " << params.size() << ") return 11;\n";
    file << "    if (" << name << "_output_size() != " << expected.size() << ") return 12;\n";
    file << "    if (" << name << "_jacobian_nnz() != " << sp.nnz() << ") return 13;\n";
    file << "    double in[" << input.size() << "] = " << c_array(input) << ";\n";
    if (!params.empty()) {
        file << "    double p[" << params.size() << "] = " << c_array(params) << ";\n";
    }
    file << "    double out[" << expected.size() << "] = {0};\n";
    file << "    " << name << "(in, " << (params.empty() ? "0" : "p") << ", out);\n";
    file << "    double expected[" << expected.size() << "] = " << c_array(expected) << ";\n";
    file << "    for (int i = 0; i < " << expected.size() << "; ++i) {\n";
    file << "        if (fabs(out[i] - expected[i]) > " << number_literal(tolerance) << ") {\n";
    file << "            printf(\"output mismatch at %d: %.17g vs %.17g\\n\", i, out[i], expected[i]);\n";
    file << "            return 20;\n";
    file << "        }\n";
    file << "    }\n";
    if (sp.nnz() > 0) {
        file << "    int rows[" << sp.nnz() << "];\n";
        file << "    int cols[" << sp.nnz() << "];\n";
        file << "    " << name << "_jacobian_sparsity(rows, cols);\n";
        int index = 0;
        for (const auto &entry : sp.entries()) {
            file << "    if (rows[" << index << "] != " << entry.first << " || cols[" << index << "] != " << entry.second << ") return " << (30 + index) << ";\n";
            ++index;
        }
    }
    file << "    return 0;\n";
    file << "}\n";
    file.close();

    const std::string compile_command = shell_quote(MOO_C_COMPILER) + " -std=c99 " +
                                        shell_quote(source_path.string()) + " -lm -o " +
                                        shell_quote(exe_path.string());
    if (std::system(compile_command.c_str()) != 0) {
        std::cerr << "compile command failed: " << compile_command << '\n';
        std::cerr << "generated source: " << source_path << '\n';
        return false;
    }
    const std::string run_command = shell_quote(exe_path.string());
    if (std::system(run_command.c_str()) != 0) {
        std::cerr << "generated executable failed: " << run_command << '\n';
        std::cerr << "generated source: " << source_path << '\n';
        return false;
    }
    return true;
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

bool test_simple_parameter_function() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Expr h = ad::parameter("h");
    ad::Function f({x, u}, h * ad::Vec{x[0] + u[0], x[1] * u[1]});
    return compile_and_run("simple_codegen", f, {1.0, 2.0, 3.0, 4.0}, {0.5});
}

bool test_nested_function_call_codegen() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::Vec{x[0] + u[0], x[1] * u[1]});

    ad::Vec z = ad::vec_variable("z", 2);
    ad::Vec v = ad::vec_variable("v", 2);
    ad::Function nested({z, v}, rhs.call({z, v}) + ad::vec_constant({1.0, -1.0}));
    const std::string expected_call = "nested_codegen_core_1";
    return compile_and_run("nested_codegen", nested, {2.0, 3.0, 4.0, 5.0}, {}, &expected_call);
}

bool test_mapped_output_modes_codegen() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, x);
    ad::Vec X = ad::vec_variable("X", 6);

    ad::Function concat_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }));
    ad::Function sum_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::OutputMode::Sum));
    ad::Function weighted_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::weighted_sum({1.0, 2.0, 3.0})));
    ad::Function scatter_fun({X}, ad::map(local, {
        ad::bind(x, ad::Map::blocks(X, 3, 2)),
    }, ad::MappedOutput::scatter({0, 2, 2, 3, 1, 3}, 4)));

    const std::vector<double> input = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const std::string loop = "for (int rep = 0; rep < 3; ++rep)";
    bool ok = true;
    ok &= compile_and_run("mapped_concat_codegen", concat_fun, input, {}, &loop);
    ok &= compile_and_run("mapped_sum_codegen", sum_fun, input, {}, &loop);
    ok &= compile_and_run("mapped_weighted_codegen", weighted_fun, input, {}, &loop);
    ok &= compile_and_run("mapped_scatter_codegen", scatter_fun, input, {}, &loop);
    return ok;
}

bool test_radau_collocation_codegen() {
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

    std::vector<double> input = {
        1.0, 2.0, 1.1, 2.1, 1.2, 2.2, 1.3, 2.3,
        1.1, 2.1, 1.2, 2.2, 1.3, 2.3,
        0.1, 0.2, 0.3, 0.4, 0.5, 0.6,
    };
    return compile_and_run("radau_codegen", defect, input, {0.25}, nullptr, 1e-8);
}

bool test_derivative_codegen() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, ad::Vec{ad::sin(x[0]) + x[0] * x[1]});

    ad::Function jvp = f.forward_function();
    ad::Function vjp = f.reverse_function();

    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Vec seed = ad::vec_variable("seed", f.input_size());
    ad::Vec hvp = f.reverse(lambda).forward_diff(f.input_vars(), seed);
    ad::Function hvp_fun({x, lambda, seed}, hvp);

    bool ok = true;
    ok &= compile_and_run("jvp_codegen", jvp, {0.5, 2.0, 1.0, -0.25}, {});
    ok &= compile_and_run("vjp_codegen", vjp, {0.5, 2.0, 3.0}, {});
    ok &= compile_and_run("hvp_codegen", hvp_fun, {0.5, 2.0, 3.0, 1.0, -0.25}, {});
    return ok;
}

bool test_invalid_function_name() {
    ad::Vec x = ad::vec_variable("x", 1);
    ad::Function f({x}, x);
    ad::CodegenOptions options;
    options.function_name = "not-valid";
    return check(throws_runtime([&] { (void)ad::generate_c(f, options); }),
                 "invalid C function name should throw");
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_simple_parameter_function();
    ok &= test_nested_function_call_codegen();
    ok &= test_mapped_output_modes_codegen();
    ok &= test_radau_collocation_codegen();
    ok &= test_derivative_codegen();
    ok &= test_invalid_function_name();
    return ok ? 0 : 1;
}
