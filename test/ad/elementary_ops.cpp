// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
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

bool near(double lhs, double rhs, double tolerance = 1e-10) {
    return std::abs(lhs - rhs) <= tolerance;
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

std::vector<double> eval_function(const ad::Function &function, const std::vector<double> &input) {
    ad::EvalWorkspace work;
    std::vector<double> out(static_cast<std::size_t>(function.output_size()), 0.0);
    function.eval(input.data(), static_cast<int>(input.size()), nullptr, 0, work, out.data(), static_cast<int>(out.size()));
    return out;
}

bool compile_and_run(const std::string &name,
                     const ad::Function &function,
                     const std::vector<double> &input,
                     double tolerance = 1e-9) {
    ad::CodegenOptions options;
    options.function_name = name;
    const ad::CodegenResult generated = ad::generate_c(function, options);
    const std::vector<double> expected = eval_function(function, input);

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_elementary_ops_tmp";
    std::filesystem::create_directories(dir);
    const std::filesystem::path source_path = dir / (name + ".c");
    const std::filesystem::path exe_path = dir / name;

    std::ofstream file(source_path);
    file << generated.source << "\n";
    file << "#include <math.h>\n";
    file << "#include <stdio.h>\n";
    file << "int main(void) {\n";
    file << "    double in[" << input.size() << "] = " << c_array(input) << ";\n";
    file << "    double out[" << expected.size() << "] = {0};\n";
    file << "    " << name << "(in, 0, out);\n";
    file << "    double expected[" << expected.size() << "] = " << c_array(expected) << ";\n";
    file << "    for (int i = 0; i < " << expected.size() << "; ++i) {\n";
    file << "        if (fabs(out[i] - expected[i]) > " << number_literal(tolerance) << ") {\n";
    file << "            printf(\"output mismatch at %d: %.17g vs %.17g\\n\", i, out[i], expected[i]);\n";
    file << "            return 20;\n";
    file << "        }\n";
    file << "    }\n";
    file << "    return 0;\n";
    file << "}\n";
    file.close();

    const std::string compile_command = shell_quote(MOO_C_COMPILER) + " -std=c99 " +
                                        shell_quote(source_path.string()) + " -lm -o " +
                                        shell_quote(exe_path.string());
    if (std::system(compile_command.c_str()) != 0) {
        std::cerr << "compile command failed: " << compile_command << '\n';
        return false;
    }
    if (std::system(shell_quote(exe_path.string()).c_str()) != 0) {
        std::cerr << "generated executable failed: " << exe_path << '\n';
        return false;
    }
    return true;
}

bool test_scalar_primal_eval() {
    ad::Vec input = ad::vec_variable("input", 2);
    ad::Expr x = input[0];
    ad::Expr y = input[1];
    ad::Vec out{
        ad::abs(x),
        ad::sqrt(y),
        ad::asin(x / 2.0),
        ad::acos(x / 2.0),
        ad::atan(x),
        ad::sinh(x),
        ad::cosh(x),
        ad::tanh(x),
        ad::log10(y),
        ad::sigmoid(x),
        ad::pow(y, x),
        ad::min(x, y),
        ad::max(x, y),
    };
    ad::Function f({input}, out);
    const std::vector<double> actual = eval_function(f, {-0.5, 2.0});

    bool ok = true;
    ok &= check(near(actual[0], 0.5), "scalar abs eval is wrong");
    ok &= check(near(actual[1], std::sqrt(2.0)), "scalar sqrt eval is wrong");
    ok &= check(near(actual[2], std::asin(-0.25)), "scalar asin eval is wrong");
    ok &= check(near(actual[3], std::acos(-0.25)), "scalar acos eval is wrong");
    ok &= check(near(actual[4], std::atan(-0.5)), "scalar atan eval is wrong");
    ok &= check(near(actual[5], std::sinh(-0.5)), "scalar sinh eval is wrong");
    ok &= check(near(actual[6], std::cosh(-0.5)), "scalar cosh eval is wrong");
    ok &= check(near(actual[7], std::tanh(-0.5)), "scalar tanh eval is wrong");
    ok &= check(near(actual[8], std::log10(2.0)), "scalar log10 eval is wrong");
    ok &= check(near(actual[9], 1.0 / (1.0 + std::exp(0.5))), "scalar sigmoid eval is wrong");
    ok &= check(near(actual[10], std::pow(2.0, -0.5)), "scalar pow(expr, expr) eval is wrong");
    ok &= check(near(actual[11], -0.5), "scalar min eval is wrong");
    ok &= check(near(actual[12], 2.0), "scalar max eval is wrong");
    return ok;
}

bool test_vector_primal_eval() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Vec Y = ad::vec_variable("Y", 3);
    ad::Function f({X, Y}, ad::concat(ad::concat(-X,
                                                 ad::sqrt(Y) + ad::tanh(X) + ad::pow(Y, X)),
                                      ad::concat(ad::min(X, Y), ad::max(X, Y))));
    const std::vector<double> actual = eval_function(f, {-1.0, 0.5, 2.0, 4.0, 9.0, 16.0});

    bool ok = true;
    ok &= check(f.output_size() == 12, "vector elementary output size is wrong");
    ok &= check(near(actual[0], 1.0) && near(actual[1], -0.5) && near(actual[2], -2.0), "vector neg eval is wrong");
    ok &= check(near(actual[3], std::sqrt(4.0) + std::tanh(-1.0) + std::pow(4.0, -1.0)), "vector expression eval row 0 is wrong");
    ok &= check(near(actual[6], -1.0) && near(actual[8], 2.0), "vector min eval is wrong");
    ok &= check(near(actual[9], 4.0) && near(actual[11], 16.0), "vector max eval is wrong");
    return ok;
}

bool test_smooth_forward_reverse_eval() {
    ad::Expr x = ad::variable("x");
    ad::Expr y = ad::variable("y");
    ad::Expr f = ad::sqrt(y) + ad::log10(y) + ad::tanh(x) + ad::pow(y, x) + ad::sigmoid(x);
    ad::Vec seed = ad::vec_variable("seed", 2);
    ad::Expr df = f.forward_diff(ad::vars(x, y), seed);
    ad::Vec grad = f.reverse_diff(ad::vars(x, y));

    ad::Values values;
    values.set(x, 0.25);
    values.set(y, 3.0);
    values.set(seed, std::vector<double>{1.5, -0.25});
    ad::EvalWorkspace work;
    double df_value = 0.0;
    df.eval(values, work, &df_value);
    double grad_value[2] = {};
    grad.eval(values, work, grad_value, 2);

    const double sig = 1.0 / (1.0 + std::exp(-0.25));
    const double dfdx = (1.0 - std::pow(std::tanh(0.25), 2.0)) + std::pow(3.0, 0.25) * std::log(3.0) + sig * (1.0 - sig);
    const double dfdy = 1.0 / (2.0 * std::sqrt(3.0)) + 1.0 / (3.0 * std::log(10.0)) + std::pow(3.0, 0.25) * 0.25 / 3.0;
    const double expected_jvp = dfdx * 1.5 + dfdy * -0.25;

    bool ok = true;
    ok &= check(near(df_value, expected_jvp), "smooth scalar forward eval is wrong");
    ok &= check(near(grad_value[0], dfdx), "smooth scalar reverse x adjoint is wrong");
    ok &= check(near(grad_value[1], dfdy), "smooth scalar reverse y adjoint is wrong");
    return ok;
}

bool test_nonsmooth_derivatives_reject() {
    ad::Expr x = ad::variable("x");
    ad::Expr y = ad::variable("y");
    ad::Vec seed = ad::vec_variable("seed", 2);

    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec Y = ad::vec_variable("Y", 2);
    ad::Vec dX = ad::vec_variable("dX", 2);
    ad::Vec lambda = ad::vec_variable("lambda", 2);

    bool ok = true;
    ok &= check(throws_runtime([&] {
        (void)ad::abs(x).forward_diff(ad::vars(x, y), seed);
    }), "forward derivative of scalar abs should reject");
    ok &= check(throws_runtime([&] {
        (void)ad::abs(x).reverse_diff(ad::vars(x));
    }), "reverse derivative of scalar abs should reject");
    ok &= check(throws_runtime([&] {
        (void)ad::min(x, y).forward_diff(ad::vars(x, y), seed);
    }), "forward derivative of scalar min should reject");
    ok &= check(throws_runtime([&] {
        (void)ad::max(x, y).reverse_diff(ad::vars(x, y));
    }), "reverse derivative of scalar max should reject");
    ok &= check(throws_runtime([&] {
        (void)ad::abs(X).forward_diff(ad::vars(X), dX);
    }), "forward derivative of vector abs should reject");
    ok &= check(throws_runtime([&] {
        (void)ad::min(X, Y).reverse_diff(ad::vars(X, Y), lambda);
    }), "reverse derivative of vector min should reject");
    return ok;
}

bool test_sparsity_and_codegen() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::Function f({X}, ad::Vec{
        ad::sqrt(X[0]) + ad::log10(X[1]) + ad::pow(X[1], X[2]),
        ad::min(X[0], X[1]) + ad::max(X[1], X[2]),
        ad::abs(X[2]) + ad::sinh(X[0]) + ad::cosh(X[1]) + ad::atan(X[2]),
    });
    ad::SparsityPattern sp = f.jacobian_sparsity();

    bool ok = true;
    ok &= check(sp.rows() == 3 && sp.cols() == 3, "elementary sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 1) && sp.contains(0, 2), "elementary first row sparsity is wrong");
    ok &= check(sp.contains(1, 0) && sp.contains(1, 1) && sp.contains(1, 2), "elementary min/max sparsity is wrong");
    ok &= check(sp.contains(2, 0) && sp.contains(2, 1) && sp.contains(2, 2), "elementary third row sparsity is wrong");
    ok &= check(compile_and_run("ad_elementary_ops", f, {4.0, 3.0, 2.0}), "generated C for elementary ops should match native eval");

    ad::Function smooth({X}, ad::Vec{
        ad::sqrt(X[0]) + ad::log10(X[1]) + ad::pow(X[1], X[2]),
        ad::sinh(X[0]) + ad::cosh(X[1]) + ad::atan(X[2]) + ad::sigmoid(X[0]),
    });
    ok &= check(compile_and_run("ad_elementary_ops_jvp", smooth.forward_function(), {4.0, 3.0, 2.0, 0.1, -0.2, 0.3}), "generated C for smooth elementary JVP should match native eval");
    ok &= check(compile_and_run("ad_elementary_ops_vjp", smooth.reverse_function(), {4.0, 3.0, 2.0, 1.0, -0.5}), "generated C for smooth elementary VJP should match native eval");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_primal_eval();
    ok &= test_vector_primal_eval();
    ok &= test_smooth_forward_reverse_eval();
    ok &= test_nonsmooth_derivatives_reject();
    ok &= test_sparsity_and_codegen();
    if (!ok) {
        return 1;
    }
    std::cout << "AD elementary operation tests passed\n";
    return 0;
}
