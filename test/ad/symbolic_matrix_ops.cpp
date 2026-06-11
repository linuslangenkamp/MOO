// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
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

bool near_array(const std::vector<double> &values, const std::vector<double> &expected) {
    if (values.size() != expected.size()) {
        return false;
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!near(values[i], expected[i])) {
            return false;
        }
    }
    return true;
}

template <typename Func>
bool throws_runtime(Func &&call) {
    try {
        call();
    } catch (const std::runtime_error &) {
        return true;
    }
    return false;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    for (const ad::GraphInfo &child : info.children) {
        if (contains_kind(child, kind)) {
            return true;
        }
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

std::vector<double> eval_function(const ad::Function &function,
                                  const std::vector<double> &input,
                                  const std::vector<double> &params = {}) {
    ad::EvalWorkspace work;
    std::vector<double> out(static_cast<std::size_t>(function.output_size()), 0.0);
    function.eval(input.empty() ? nullptr : input.data(),
                  static_cast<int>(input.size()),
                  params.empty() ? nullptr : params.data(),
                  static_cast<int>(params.size()),
                  work,
                  out.data(),
                  static_cast<int>(out.size()));
    return out;
}

bool compile_and_run(const std::string &name,
                     const ad::Function &function,
                     const std::vector<double> &input,
                     const std::vector<double> &params = {}) {
    ad::CodegenOptions options;
    options.function_name = name;
    const ad::CodegenResult generated = ad::generate_c(function, options);
    const std::vector<double> expected = eval_function(function, input, params);

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_symbolic_matrix_tmp";
    std::filesystem::create_directories(dir);
    const std::filesystem::path source_path = dir / (name + ".c");
    const std::filesystem::path exe_path = dir / name;

    std::ofstream file(source_path);
    file << generated.source << "\n";
    file << "#include <math.h>\n";
    file << "#include <stdio.h>\n";
    file << "int main(void) {\n";
    file << "    double in[" << input.size() << "] = " << c_array(input) << ";\n";
    file << "    double p[" << params.size() << "] = " << c_array(params) << ";\n";
    file << "    double out[" << expected.size() << "] = {0};\n";
    file << "    " << name << "(in, p, out);\n";
    file << "    double expected[" << expected.size() << "] = " << c_array(expected) << ";\n";
    file << "    for (int i = 0; i < " << expected.size() << "; ++i) {\n";
    file << "        if (fabs(out[i] - expected[i]) > 1e-9) {\n";
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

bool test_symbolic_matvec() {
    ad::Mat A = ad::mat_variable("A", 2, 3);
    ad::Vec x = ad::vec_variable("x", 3);
    ad::Vec y = A * x;
    ad::Function f({A.vec(), x}, y);
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();
    const std::vector<double> input = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0};
    const std::vector<double> lambda_input = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 0.5, -1.0};
    const ad::SparsityPattern sp = f.jacobian_sparsity();
    const std::vector<double> vjp = eval_function(rf, lambda_input);

    bool ok = true;
    ok &= check(y.size() == 2 && contains_kind(ad::inspect(y), ad::GraphNodeKind::SymbolicMatVec), "symbolic matvec graph shape is wrong");
    ok &= check(near_array(eval_function(f, input), {76.0, 100.0}), "symbolic matvec eval is wrong");
    ok &= check(sp.rows() == 2 && sp.cols() == 9, "symbolic matvec sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2) && sp.contains(0, 4), "symbolic matvec matrix-row sparsity is wrong");
    ok &= check(sp.contains(0, 6) && sp.contains(0, 7) && sp.contains(0, 8), "symbolic matvec rhs sparsity is wrong");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 3) && sp.contains(1, 5), "symbolic matvec second-row matrix sparsity is wrong");
    ok &= check(near_array(vjp, {3.5, -7.0, 4.0, -8.0, 4.5, -9.0, -1.5, -2.5, -3.5}), "symbolic matvec reverse eval is wrong");
    ok &= check(compile_and_run("ad_symbolic_matvec", f, input), "generated C for symbolic matvec should match native eval");
    ok &= check(compile_and_run("ad_symbolic_matvec_jvp", df, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0,
                                                               0.1, 0.2, 0.3, 0.4, 0.5, 0.6, -1.0, 2.0, -3.0}), "generated C for symbolic matvec JVP should match native eval");
    ok &= check(compile_and_run("ad_symbolic_matvec_vjp", rf, lambda_input), "generated C for symbolic matvec VJP should match native eval");
    return ok;
}

bool test_symbolic_matmul() {
    ad::Mat A = ad::mat_variable("A", 2, 3);
    ad::Mat B = ad::mat_variable("B", 3, 2, ad::MatrixLayout::RowMajor);
    ad::Mat C = A * B;
    ad::Function f({A.vec(), B.vec()}, C.vec());
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();
    const std::vector<double> input = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
    const ad::SparsityPattern sp = f.jacobian_sparsity();

    bool ok = true;
    ok &= check(C.rows() == 2 && C.cols() == 2 && contains_kind(ad::inspect(C.vec()), ad::GraphNodeKind::SymbolicMatMul), "symbolic matmul graph shape is wrong");
    ok &= check(near_array(eval_function(f, input), {89.0, 116.0, 98.0, 128.0}), "symbolic matmul eval is wrong");
    ok &= check(sp.rows() == 4 && sp.cols() == 12, "symbolic matmul sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2) && sp.contains(0, 4), "symbolic matmul lhs sparsity is wrong");
    ok &= check(sp.contains(0, 6) && sp.contains(0, 8) && sp.contains(0, 10), "symbolic matmul rhs sparsity is wrong");
    ok &= check(throws_runtime([&] { (void)(A * A); }), "invalid matrix multiplication shape should throw");
    ok &= check(compile_and_run("ad_symbolic_matmul", f, input), "generated C for symbolic matmul should match native eval");
    ok &= check(compile_and_run("ad_symbolic_matmul_jvp", df, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                                               0.1, 0.2, 0.3, 0.4, 0.5, 0.6, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0}), "generated C for symbolic matmul JVP should match native eval");
    ok &= check(compile_and_run("ad_symbolic_matmul_vjp", rf, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0,
                                                               0.5, -1.0, 1.5, -2.0}), "generated C for symbolic matmul VJP should match native eval");
    return ok;
}

bool test_outer_bilinear_quadratic() {
    ad::Vec a = ad::vec_variable("a", 2);
    ad::Vec b = ad::vec_variable("b", 3);
    ad::Mat O = ad::outer(a, b);
    ad::Function outer_fun({a, b}, O.vec());

    ad::Vec x = ad::vec_variable("x", 2);
    ad::Mat A = ad::mat_variable("A", 2, 2);
    ad::Vec y = ad::vec_variable("y", 2);
    ad::Expr bilin = ad::bilinear_form(x, A, y);
    ad::Expr quad = ad::quadratic_form(x, A);
    ad::Function f({x, A.vec(), y}, ad::Vec{bilin, quad});
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();
    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Vec seed = ad::vec_variable("seed", f.input_size());
    ad::Vec hvp = f.reverse(lambda).forward_diff(f.input_vars(), seed);

    bool ok = true;
    ok &= check(contains_kind(ad::inspect(O.vec()), ad::GraphNodeKind::OuterProduct), "outer product should be structured");
    ok &= check(near_array(eval_function(outer_fun, {2.0, 3.0, 5.0, 7.0, 11.0}), {10.0, 15.0, 14.0, 21.0, 22.0, 33.0}), "outer product eval is wrong");
    ok &= check(contains_kind(f.info().output_graph, ad::GraphNodeKind::SymbolicMatVec), "bilinear/quadratic forms should contain symbolic matvec");
    ok &= check(near_array(eval_function(f, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}), {213.0, 45.0}), "bilinear/quadratic eval is wrong");
    ok &= check(hvp.size() == f.input_size(), "bilinear/quadratic HVP composition size is wrong");
    ok &= check(compile_and_run("ad_outer_product", outer_fun, {2.0, 3.0, 5.0, 7.0, 11.0}), "generated C for outer product should match native eval");
    ok &= check(compile_and_run("ad_bilinear_quad", f, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0}), "generated C for bilinear/quadratic function should match native eval");
    ok &= check(compile_and_run("ad_bilinear_quad_jvp", df, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0,
                                                            0.1, -0.2, 0.3, -0.4, 0.5, -0.6, 0.7, -0.8}), "generated C for bilinear/quadratic JVP should match native eval");
    ok &= check(compile_and_run("ad_bilinear_quad_vjp", rf, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 0.25, -1.5}), "generated C for bilinear/quadratic VJP should match native eval");
    return ok;
}

bool test_explicit_elementwise_helpers() {
    ad::Mat A = ad::mat_variable("A", 2, 2);
    ad::Mat B = A.as_layout(ad::MatrixLayout::RowMajor);
    ad::Mat H = ad::hadamard(A, B);
    ad::Mat D = ad::elem_div(A + 1.0, B + 2.0);
    ad::Function f({A.vec()}, std::vector<ad::Vec>{H.vec(), D.vec()});

    bool ok = true;
    ok &= check(H.rows() == 2 && H.cols() == 2 && H.layout() == A.layout(), "hadamard metadata is wrong");
    ok &= check(near_array(eval_function(f, {1.0, 2.0, 3.0, 4.0}), {1.0, 4.0, 9.0, 16.0, 2.0 / 3.0, 3.0 / 4.0, 4.0 / 5.0, 5.0 / 6.0}), "explicit elementwise matrix helpers evaluate incorrectly");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_symbolic_matvec();
    ok &= test_symbolic_matmul();
    ok &= test_outer_bilinear_quadratic();
    ok &= test_explicit_elementwise_helpers();
    if (!ok) {
        return 1;
    }
    std::cout << "test_ad_symbolic_matrix_ops passed\n";
    return 0;
}
