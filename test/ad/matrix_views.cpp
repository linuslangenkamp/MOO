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

std::vector<double> eval_vec(const ad::Vec &vec, const ad::Values &values) {
    ad::EvalWorkspace work;
    std::vector<double> out(static_cast<std::size_t>(vec.size()), 0.0);
    vec.eval(values, work, out.data(), static_cast<int>(out.size()));
    return out;
}

double eval_expr(const ad::Expr &expr, const ad::Values &values) {
    ad::EvalWorkspace work;
    double out = 0.0;
    expr.eval(values, work, &out);
    return out;
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

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_matrix_views_tmp";
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

bool test_construction_layout_and_validation() {
    ad::Mat A = ad::mat_variable("A", 2, 3);
    ad::Mat R = ad::mat_variable("R", 2, 3, ad::MatrixLayout::RowMajor);
    ad::Mat P = ad::mat_parameter("P", 2, 2);
    ad::Mat C = ad::mat_constant(2, 2, {1.0, 2.0, 3.0, 4.0});

    ad::Values values;
    values.set(A.vec(), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    values.set(R.vec(), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    values.set(P.vec(), {7.0, 8.0, 9.0, 10.0});

    bool ok = true;
    ok &= check(A.valid() && A.rows() == 2 && A.cols() == 3 && A.size() == 6, "column-major matrix metadata is wrong");
    ok &= check(A.layout() == ad::MatrixLayout::ColumnMajor, "default matrix layout should be column-major");
    ok &= check(R.layout() == ad::MatrixLayout::RowMajor, "row-major matrix layout is wrong");
    ok &= check(near(eval_expr(A(0, 1), values), 3.0) && near(eval_expr(A(1, 2), values), 6.0), "column-major element indexing is wrong");
    ok &= check(near(eval_expr(R(0, 1), values), 2.0) && near(eval_expr(R(1, 0), values), 4.0), "row-major element indexing is wrong");
    ok &= check(near(eval_expr(P(1, 0), values), 8.0), "matrix parameter values are wrong");
    ok &= check(near_array(eval_vec(C.vec(), values), {1.0, 2.0, 3.0, 4.0}), "matrix constant values are wrong");
    ok &= check(throws_runtime([] { (void)ad::mat(ad::vec_variable("bad", 3), 2, 2); }), "matrix dimension mismatch should throw");
    ok &= check(throws_runtime([&] { (void)A(2, 0); }), "matrix row out-of-bounds should throw");
    ok &= check(throws_runtime([&] { (void)A.col(3); }), "matrix column out-of-bounds should throw");
    return ok;
}

bool test_views_and_layout_conversion() {
    ad::Mat A = ad::mat_variable("A", 2, 3);
    ad::Values values;
    values.set(A.vec(), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

    ad::Vec row = A.row(1);
    ad::Vec col = A.col(2);
    ad::Mat block = A.block(0, 1, 2, 2);
    ad::Mat transpose = A.transpose();
    ad::Mat row_major = A.as_layout(ad::MatrixLayout::RowMajor);
    ad::Mat reshaped = A.reshape(3, 2);

    bool ok = true;
    ok &= check(near_array(eval_vec(row, values), {2.0, 4.0, 6.0}), "matrix row view is wrong");
    ok &= check(near_array(eval_vec(col, values), {5.0, 6.0}), "matrix column view is wrong");
    ok &= check(block.rows() == 2 && block.cols() == 2 && near_array(eval_vec(block.vec(), values), {3.0, 4.0, 5.0, 6.0}), "matrix block view is wrong");
    ok &= check(transpose.rows() == 3 && transpose.cols() == 2 && near_array(eval_vec(transpose.vec(), values), {1.0, 3.0, 5.0, 2.0, 4.0, 6.0}), "matrix transpose view is wrong");
    ok &= check(row_major.layout() == ad::MatrixLayout::RowMajor && near_array(eval_vec(row_major.vec(), values), {1.0, 3.0, 5.0, 2.0, 4.0, 6.0}), "matrix layout conversion is wrong");
    ok &= check(reshaped.rows() == 3 && reshaped.cols() == 2 && near_array(eval_vec(reshaped.vec(), values), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0}), "matrix reshape should preserve flat payload");
    ok &= check(contains_kind(ad::inspect(row_major.vec()), ad::GraphNodeKind::Gather), "matrix layout conversion should be a structured gather");
    ok &= check(throws_runtime([&] { (void)A.block(1, 2, 2, 1); }), "out-of-bounds matrix block should throw");
    ok &= check(throws_runtime([&] { (void)A.reshape(4, 2); }), "size-changing matrix reshape should throw");
    return ok;
}

bool test_elementwise_operations() {
    ad::Mat A = ad::mat_variable("A", 2, 3);
    ad::Mat B = A.as_layout(ad::MatrixLayout::RowMajor);
    ad::Expr alpha = ad::parameter("alpha");
    ad::Mat sum = A + B;
    ad::Mat scaled = alpha * A;
    ad::Mat unary = ad::sqrt(ad::hadamard(A, A) + 1.0);
    ad::Mat powered = ad::pow(A + 2.0, 2.0);

    ad::Values values;
    values.set(A.vec(), {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    values.set(alpha.param(), 0.5);

    bool ok = true;
    ok &= check(sum.layout() == ad::MatrixLayout::ColumnMajor, "mixed-layout elementwise op should return lhs layout");
    ok &= check(near_array(eval_vec(sum.vec(), values), {2.0, 4.0, 6.0, 8.0, 10.0, 12.0}), "mixed-layout elementwise addition is wrong");
    ok &= check(near_array(eval_vec(scaled.vec(), values), {0.5, 1.0, 1.5, 2.0, 2.5, 3.0}), "Expr-scaled matrix is wrong");
    ok &= check(near_array(eval_vec(unary.vec(), values), {
        std::sqrt(2.0), std::sqrt(5.0), std::sqrt(10.0), std::sqrt(17.0), std::sqrt(26.0), std::sqrt(37.0)}), "matrix unary wrapper is wrong");
    ok &= check(near_array(eval_vec(powered.vec(), values), {9.0, 16.0, 25.0, 36.0, 49.0, 64.0}), "matrix pow wrapper is wrong");
    ok &= check(contains_kind(ad::inspect(B.vec()), ad::GraphNodeKind::Gather), "explicit matrix layout conversion should be a structured gather");
    ok &= check(throws_runtime([&] { (void)(A + ad::mat_variable("C", 3, 2)); }), "matrix shape mismatch should throw");
    return ok;
}

bool test_function_ad_sparsity_eval_and_codegen() {
    ad::Mat X = ad::mat_variable("X", 2, 2);
    ad::Mat Y = ad::sin(X) + X.transpose();
    ad::Function f({X.vec()}, std::vector<ad::Vec>{Y.vec(), X.row(0)});
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();
    ad::SparsityPattern sp = f.jacobian_sparsity();
    const std::vector<double> out = eval_function(f, {1.0, 2.0, 3.0, 4.0});

    bool ok = true;
    ok &= check(f.output_count() == 2 && f.output_size() == 6, "matrix-view function output grouping is wrong");
    ok &= check(df.output_count() == 2 && df.output_size() == f.output_size(), "matrix-view forward function grouping is wrong");
    ok &= check(rf.output_count() == 1 && rf.output_size() == f.input_size(), "matrix-view reverse function grouping is wrong");
    ok &= check(near_array(out, {
        std::sin(1.0) + 1.0,
        std::sin(2.0) + 3.0,
        std::sin(3.0) + 2.0,
        std::sin(4.0) + 4.0,
        1.0,
        3.0}), "matrix-view function eval is wrong");
    ok &= check(sp.rows() == 6 && sp.cols() == 4, "matrix-view function sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0), "matrix-view sparsity row 0 is wrong");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 2), "matrix transpose sparsity row 1 is wrong");
    ok &= check(sp.contains(2, 1) && sp.contains(2, 2), "matrix transpose sparsity row 2 is wrong");
    ok &= check(sp.contains(4, 0) && sp.contains(5, 2), "matrix row output sparsity is wrong");
    ok &= check(contains_kind(f.info().output_graph, ad::GraphNodeKind::Gather), "matrix-view function should contain structured gather views");
    ok &= check(compile_and_run("ad_matrix_view", f, {1.0, 2.0, 3.0, 4.0}), "generated C for matrix-view function should match native eval");
    ok &= check(compile_and_run("ad_matrix_view_jvp", df, {1.0, 2.0, 3.0, 4.0, 0.1, 0.2, 0.3, 0.4}), "generated C for matrix-view JVP should match native eval");
    ok &= check(compile_and_run("ad_matrix_view_vjp", rf, {1.0, 2.0, 3.0, 4.0, 0.5, 0.25, -0.75, 1.5, -2.0, 3.0}), "generated C for matrix-view VJP should match native eval");
    return ok;
}

bool test_multi_output_matrix_reconstruction() {
    ad::Mat X = ad::mat_variable("X", 2, 2);
    ad::Function f({X.vec()}, std::vector<ad::Vec>{X.transpose().vec(), X.as_layout(ad::MatrixLayout::RowMajor).vec()});
    ad::Mat transpose_out = ad::mat(f.output(0), 2, 2);
    ad::Mat row_major_out = ad::mat(f.output(1), 2, 2, ad::MatrixLayout::RowMajor);

    ad::Values values;
    values.set(X.vec(), {1.0, 2.0, 3.0, 4.0});

    bool ok = true;
    ok &= check(f.output_count() == 2 && f.output_offset(1) == 4, "matrix multi-output offsets are wrong");
    ok &= check(near(eval_expr(transpose_out(0, 1), values), 2.0), "matrix reconstruction from first output group is wrong");
    ok &= check(near(eval_expr(row_major_out(1, 0), values), 2.0), "row-major matrix reconstruction from second output group is wrong");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_construction_layout_and_validation();
    ok &= test_views_and_layout_conversion();
    ok &= test_elementwise_operations();
    ok &= test_function_ad_sparsity_eval_and_codegen();
    ok &= test_multi_output_matrix_reconstruction();
    if (!ok) {
        return 1;
    }
    std::cout << "test_ad_matrix_views passed\n";
    return 0;
}
