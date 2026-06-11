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

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_symbolic_sparse_matrix_tmp";
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

bool test_construction_and_dense_conversion() {
    ad::SparseMat A = ad::sparse_mat_variable("A", 3, 3, {2, 0, 1}, {1, 2, 0});
    ad::SparseMat C = ad::sparse_mat_constant(3, 3, {2, 0, 1}, {1, 2, 0}, {10.0, 20.0, 30.0});
    ad::Mat dense = C.to_dense();
    ad::Function dense_fun({}, dense.vec());

    bool ok = true;
    ok &= check(A.nnz() == 3 && A.values().node_kind() == ad::GraphNodeKind::VectorVariable, "sparse variable values should stay a direct variable group");
    ok &= check(A.row_indices() == std::vector<int>({0, 1, 2}), "sparse variable rows should be canonical");
    ok &= check(A.col_indices() == std::vector<int>({2, 0, 1}), "sparse variable cols should be canonical");
    ok &= check(near_array(eval_function(dense_fun, {}), {0.0, 30.0, 0.0, 0.0, 0.0, 10.0, 20.0, 0.0, 0.0}), "sparse to_dense eval is wrong");
    ok &= check(contains_kind(ad::inspect(A.to_dense().vec()), ad::GraphNodeKind::ScatterAdd), "symbolic sparse to_dense should use structured scatter_add");
    ok &= check(throws_runtime([&] { (void)ad::sparse_mat_variable("bad", 2, 2, {0, 0}, {1, 1}); }), "duplicate sparse entries should throw");
    ok &= check(throws_runtime([&] { (void)ad::sparse_mat_variable("bad", 2, 2, {0, 2}, {1, 1}); }), "out-of-bounds sparse entries should throw");
    return ok;
}

bool test_sparse_symbolic_matvec() {
    ad::SparseMat A = ad::sparse_mat_variable("A", 2, 3, {0, 1, 1}, {0, 1, 2});
    ad::Vec x = ad::vec_variable("x", 3);
    ad::Vec y = A * x;
    ad::Function f({A.values(), x}, y);
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();
    const ad::SparsityPattern sp = f.jacobian_sparsity();
    const std::vector<double> input = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
    const std::vector<double> reverse_input = {2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.5, -1.0};

    bool ok = true;
    ok &= check(contains_kind(ad::inspect(y), ad::GraphNodeKind::SymbolicSparseMatVec), "sparse symbolic matvec should preserve sparse matvec node");
    ok &= check(near_array(eval_function(f, input), {10.0, 46.0}), "sparse symbolic matvec eval is wrong");
    ok &= check(sp.rows() == 2 && sp.cols() == 6, "sparse symbolic matvec sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 3), "row 0 sparse matvec sparsity is wrong");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 2) && sp.contains(1, 4) && sp.contains(1, 5), "row 1 sparse matvec sparsity is wrong");
    ok &= check(near_array(eval_function(rf, reverse_input), {2.5, -6.0, -7.0, 1.0, -3.0, -4.0}), "sparse symbolic matvec reverse eval is wrong");
    ok &= check(compile_and_run("ad_symbolic_sparse_matvec", f, input), "generated C for sparse symbolic matvec should match native eval");
    ok &= check(compile_and_run("ad_symbolic_sparse_matvec_jvp", df, {2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                                                                       0.1, 0.2, 0.3, -1.0, 2.0, -3.0}), "generated C for sparse symbolic matvec JVP should match native eval");
    ok &= check(compile_and_run("ad_symbolic_sparse_matvec_vjp", rf, reverse_input), "generated C for sparse symbolic matvec VJP should match native eval");
    return ok;
}

bool test_sparse_dense_matmul() {
    ad::SparseMat A = ad::sparse_mat_variable("A", 2, 3, {0, 1, 1}, {0, 1, 2});
    ad::Mat B = ad::mat_variable("B", 3, 2);
    ad::Mat C = A * B;
    ad::Function f({A.values(), B.vec()}, C.vec());
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();
    const std::vector<double> input = {2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const std::vector<double> reverse_input = {2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.5, -1.0, 1.5, -2.0};

    bool ok = true;
    ok &= check(C.rows() == 2 && C.cols() == 2 && contains_kind(ad::inspect(C.vec()), ad::GraphNodeKind::SymbolicSparseMatMul), "sparse*dense matmul graph is wrong");
    ok &= check(near_array(eval_function(f, input), {2.0, 18.0, 8.0, 39.0}), "sparse*dense matmul eval is wrong");
    ok &= check(near_array(eval_function(rf, reverse_input), {6.5, -12.0, -15.0, 1.0, -3.0, -4.0, 3.0, -6.0, -8.0}), "sparse*dense matmul reverse eval is wrong");
    ok &= check(compile_and_run("ad_symbolic_sparse_dense_matmul", f, input), "generated C for sparse*dense matmul should match native eval");
    ok &= check(compile_and_run("ad_symbolic_sparse_dense_matmul_jvp", df, {2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
                                                                            0.1, 0.2, 0.3, -1.0, 2.0, -3.0, 4.0, -5.0, 6.0}), "generated C for sparse*dense matmul JVP should match native eval");
    ok &= check(compile_and_run("ad_symbolic_sparse_dense_matmul_vjp", rf, reverse_input), "generated C for sparse*dense matmul VJP should match native eval");
    return ok;
}

bool test_dense_sparse_matmul_and_sparse_ops() {
    ad::Mat D = ad::mat_variable("D", 2, 2);
    ad::SparseMat S = ad::sparse_mat_variable("S", 2, 3, {0, 1}, {1, 2});
    ad::Mat Y = D * S;
    ad::Function f({D.vec(), S.values()}, Y.vec());

    ad::SparseMat P = ad::sparse_mat_parameter("P", 2, 3, {0, 1}, {1, 2});
    ad::SparseMat H = ad::hadamard(S + P, S);
    ad::Function h({S.values()}, H.values());
    const std::vector<double> input = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

    bool ok = true;
    ok &= check(Y.rows() == 2 && Y.cols() == 3 && contains_kind(ad::inspect(Y.vec()), ad::GraphNodeKind::SymbolicSparseMatMul), "dense*sparse matmul graph is wrong");
    ok &= check(near_array(eval_function(f, input), {0.0, 0.0, 5.0, 10.0, 18.0, 24.0}), "dense*sparse matmul eval is wrong");
    ok &= check(H.nnz() == 2, "sparse hadamard should keep intersection pattern");
    ok &= check(h.parameters().size() == 2, "sparse parameter values should be collected through sparse ops");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    const auto run = [&](const char *name, bool (*test)()) {
        try {
            ok &= test();
        } catch (const std::exception &ex) {
            std::cerr << "FAIL: " << name << " threw: " << ex.what() << '\n';
            ok = false;
        }
    };
    run("construction_and_dense_conversion", test_construction_and_dense_conversion);
    run("sparse_symbolic_matvec", test_sparse_symbolic_matvec);
    run("sparse_dense_matmul", test_sparse_dense_matmul);
    run("dense_sparse_matmul_and_sparse_ops", test_dense_sparse_matmul_and_sparse_ops);
    return ok ? 0 : 1;
}
