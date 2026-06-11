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

bool near(double lhs, double rhs, double tolerance = 1e-9) {
    return std::fabs(lhs - rhs) <= tolerance;
}

bool near_array(const std::vector<double> &values, const std::vector<double> &expected, double tolerance = 1e-9) {
    if (values.size() != expected.size()) {
        return false;
    }
    for (std::size_t i = 0; i < values.size(); ++i) {
        if (!near(values[i], expected[i], tolerance)) {
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

std::vector<double> eval_function(const ad::Function &function, const std::vector<double> &input) {
    ad::EvalWorkspace work;
    std::vector<double> out(static_cast<std::size_t>(function.output_size()), 0.0);
    function.eval(input.empty() ? nullptr : input.data(),
                  static_cast<int>(input.size()),
                  nullptr,
                  0,
                  work,
                  out.empty() ? nullptr : out.data(),
                  static_cast<int>(out.size()));
    return out;
}

bool compile_and_run(const std::string &name, const ad::Function &function, const std::vector<double> &input) {
    ad::CodegenOptions options;
    options.function_name = name;
    const ad::CodegenResult generated = ad::generate_c(function, options);
    const std::vector<double> expected = eval_function(function, input);
    if (generated.source.find("moo_ad_lapack_linear_solve") == std::string::npos) {
        std::cerr << "FAIL: generated solve code is missing LAPACK linear solver backend call\n";
        return false;
    }

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_linear_solve_tmp";
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
    file << "        if (fabs(out[i] - expected[i]) > 1e-9) {\n";
    file << "            printf(\"mismatch %d %.17g %.17g\\n\", i, out[i], expected[i]);\n";
    file << "            return 20;\n";
    file << "        }\n";
    file << "    }\n";
    file << "    return 0;\n";
    file << "}\n";
    file.close();

    const std::string compile_command = shell_quote(MOO_CXX_COMPILER) + " " +
                                        shell_quote(source_path.string()) + " -lm -o " +
                                        shell_quote(exe_path.string()) + " " +
                                        shell_quote(MOO_AD_LIBRARY) + " -Wl,-rpath," +
                                        shell_quote(MOO_AD_LIBRARY_DIR);
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

bool test_primal_eval_transpose_and_validation() {
    ad::Mat A = ad::mat_variable("A", 2, 2);
    ad::Vec b = ad::vec_variable("b", 2);
    ad::Vec y = ad::solve(A, b);
    ad::Vec yt = ad::solve_transpose(A, b);
    ad::Function f({A.vec(), b}, ad::concat(y, yt));

    const std::vector<double> input = {
        2.0, 4.0, 1.0, 3.0,
        1.0, 2.0,
    };
    const std::vector<double> out = eval_function(f, input);

    bool ok = true;
    ok &= check(near_array(out, {0.5, 0.0, -2.5, 1.5}), "linear solve or transpose solve eval is wrong");
    ok &= check(contains_kind(ad::inspect(y), ad::GraphNodeKind::LinearSolve), "solve graph should contain LinearSolve");
    ok &= check(ad::inspect(yt).op == "solve_transpose", "transpose solve inspection op is wrong");
    ok &= check(throws_runtime([&] {
        (void)ad::solve(ad::mat_variable("bad", 2, 3), ad::vec_variable("rhs", 2));
    }), "non-square solve matrix should throw");
    ok &= check(throws_runtime([&] {
        (void)ad::solve(A, ad::vec_variable("wrong", 3));
    }), "solve rhs size mismatch should throw");
    return ok;
}

bool test_forward_reverse_and_sparsity() {
    ad::Mat A = ad::mat_variable("A", 2, 2);
    ad::Vec b = ad::vec_variable("b", 2);
    ad::Vec seed = ad::vec_variable("seed", 6);
    ad::Vec lambda = ad::vec_variable("lambda", 2);
    ad::Vec y = ad::solve(A, b);
    ad::Vec dy = y.forward_diff(ad::vars(A.vec(), b), seed);
    ad::Vec adj = y.reverse_diff(ad::vars(A.vec(), b), lambda);

    ad::Function jvp({A.vec(), b, seed}, dy);
    ad::Function vjp({A.vec(), b, lambda}, adj);

    const std::vector<double> primal = {
        2.0, 4.0, 1.0, 3.0,
        1.0, 2.0,
    };
    const std::vector<double> tangent = {
        0.1, 0.3, 0.2, 0.4,
        0.5, -0.25,
    };
    std::vector<double> jvp_input = primal;
    jvp_input.insert(jvp_input.end(), tangent.begin(), tangent.end());
    const std::vector<double> jvp_out = eval_function(jvp, jvp_input);

    std::vector<double> vjp_input = primal;
    vjp_input.push_back(1.0);
    vjp_input.push_back(-2.0);
    const std::vector<double> vjp_out = eval_function(vjp, vjp_input);

    const ad::SparsityPattern sp = ad::sparsity(y, ad::vars(A.vec(), b));

    bool ok = true;
    ok &= check(near_array(jvp_out, {0.875, -1.3}), "linear solve JVP eval is wrong");
    ok &= check(near_array(vjp_out, {-2.75, 1.25, 0.0, 0.0, 5.5, -2.5}), "linear solve VJP eval is wrong");
    ok &= check(contains_kind(ad::inspect(dy), ad::GraphNodeKind::LinearSolve), "linear solve JVP should contain LinearSolve");
    ok &= check(contains_kind(ad::inspect(adj), ad::GraphNodeKind::LinearSolve), "linear solve VJP should contain LinearSolve");
    ok &= check(sp.rows() == 2 && sp.cols() == 6 && sp.nnz() == 12, "linear solve sparsity should be conservative dense");
    for (int row = 0; row < 2; ++row) {
        for (int col = 0; col < 6; ++col) {
            ok &= check(sp.contains(row, col), "linear solve conservative sparsity is missing an entry");
        }
    }
    return ok;
}

bool test_transpose_derivatives_hvp_and_codegen() {
    ad::Mat A = ad::mat_variable("A", 2, 2);
    ad::Vec b = ad::vec_variable("b", 2);
    ad::Vec y = ad::solve_transpose(A, b);
    ad::Vec seed = ad::vec_variable("seed", 6);
    ad::Vec lambda = ad::vec_variable("lambda", 2);
    ad::Vec dy = y.forward_diff(ad::vars(A.vec(), b), seed);
    ad::Vec adj = y.reverse_diff(ad::vars(A.vec(), b), lambda);
    ad::Vec hvp = adj.forward_diff(ad::vars(A.vec(), b), seed);

    ad::Function jvp({A.vec(), b, seed}, dy);
    ad::Function h({A.vec(), b, lambda, seed}, hvp);
    ad::Function generated({A.vec(), b}, ad::concat(ad::solve(A, b), y));

    const std::vector<double> primal = {
        2.0, 4.0, 1.0, 3.0,
        1.0, 2.0,
    };
    const std::vector<double> tangent = {
        0.1, 0.3, 0.2, 0.4,
        0.5, -0.25,
    };
    std::vector<double> jvp_input = primal;
    jvp_input.insert(jvp_input.end(), tangent.begin(), tangent.end());
    std::vector<double> hvp_input = primal;
    hvp_input.push_back(1.0);
    hvp_input.push_back(-2.0);
    hvp_input.insert(hvp_input.end(), tangent.begin(), tangent.end());

    bool ok = true;
    ok &= check(near_array(eval_function(jvp, jvp_input), {1.15, -0.5}), "transpose solve JVP eval is wrong");
    ok &= check(h.output_size() == 6 && eval_function(h, hvp_input).size() == 6, "linear solve HVP graph/eval shape is wrong");
    ok &= check(contains_kind(ad::inspect(hvp), ad::GraphNodeKind::LinearSolve), "linear solve HVP should contain LinearSolve");
    ok &= check(compile_and_run("linear_solve_codegen", generated, primal), "linear solve codegen failed");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_primal_eval_transpose_and_validation();
    ok &= test_forward_reverse_and_sparsity();
    ok &= test_transpose_derivatives_hvp_and_codegen();
    if (!ok) {
        return 1;
    }
    std::cout << "linear solve tests passed\n";
    return 0;
}
