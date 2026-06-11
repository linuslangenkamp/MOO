// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
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

bool throws_runtime(void (*call)()) {
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

void collect_unique_kind_ids(const ad::GraphInfo &info, ad::GraphNodeKind kind, std::vector<ad::NodeId> &ids) {
    if (info.kind == kind) {
        bool seen = false;
        for (ad::NodeId id : ids) {
            if (id == info.id) {
                seen = true;
                break;
            }
        }
        if (!seen) {
            ids.push_back(info.id);
        }
    }
    for (const ad::GraphInfo &child : info.children) {
        collect_unique_kind_ids(child, kind, ids);
    }
}

int count_unique_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    std::vector<ad::NodeId> ids;
    collect_unique_kind_ids(info, kind, ids);
    return static_cast<int>(ids.size());
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

std::vector<double> eval_function(const ad::Function &function, const std::vector<double> &input, const std::vector<double> &params = {}) {
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

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_multi_output_tmp";
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

bool test_single_output_compatibility() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Function f({X}, ad::sin(X));
    ad::FunctionInfo info = f.info();

    bool ok = true;
    ok &= check(f.output_count() == 1, "single-output function should report one output group");
    ok &= check(f.outputs().size() == 1, "single-output function outputs() size is wrong");
    ok &= check(f.output_size() == 2 && f.output(0).size() == 2, "single-output function sizes are wrong");
    ok &= check(f.output_offset(0) == 0, "single-output offset should be zero");
    ok &= check(f.output().node_id() == f.output(0).node_id(), "single-output flat output should reuse the group output");
    ok &= check(info.output_count == 1 && info.outputs.size() == 1, "single-output FunctionInfo output count is wrong");
    ok &= check(info.outputs[0].offset == 0 && info.outputs[0].size == 2, "single-output FunctionInfo output layout is wrong");
    return ok;
}

bool test_multi_output_metadata_eval_sparsity() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec U = ad::vec_variable("U", 1);
    ad::Expr h = ad::parameter("h");
    ad::Function f({X, U}, ad::params(h), std::vector<ad::Vec>{
        ad::Vec{X[0] + U[0], X[1] + U[0]},
        ad::Vec{h * X[0], X[1] * U[0]},
        ad::vec_constant({7.0}),
    });
    ad::FunctionInfo info = f.info();
    const std::vector<double> out = eval_function(f, {1.0, 2.0, 3.0}, {0.5});
    ad::SparsityPattern sp = f.jacobian_sparsity();

    bool ok = true;
    ok &= check(f.output_count() == 3 && f.output_size() == 5, "multi-output layout size is wrong");
    ok &= check(info.output_count == 3 && info.outputs.size() == 3, "multi-output FunctionInfo count is wrong");
    ok &= check(info.outputs[0].offset == 0 && info.outputs[0].size == 2, "first output layout is wrong");
    ok &= check(info.outputs[1].offset == 2 && info.outputs[1].size == 2, "second output layout is wrong");
    ok &= check(info.outputs[2].offset == 4 && info.outputs[2].size == 1, "third output layout is wrong");
    ok &= check(f.output(0).size() == 2 && f.output(1).size() == 2 && f.output(2).size() == 1, "grouped output access size is wrong");
    ok &= check(f.parameters().size() == 1 && f.parameters()[0] == h.param(), "multi-output parameter collection is wrong");
    ok &= check(out.size() == 5 && near(out[0], 4.0) && near(out[1], 5.0) && near(out[2], 0.5) && near(out[3], 6.0) && near(out[4], 7.0), "multi-output eval result is wrong");
    ok &= check(sp.rows() == 5 && sp.cols() == 3, "multi-output sparsity dimensions are wrong");
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2), "first output row sparsity is wrong");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 2), "second output row sparsity is wrong");
    ok &= check(sp.contains(2, 0), "parameter-scaled output sparsity is wrong");
    ok &= check(sp.contains(3, 1) && sp.contains(3, 2), "product output sparsity is wrong");
    ok &= check(!sp.contains(4, 0) && !sp.contains(4, 1) && !sp.contains(4, 2), "constant output row should be sparse-empty");
    ok &= check(compile_and_run("ad_multi_output", f, {1.0, 2.0, 3.0}, {0.5}), "generated C for multi-output function should match native eval");
    return ok;
}

bool test_call_outputs_single_boundary() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec U = ad::vec_variable("U", 1);
    ad::Function f({X, U}, std::vector<ad::Vec>{ad::sin(X), ad::Vec{X[0] * U[0]}});

    ad::Vec Y = ad::vec_variable("Y", 2);
    ad::Vec V = ad::vec_variable("V", 1);
    std::vector<ad::Vec> outputs = f.call_outputs({Y, V});
    ad::Function outer({Y, V}, ad::concat(outputs[0], outputs[1]));
    ad::GraphInfo info = outer.info().output_graph;

    bool ok = true;
    ok &= check(outputs.size() == 2 && outputs[0].size() == 2 && outputs[1].size() == 1, "call_outputs group sizes are wrong");
    ok &= check(count_unique_kind(info, ad::GraphNodeKind::FunctionCall) == 1, "call_outputs should slice one FunctionCall boundary");
    ok &= check(contains_kind(info, ad::GraphNodeKind::Slice), "call_outputs should expose grouped slices");
    ok &= check(outer.input_size() == 3 && outer.output_size() == 3, "outer grouped-call function layout is wrong");
    return ok;
}

bool test_transformed_grouping_and_codegen() {
    ad::Vec X = ad::vec_variable("X", 2);
    ad::Vec U = ad::vec_variable("U", 1);
    ad::Function f({X, U}, std::vector<ad::Vec>{
        ad::sigmoid(X),
        ad::Vec{X[0] * U[0], X[1] + U[0]},
    });

    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();

    bool ok = true;
    ok &= check(df.output_count() == f.output_count() && df.output(0).size() == 2 && df.output(1).size() == 2, "forward_function should preserve output grouping");
    ok &= check(rf.output_count() == f.input_count() && rf.output(0).size() == 2 && rf.output(1).size() == 1, "reverse_function should group adjoints by primal input groups");
    ok &= check(compile_and_run("ad_multi_output_jvp", df, {1.0, 2.0, 3.0, 0.1, -0.2, 0.3}), "generated C for grouped forward function should match native eval");
    ok &= check(compile_and_run("ad_multi_output_vjp", rf, {1.0, 2.0, 3.0, 0.5, -0.25, 1.5, 2.0}), "generated C for grouped reverse function should match native eval");
    return ok;
}

bool test_validation_across_outputs() {
    bool ok = true;
    ok &= check(throws_runtime([] {
        ad::Vec X = ad::vec_variable("X", 1);
        ad::Vec Z = ad::vec_variable("Z", 1);
        (void)ad::Function({X}, std::vector<ad::Vec>{X, Z});
    }), "undeclared variables in any output group should throw");
    ok &= check(throws_runtime([] {
        ad::Vec X = ad::vec_variable("X", 1);
        ad::Expr p = ad::parameter("p");
        (void)ad::Function({X}, ad::Params{}, std::vector<ad::Vec>{ad::Vec{X[0] + p}});
    }), "undeclared explicit parameters in any output group should throw");
    ok &= check(throws_runtime([] {
        ad::Vec X = ad::vec_variable("X", 1);
        ad::Vec invalid;
        (void)ad::Function({X}, std::vector<ad::Vec>{invalid});
    }), "invalid output group should throw");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_single_output_compatibility();
    ok &= test_multi_output_metadata_eval_sparsity();
    ok &= test_call_outputs_single_boundary();
    ok &= test_transformed_grouping_and_codegen();
    ok &= test_validation_across_outputs();
    if (!ok) {
        return 1;
    }
    std::cout << "AD multi-output function tests passed\n";
    return 0;
}
