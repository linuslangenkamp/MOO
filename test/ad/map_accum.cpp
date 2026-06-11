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

bool near(double lhs, double rhs, double tol = 1e-9) {
    return std::fabs(lhs - rhs) <= tol;
}

bool throws_runtime(const std::function<void()> &fn) {
    try {
        fn();
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
    function.eval(input.empty() ? nullptr : input.data(), static_cast<int>(input.size()),
                  params.empty() ? nullptr : params.data(), static_cast<int>(params.size()),
                  work,
                  out.empty() ? nullptr : out.data(), static_cast<int>(out.size()));
    return out;
}

bool compile_and_run(const std::string &name,
                     const ad::Function &function,
                     const std::vector<double> &input) {
    ad::CodegenOptions options;
    options.function_name = name;
    const ad::CodegenResult generated = ad::generate_c(function, options);
    const std::vector<double> expected = eval_function(function, input);
    const bool has_loop = generated.source.find("for (int rep = 0; rep < 3; ++rep)") != std::string::npos;
    if (!has_loop) {
        std::cerr << "FAIL: generated map-accum code is missing recurrence loop\n";
        return false;
    }

    const std::filesystem::path dir = std::filesystem::path(MOO_BINARY_DIR) / "ad_mapaccum_tmp";
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

ad::Function make_step(ad::Vec &carry, ad::Vec &u) {
    carry = ad::vec_variable("carry", 2);
    u = ad::vec_variable("u", 2);
    ad::Vec next = carry + u;
    ad::Vec emit = carry * u;
    return ad::Function({carry, u}, std::vector<ad::Vec>{next, emit});
}

bool test_construction_eval_and_sparsity() {
    ad::Vec carry;
    ad::Vec u;
    ad::Function step = make_step(carry, u);
    ad::Vec c0 = ad::vec_variable("c0", 2);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::MapAccumResult result = ad::map_accum(step, carry, c0, 3, {
        ad::bind(u, ad::Map::blocks(U, 3, 2)),
    });

    ad::Vec full = ad::concat(result.final_carry, ad::concat(result.carry_trajectory, result.outputs[0]));
    ad::Function f({c0, U}, full);
    const std::vector<double> out = eval_function(f, {1.0, 10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    const std::vector<double> expected = {
        13.0, 25.0,
        3.0, 13.0, 7.0, 18.0, 13.0, 25.0,
        2.0, 30.0, 12.0, 65.0, 42.0, 126.0,
    };

    bool ok = true;
    ok &= check(out.size() == expected.size(), "map_accum output size is wrong");
    for (std::size_t i = 0; i < expected.size(); ++i) {
        ok &= check(near(out[i], expected[i]), "map_accum eval mismatch at row " + std::to_string(i));
    }
    ok &= check(contains_kind(ad::inspect(full), ad::GraphNodeKind::MapAccumCall), "map_accum graph should contain MapAccumCall");
    ok &= check(ad::count_nodes(ad::inspect(full), ad::GraphNodeKind::VectorBinary) == 0, "map_accum inspection should not inline step body");

    const ad::SparsityPattern sp = f.jacobian_sparsity();
    ok &= check(sp.contains(0, 0) && sp.contains(0, 2) && sp.contains(0, 4) && sp.contains(0, 6), "final carry row 0 sparsity is wrong");
    ok &= check(sp.contains(1, 1) && sp.contains(1, 3) && sp.contains(1, 5) && sp.contains(1, 7), "final carry row 1 sparsity is wrong");
    ok &= check(!sp.contains(0, 1) && !sp.contains(1, 0), "final carry sparsity should preserve components");
    return ok;
}

bool test_fold_scan_validation_and_derivatives() {
    ad::Vec carry;
    ad::Vec u;
    ad::Function step = make_step(carry, u);
    ad::Vec c0 = ad::vec_variable("c0", 2);
    ad::Vec U = ad::vec_variable("U", 6);

    ad::Vec folded = ad::fold(step, carry, c0, 3, {ad::bind(u, ad::Map::blocks(U, 3, 2))});
    ad::Vec scanned = ad::scan(step, carry, c0, 3, {ad::bind(u, ad::Map::blocks(U, 3, 2))});
    ad::Function f({c0, U}, folded);
    ad::Function scan_fun({c0, U}, scanned);

    std::vector<double> final = eval_function(f, {1.0, 10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
    std::vector<double> scan = eval_function(scan_fun, {1.0, 10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});

    ad::Vec seed = ad::vec_variable("seed", f.input_size());
    ad::Function jvp_fun({c0, U, seed}, f.forward(seed));
    std::vector<double> jvp = eval_function(jvp_fun, {1.0, 10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
                                                      1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});

    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Function vjp_fun({c0, U, lambda}, f.reverse(lambda));
    std::vector<double> vjp = eval_function(vjp_fun, {1.0, 10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 2.0, 3.0});

    bool ok = true;
    ok &= check(final.size() == 2 && near(final[0], 13.0) && near(final[1], 25.0), "fold final carry eval is wrong");
    ok &= check(scan.size() == 6 && near(scan[4], 13.0) && near(scan[5], 25.0), "scan eval is wrong");
    ok &= check(jvp.size() == 2 && near(jvp[0], 4.0) && near(jvp[1], 4.0), "fold JVP eval is wrong");
    const std::vector<double> expected_vjp = {2.0, 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 3.0};
    for (std::size_t i = 0; i < expected_vjp.size(); ++i) {
        ok &= check(near(vjp[i], expected_vjp[i]), "fold VJP eval mismatch at row " + std::to_string(i));
    }

    ok &= check(throws_runtime([&] {
        (void)ad::fold(step, carry, c0, 0, {ad::bind(u, ad::Map::blocks(U, 3, 2))});
    }), "zero reps should throw");
    ok &= check(throws_runtime([&] {
        (void)ad::fold(step, carry, c0, 3, {});
    }), "missing sequence binding should throw");
    ok &= check(throws_runtime([&] {
        (void)ad::fold(step, carry, c0, 3, {ad::bind(carry, ad::Map::blocks(U, 3, 2))});
    }), "carry sequence binding should throw");
    return ok;
}

bool test_duplicate_reverse_accumulation_and_hvp_shape() {
    ad::Vec carry = ad::vec_variable("carry", 1);
    ad::Vec u = ad::vec_variable("u", 1);
    ad::Function step({carry, u}, carry + u);
    ad::Vec c0 = ad::vec_variable("c0", 1);
    ad::Vec U = ad::vec_variable("U", 1);
    ad::Vec folded = ad::fold(step, carry, c0, 3, {
        ad::bind(u, ad::Map::stride(U, 3, 1, 0, 0, 1)),
    });
    ad::Function f({c0, U}, folded);

    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Function vjp_fun({c0, U, lambda}, f.reverse(lambda));
    std::vector<double> vjp = eval_function(vjp_fun, {5.0, 7.0, 2.0});

    ad::Vec seed = ad::vec_variable("seed", f.input_size());
    ad::Vec hvp = f.reverse(lambda).forward_diff(f.input_vars(), seed);

    bool ok = true;
    ok &= check(vjp.size() == 2 && near(vjp[0], 2.0) && near(vjp[1], 6.0), "duplicate reverse accumulation through repeated map index is wrong");
    ok &= check(hvp.size() == f.input_size(), "map-accum HVP composition has wrong size");
    ok &= check(ad::sparsity(hvp, ad::vars(c0, U)).rows() == f.input_size(), "map-accum HVP sparsity construction failed");
    return ok;
}

bool test_codegen() {
    ad::Vec carry;
    ad::Vec u;
    ad::Function step = make_step(carry, u);
    ad::Vec c0 = ad::vec_variable("c0", 2);
    ad::Vec U = ad::vec_variable("U", 6);
    ad::Vec folded = ad::fold(step, carry, c0, 3, {ad::bind(u, ad::Map::blocks(U, 3, 2))});
    ad::Function f({c0, U}, folded);
    return compile_and_run("mapaccum_codegen", f, {1.0, 10.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0});
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_construction_eval_and_sparsity();
    ok &= test_fold_scan_validation_and_derivatives();
    ok &= test_duplicate_reverse_accumulation_and_hvp_shape();
    ok &= test_codegen();
    if (!ok) {
        return 1;
    }
    std::cout << "map_accum tests passed\n";
    return 0;
}
