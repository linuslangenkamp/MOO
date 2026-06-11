// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
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

template <typename Fn>
bool throws_cleanly(Fn &&fn) {
    try {
        fn();
    } catch (const std::exception &) {
        return true;
    }
    return false;
}

std::string read_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open " + path);
    }
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

std::string lower_ascii(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char ch) {
        return static_cast<char>(std::tolower(ch));
    });
    return value;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool test_simple_function() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, x + 1.0);
    ad::FunctionInfo info = f.info();

    bool ok = true;
    ok &= check(f.valid(), "simple function should be valid");
    ok &= check(f.input_count() == 1, "simple function input count is wrong");
    ok &= check(f.input_size() == 2, "simple function input size is wrong");
    ok &= check(f.output_size() == 2, "simple function output size is wrong");
    ok &= check(f.input_vars().size() == 2, "simple function flattened input vars are wrong");
    ok &= check(f.parameters().empty(), "simple function should not collect parameters");
    ok &= check(info.inputs.size() == 1 && info.inputs[0].label == "x" && info.inputs[0].size == 2, "simple function input info is wrong");
    ok &= check(info.output_graph.kind == ad::GraphNodeKind::VectorBinary, "simple function output graph kind is wrong");
    return ok;
}

bool test_two_input_groups_order() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 1);
    ad::Function f({x, u}, ad::concat(x, u + 2.0));
    ad::Vars vars = f.input_vars();

    bool ok = true;
    ok &= check(f.input_count() == 2, "two-input function count is wrong");
    ok &= check(f.input_size() == 3, "two-input flattened size is wrong");
    ok &= check(vars.size() == 3, "two-input flattened vars size is wrong");
    ok &= check(vars[0] == x.variables()[0] && vars[1] == x.variables()[1] && vars[2] == u.variables()[0], "flattened input variable order should preserve {x, u}");
    ok &= check(f.info().inputs[0].label == "x" && f.info().inputs[1].label == "u", "input labels should preserve group order");
    return ok;
}

bool test_reject_invalid_inputs_and_outputs() {
    bool ok = true;
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Function f({x + 1.0}, x);
        (void)f;
    }), "non-variable input group should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Function f({x, x}, x);
        (void)f;
    }), "duplicate input group should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec z = ad::vec_variable("z", 0);
        ad::Function f({z}, ad::vec_constant({}));
        (void)f;
    }), "zero-size input group should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Vec y = ad::vec_variable("y", 2);
        ad::Function f({x}, x + y);
        (void)f;
    }), "undeclared output variable should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec X1 = ad::vec_variable("X", 2);
        ad::Vec X2 = ad::vec_variable("X", 2);
        ad::Function f({X1}, X2);
        (void)f;
    }), "same-name but different-ID variable should throw when undeclared");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 1);
        ad::Function f({}, x);
        (void)f;
    }), "empty input list should reject variable-dependent output");
    return ok;
}

bool test_parameter_collection() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Expr h = ad::parameter("h");
    ad::Function scalar_param({x}, h * x);

    ad::Vec p = ad::vec_parameter("p", 2);
    ad::Function vector_param({x}, x + p);

    bool ok = true;
    ok &= check(scalar_param.parameters().size() == 1, "scalar parameter dependency count is wrong");
    ok &= check(scalar_param.parameters()[0] == h.param(), "scalar parameter dependency ID is wrong");
    ok &= check(vector_param.parameters().size() == 2, "vector parameter dependency count should be flattened components");
    ok &= check(vector_param.parameters()[0] == p.parameters()[0] && vector_param.parameters()[1] == p.parameters()[1], "vector parameter dependency IDs are wrong");
    ok &= check(vector_param.info().parameter_count == 2, "FunctionInfo parameter count is wrong");
    return ok;
}

bool test_explicit_parameter_layout() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Expr a = ad::parameter("a");
    ad::Expr b = ad::parameter("b");
    ad::Expr unused = ad::parameter("unused");

    ad::Function f({x}, ad::params(b, unused, a), b * x + a * x);
    ad::Function df = f.forward_function();
    ad::Function rf = f.reverse_function();

    bool ok = true;
    ok &= check(f.parameters().size() == 3, "explicit function parameter count is wrong");
    ok &= check(f.parameters()[0] == b.param() && f.parameters()[1] == unused.param() && f.parameters()[2] == a.param(), "explicit function parameter order should be preserved");
    ok &= check(df.parameters().size() == 3 && df.parameters()[0] == b.param() && df.parameters()[1] == unused.param() && df.parameters()[2] == a.param(), "forward_function should preserve explicit parameter layout");
    ok &= check(rf.parameters().size() == 3 && rf.parameters()[0] == b.param() && rf.parameters()[1] == unused.param() && rf.parameters()[2] == a.param(), "reverse_function should preserve explicit parameter layout");
    ok &= check(throws_cleanly([&] {
        ad::Function missing({x}, ad::params(a), b * x);
        (void)missing;
    }), "explicit parameter layout should reject undeclared output parameter dependencies");
    ok &= check(throws_cleanly([&] {
        ad::Function duplicate({x}, ad::Params({a.param(), a.param()}), a * x);
        (void)duplicate;
    }), "explicit parameter layout should reject duplicate parameter IDs");
    return ok;
}

bool test_structured_output_traversal() {
    ad::Vec x = ad::vec_variable("x", 3);
    ad::Vec p = ad::vec_parameter("p", 2);
    ad::DenseMatrix A(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});

    ad::Function f({x}, ad::sigmoid(A * x + p));
    ad::GraphInfo graph = f.info().output_graph;

    bool ok = true;
    ok &= check(f.parameters().size() == 2, "structured output should collect vector parameter components");
    ok &= check(contains_kind(graph, ad::GraphNodeKind::DenseMatVec), "structured output should retain DenseMatVec");
    ok &= check(contains_kind(graph, ad::GraphNodeKind::VectorUnary), "structured output should retain VectorUnary");
    ok &= check(ad::count_nodes(graph, ad::GraphNodeKind::VectorElement) == 0, "structured output traversal should not scalar-lower");
    return ok;
}

bool test_local_defect_like_function() {
    ad::Vec z = ad::vec_variable("z", 4);
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(2, 4, {1.0, 0.0, -1.0, 0.0,
                                0.0, 1.0, 0.0, -1.0});

    ad::Vec rhs = ad::sigmoid(x + u);
    ad::Vec defect = Dloc * z - h * rhs;
    ad::Function f({z, x, u}, defect);

    bool ok = true;
    ok &= check(f.input_count() == 3, "defect function input count is wrong");
    ok &= check(f.input_size() == 8, "defect function flattened input size is wrong");
    ok &= check(f.output_size() == 2, "defect function output size is wrong");
    ok &= check(f.parameters().size() == 1 && f.parameters()[0] == h.param(), "defect function parameter dependency is wrong");
    ok &= check(contains_kind(f.info().output_graph, ad::GraphNodeKind::DenseMatVec), "defect function should retain DenseMatVec");
    ok &= check(contains_kind(f.info().output_graph, ad::GraphNodeKind::VectorScale), "defect function should retain scalar/vector scale");
    return ok;
}

bool test_copy_and_empty_decisions() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, x + 1.0);
    ad::Function copy = f;

    ad::Expr h = ad::parameter("h");
    ad::Function parameter_only({}, ad::Vec{h});
    ad::Function constant_empty({}, ad::vec_constant({}));

    bool ok = true;
    ok &= check(copy.valid() && copy.input_size() == 2 && copy.output_size() == 2, "copied function metadata should remain valid");
    ok &= check(parameter_only.input_count() == 0 && parameter_only.input_size() == 0, "parameter-only function should allow empty input list");
    ok &= check(parameter_only.parameters().size() == 1 && parameter_only.parameters()[0] == h.param(), "parameter-only function should collect scalar parameter");
    ok &= check(constant_empty.input_count() == 0 && constant_empty.output_size() == 0, "constant empty output function should be allowed");
    return ok;
}

bool test_public_headers_do_not_expose_deferred_concepts() {
    const std::vector<std::string> headers{
        "ad.h",
        "codegen_options.h",
        "expr.h",
        "function.h",
        "graph_info.h",
        "matrix.h",
        "map_kind.h",
        "symbol.h",
        "vec.h",
    };
    const std::vector<std::string> forbidden{
        "builder",
        "plan",
        "emitter",
        "callback",
        "vm",
        "python",
    };

    bool ok = true;
    for (const std::string &header : headers) {
        const std::string body = lower_ascii(read_file(std::string(MOO_AD_SOURCE_DIR) + "/" + header));
        for (const std::string &word : forbidden) {
            ok &= check(body.find(word) == std::string::npos, "forbidden public concept found in " + header + ": " + word);
        }
    }
    return ok;
}

bool test_codegen_options_are_metadata_only() {
    ad::CodegenOptions options;
    options.function_name = "rhs_eval";
    options.scalar_type = "double";
    options.emit_sparsity_metadata = true;

    bool ok = true;
    ok &= check(options.function_name == "rhs_eval", "CodegenOptions function name metadata is wrong");
    ok &= check(options.scalar_type == "double", "CodegenOptions scalar type metadata is wrong");
    ok &= check(options.emit_sparsity_metadata, "CodegenOptions sparsity metadata flag is wrong");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_simple_function();
    ok &= test_two_input_groups_order();
    ok &= test_reject_invalid_inputs_and_outputs();
    ok &= test_parameter_collection();
    ok &= test_explicit_parameter_layout();
    ok &= test_structured_output_traversal();
    ok &= test_local_defect_like_function();
    ok &= test_copy_and_empty_decisions();
    ok &= test_public_headers_do_not_expose_deferred_concepts();
    ok &= test_codegen_options_are_metadata_only();
    return ok ? 0 : 1;
}
