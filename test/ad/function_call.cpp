// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>

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

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool test_basic_call_node() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, ad::sin(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Vec out = f.call({y});
    ad::GraphInfo info = ad::inspect(out);

    bool ok = true;
    ok &= check(out.size() == 2, "function call output size is wrong");
    ok &= check(info.kind == ad::GraphNodeKind::FunctionCall, "function call should create a FunctionCall node");
    ok &= check(info.children.size() == 1 && info.children[0].kind == ad::GraphNodeKind::VectorVariable, "FunctionCall inspection should show argument graph");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "FunctionCall inspection should not inline callee body");
    return ok;
}

bool test_argument_validation() {
    bool ok = true;
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Vec u = ad::vec_variable("u", 2);
        ad::Function f({x, u}, x + u);
        (void)f.call({x});
    }), "wrong function call argument count should throw");
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Function f({x}, ad::sin(x));
        ad::Vec y = ad::vec_variable("y", 3);
        (void)f.call({y});
    }), "wrong function call argument size should throw");
    ok &= check(throws_cleanly([] {
        ad::Function invalid;
        ad::Vec y = ad::vec_variable("y", 2);
        (void)invalid.call({y});
    }), "calling invalid function should throw");
    return ok;
}

bool test_call_inside_outer_function() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, ad::sin(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Function g({y}, f.call({y}));

    bool ok = true;
    ok &= check(g.input_vars().size() == 2, "outer function input variable count is wrong");
    ok &= check(g.input_vars()[0] == y.variables()[0] && g.input_vars()[1] == y.variables()[1], "outer function should depend on argument variables");
    ok &= check(g.input_vars()[0] != x.variables()[0], "callee local variables should not leak into outer function");
    ok &= check(g.info().output_graph.kind == ad::GraphNodeKind::FunctionCall, "outer function output should retain FunctionCall boundary");
    return ok;
}

bool test_call_argument_expression_and_parameters() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function f({x}, ad::sin(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Vec p = ad::vec_parameter("p", 2);
    ad::Function g({y}, f.call({y + p}));

    ad::Expr a = ad::parameter("a");
    ad::Function scaled({x}, a * ad::sin(x));
    ad::Function h({y}, scaled.call({y}));

    bool ok = true;
    ok &= check(g.input_vars().size() == 2 && g.input_vars()[0] == y.variables()[0], "argument expression should expose outer variable dependencies");
    ok &= check(g.parameters().size() == 2 && g.parameters()[0] == p.parameters()[0] && g.parameters()[1] == p.parameters()[1], "argument expression should expose parameter dependencies");
    ok &= check(h.parameters().size() == 1 && h.parameters()[0] == a.param(), "callee parameter dependency should propagate through FunctionCall");
    return ok;
}

bool test_call_argument_undeclared_variable_rejected() {
    bool ok = true;
    ok &= check(throws_cleanly([] {
        ad::Vec x = ad::vec_variable("x", 2);
        ad::Function f({x}, ad::sin(x));
        ad::Vec y = ad::vec_variable("y", 2);
        ad::Vec z = ad::vec_variable("z", 2);
        ad::Function g({y}, f.call({y + z}));
        (void)g;
    }), "outer function should reject undeclared variables inside call arguments");
    return ok;
}

bool test_nested_defect_kernel() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec z = ad::vec_variable("z", 4);
    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(2, 4, {1.0, 0.0, -1.0, 0.0,
                                0.0, 1.0, 0.0, -1.0});

    ad::Vec defect_expr = Dloc * z - h * rhs.call({xc, uc});
    ad::Function defect({z, xc, uc}, defect_expr);
    ad::GraphInfo graph = defect.info().output_graph;

    bool ok = true;
    ok &= check(defect.input_count() == 3, "defect input count is wrong");
    ok &= check(defect.parameters().size() == 1 && defect.parameters()[0] == h.param(), "defect parameter dependency is wrong");
    ok &= check(contains_kind(graph, ad::GraphNodeKind::DenseMatVec), "defect graph should retain DenseMatVec");
    ok &= check(contains_kind(graph, ad::GraphNodeKind::FunctionCall), "defect graph should retain FunctionCall");
    ok &= check(defect.input_vars().size() == 8, "callee local variables should not leak into defect inputs");
    return ok;
}

bool test_same_name_binding_and_lifetime() {
    ad::Vec x_local = ad::vec_variable("x", 2);
    ad::Function f({x_local}, ad::sin(x_local));

    ad::Vec x_global = ad::vec_variable("x", 2);
    ad::Function g({x_global}, f.call({x_global}));

    ad::Function copy = f;
    ad::Vec y = ad::vec_variable("y", 2);
    ad::Function h({y}, copy.call({y}));

    bool ok = true;
    ok &= check(g.input_vars().size() == 2 && g.input_vars()[0] == x_global.variables()[0], "same-name binding should use argument IDs");
    ok &= check(g.input_vars()[0] != x_local.variables()[0], "same-name local variable should not leak");
    ok &= check(h.valid() && h.output_size() == 2 && h.info().output_graph.kind == ad::GraphNodeKind::FunctionCall, "copied function should remain callable");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_basic_call_node();
    ok &= test_argument_validation();
    ok &= test_call_inside_outer_function();
    ok &= test_call_argument_expression_and_parameters();
    ok &= test_call_argument_undeclared_variable_rejected();
    ok &= test_nested_defect_kernel();
    ok &= test_same_name_binding_and_lifetime();
    return ok ? 0 : 1;
}
