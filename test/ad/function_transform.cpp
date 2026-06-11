// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
#include <iostream>
#include <string>

namespace {

bool check(bool condition, const std::string &message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        return false;
    }
    return true;
}

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

bool test_forward_function_layout_body_and_cache() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Function df1 = rhs.forward_function();
    ad::Function df2 = rhs.forward_function();
    ad::FunctionInfo info1 = df1.info();
    ad::FunctionInfo info2 = df2.info();

    bool ok = true;
    ok &= check(df1.input_count() == 4, "forward_function should expose primal plus tangent input groups");
    ok &= check(df1.input_size() == 8, "forward_function flattened input size is wrong");
    ok &= check(df1.output_size() == rhs.output_size(), "forward_function output size is wrong");
    ok &= check(info1.inputs[0].size == 2 && info1.inputs[1].size == 2 && info1.inputs[2].size == 2 && info1.inputs[3].size == 2, "forward_function input group sizes are wrong");
    ok &= check(contains_kind(info1.output_graph, ad::GraphNodeKind::VectorUnary), "forward_function body should contain unary derivative structure");
    ok &= check(contains_kind(info1.output_graph, ad::GraphNodeKind::Slice), "forward_function body should slice tangent layout");
    ok &= check(info1.output_graph.id == info2.output_graph.id, "forward_function should reuse cached transformed output graph");
    ok &= check(info1.inputs[2].node_id == info2.inputs[2].node_id && info1.inputs[3].node_id == info2.inputs[3].node_id, "forward_function should reuse cached tangent groups");
    return ok;
}

bool test_reverse_function_layout_body_and_cache() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Function rf1 = rhs.reverse_function();
    ad::Function rf2 = rhs.reverse_function();
    ad::FunctionInfo info1 = rf1.info();
    ad::FunctionInfo info2 = rf2.info();

    bool ok = true;
    ok &= check(rf1.input_count() == 3, "reverse_function should expose primal input groups plus lambda");
    ok &= check(rf1.input_size() == 6, "reverse_function flattened input size is wrong");
    ok &= check(rf1.output_size() == rhs.input_size(), "reverse_function output size should match primal input size");
    ok &= check(info1.inputs[0].size == 2 && info1.inputs[1].size == 2 && info1.inputs[2].size == 2, "reverse_function input group sizes are wrong");
    ok &= check(contains_kind(info1.output_graph, ad::GraphNodeKind::VectorUnary), "reverse_function body should contain unary derivative structure");
    ok &= check(contains_kind(info1.output_graph, ad::GraphNodeKind::ScatterSlice), "reverse_function body should scatter into flattened input layout");
    ok &= check(info1.output_graph.id == info2.output_graph.id, "reverse_function should reuse cached transformed output graph");
    ok &= check(info1.inputs[2].node_id == info2.inputs[2].node_id, "reverse_function should reuse cached lambda group");
    return ok;
}

bool test_convenience_wrappers_call_transformed_functions() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));
    ad::Vec seed = ad::vec_variable("seed", rhs.input_size());
    ad::Vec lambda = ad::vec_variable("lambda", rhs.output_size());

    ad::Vec jvp = rhs.forward(seed);
    ad::Vec vjp = rhs.reverse(lambda);
    ad::GraphInfo jvp_info = ad::inspect(jvp);
    ad::GraphInfo vjp_info = ad::inspect(vjp);

    bool ok = true;
    ok &= check(jvp.size() == rhs.output_size(), "forward convenience output size is wrong");
    ok &= check(vjp.size() == rhs.input_size(), "reverse convenience output size is wrong");
    ok &= check(jvp_info.kind == ad::GraphNodeKind::FunctionCall, "forward convenience should be a generic FunctionCall");
    ok &= check(vjp_info.kind == ad::GraphNodeKind::FunctionCall, "reverse convenience should be a generic FunctionCall");
    ok &= check(ad::count_nodes(jvp_info, ad::GraphNodeKind::VectorUnary) == 0, "forward convenience should not inline transformed function body");
    ok &= check(ad::count_nodes(vjp_info, ad::GraphNodeKind::VectorUnary) == 0, "reverse convenience should not inline transformed function body");
    return ok;
}

bool test_hvp_transformed_function_is_ordinary_function() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function rhs({x}, ad::sigmoid(x));

    ad::Function grad = rhs.reverse_function();
    ad::Function hv = grad.forward_function();
    ad::FunctionInfo info = hv.info();

    bool ok = true;
    ok &= check(hv.output_size() == rhs.input_size(), "forward(reverse(function)) output size is wrong");
    ok &= check(hv.input_count() == 4, "forward(reverse(function)) should expose reverse inputs plus their tangent groups");
    ok &= check(contains_kind(info.output_graph, ad::GraphNodeKind::VectorUnary), "forward(reverse(function)) body should remain symbolic");
    ok &= check(!contains_kind(info.output_graph, ad::GraphNodeKind::FunctionCall), "simple local HVP transform should not require a call boundary");
    return ok;
}

bool test_nested_defect_transform_bodies_preserve_structure() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Vec u = ad::vec_variable("u", 2);
    ad::Function rhs({x, u}, ad::sigmoid(x + u));

    ad::Vec z = ad::vec_variable("z", 4);
    ad::Vec xc = ad::vec_variable("xc", 2);
    ad::Vec uc = ad::vec_variable("uc", 2);
    ad::Expr h = ad::parameter("h");
    ad::DenseMatrix Dloc(2, 4, {1.0, 0.0, -1.0, 0.0,
                                0.0, 1.0, 0.0, -1.0});

    ad::Function defect({z, xc, uc}, Dloc * z - h * rhs.call({xc, uc}));
    ad::Function df = defect.forward_function();
    ad::Function rf = defect.reverse_function();

    bool ok = true;
    ok &= check(contains_kind(df.info().output_graph, ad::GraphNodeKind::DenseMatVec), "defect forward_function body should preserve DenseMatVec");
    ok &= check(contains_kind(df.info().output_graph, ad::GraphNodeKind::FunctionCall), "defect forward_function body should preserve transformed RHS call");
    ok &= check(df.parameters().size() == 1 && df.parameters()[0] == h.param(), "defect forward_function should retain parameter dependency");
    ok &= check(contains_kind(rf.info().output_graph, ad::GraphNodeKind::DenseMatVec), "defect reverse_function body should preserve transpose DenseMatVec");
    ok &= check(contains_kind(rf.info().output_graph, ad::GraphNodeKind::FunctionCall), "defect reverse_function body should preserve transformed RHS call");
    ok &= check(rf.parameters().size() == 1 && rf.parameters()[0] == h.param(), "defect reverse_function should retain parameter dependency");
    return ok;
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_forward_function_layout_body_and_cache();
    ok &= test_reverse_function_layout_body_and_cache();
    ok &= test_convenience_wrappers_call_transformed_functions();
    ok &= test_hvp_transformed_function_is_ordinary_function();
    ok &= test_nested_defect_transform_bodies_preserve_structure();
    return ok ? 0 : 1;
}
