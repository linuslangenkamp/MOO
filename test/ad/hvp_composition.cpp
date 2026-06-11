// SPDX-License-Identifier: LGPL-3.0-or-later
#include <ad.h>

#include <algorithm>
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

bool contains_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return true;
    }
    return std::any_of(info.children.begin(), info.children.end(), [&](const ad::GraphInfo &child) {
        return contains_kind(child, kind);
    });
}

const ad::GraphInfo *find_kind(const ad::GraphInfo &info, ad::GraphNodeKind kind) {
    if (info.kind == kind) {
        return &info;
    }
    for (const ad::GraphInfo &child : info.children) {
        if (const ad::GraphInfo *found = find_kind(child, kind)) {
            return found;
        }
    }
    return nullptr;
}

std::string read_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("failed to open " + path);
    }
    return std::string(std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>());
}

bool test_scalar_nonlinear_hvp_composition() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Expr f_expr = ad::sin(x[0]) + x[0] * x[1];
    ad::Function f({x}, ad::Vec{f_expr});

    ad::Vec lambda = ad::vec_variable("lambda", 1);
    ad::Vec seed = ad::vec_variable("seed", f.input_size());

    ad::Vec grad = f.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(f.input_vars(), seed);
    ad::GraphInfo info = ad::inspect(hvp);

    bool ok = true;
    ok &= check(grad.size() == f.input_size(), "scalar nonlinear grad size is wrong");
    ok &= check(hvp.size() == f.input_size(), "scalar nonlinear HVP size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "scalar nonlinear HVP should use transformed FunctionCall composition");
    ok &= check(!contains_kind(info, ad::GraphNodeKind::DenseMatVec), "scalar nonlinear HVP should not introduce dense Hessian matvec");
    return ok;
}

bool test_vector_lagrangian_hvp_composition() {
    ad::Vec x = ad::vec_variable("x", 3);
    ad::Function f({x}, ad::sigmoid(x));

    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Vec seed = ad::vec_variable("seed", f.input_size());

    ad::Vec grad = f.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(f.input_vars(), seed);
    ad::GraphInfo info = ad::inspect(hvp);

    bool ok = true;
    ok &= check(grad.size() == 3, "vector lagrangian grad size is wrong");
    ok &= check(hvp.size() == 3, "vector lagrangian HVP size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "vector HVP should use transformed FunctionCall composition");
    ok &= check(!contains_kind(info, ad::GraphNodeKind::DenseMatVec), "vector HVP should not introduce dense Hessian matvec");
    return ok;
}

bool test_linear_matvec_hvp_shape() {
    ad::Vec X = ad::vec_variable("X", 3);
    ad::DenseMatrix D(2, 3, {1.0, 0.0, -2.0,
                             0.0, 3.0, 0.5});
    ad::Function f({X}, D * X);

    ad::Vec lambda = ad::vec_variable("lambda", f.output_size());
    ad::Vec seed = ad::vec_variable("seed", f.input_size());

    ad::Vec grad = f.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(f.input_vars(), seed);
    ad::GraphInfo info = ad::inspect(hvp);

    bool ok = true;
    ok &= check(grad.size() == 3, "linear matvec grad size is wrong");
    ok &= check(hvp.size() == 3, "linear matvec HVP size is wrong");
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "linear matvec HVP should use transformed FunctionCall composition");
    ok &= check(!contains_kind(info, ad::GraphNodeKind::VectorElement), "linear matvec HVP should not construct a dense Hessian by scalar element extraction");
    return ok;
}

bool test_collocation_like_local_defect_hvp_shape() {
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
    ad::Vec lambda = ad::vec_variable("lambda", defect.output_size());
    ad::Vec seed = ad::vec_variable("seed", defect.input_size());

    ad::Vec grad = defect.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(defect.input_vars(), seed);
    ad::GraphInfo info = ad::inspect(hvp);

    bool ok = true;
    ok &= check(grad.size() == defect.input_size(), "defect grad size is wrong");
    ok &= check(hvp.size() == defect.input_size(), "defect HVP size is wrong");
    ok &= check(contains_kind(ad::inspect(grad), ad::GraphNodeKind::FunctionCall), "defect reverse should retain generic FunctionCall boundary");
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "defect HVP should retain generic transformed FunctionCall boundary");
    ok &= check(hvp.parameters().size() == 1 && hvp.parameters()[0] == h.param(), "defect HVP should not introduce parameter tangent/adjoint dependencies");
    return ok;
}

bool test_function_call_boundary_survives_hvp_composition() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Function outer({y}, local.call({y}));

    ad::Vec lambda = ad::vec_variable("lambda", outer.output_size());
    ad::Vec seed = ad::vec_variable("seed", outer.input_size());
    ad::Vec grad = outer.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(outer.input_vars(), seed);
    ad::GraphInfo grad_info = ad::inspect(grad);
    ad::GraphInfo hvp_info = ad::inspect(hvp);

    bool ok = true;
    ok &= check(contains_kind(grad_info, ad::GraphNodeKind::FunctionCall), "reverse through outer call should retain generic FunctionCall");
    ok &= check(contains_kind(hvp_info, ad::GraphNodeKind::FunctionCall), "forward over reverse should retain generic transformed FunctionCall");
    ok &= check(ad::count_nodes(hvp_info, ad::GraphNodeKind::VectorUnary) == 0, "HVP call boundary inspection should not inline callee body");
    return ok;
}

bool test_lambda_tangent_behavior() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function local({x}, ad::sigmoid(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Vec lambda = ad::vec_variable("lambda", 2);
    ad::Vec call = local.call({y});
    ad::Vec grad = call.reverse_diff(ad::vars(y), lambda);

    ad::Vec seed_no_lambda = ad::vec_variable("seed_no_lambda", y.size());
    ad::Vec hvp_no_lambda = grad.forward_diff(ad::vars(y), seed_no_lambda);
    const ad::GraphInfo *no_lambda_boundary = find_kind(ad::inspect(hvp_no_lambda), ad::GraphNodeKind::FunctionCall);
    ad::Vars no_lambda_deps = hvp_no_lambda.variables();

    ad::Vec seed_with_lambda = ad::vec_variable("seed_with_lambda", y.size() + lambda.size());
    ad::Vec hvp_with_lambda = grad.forward_diff(ad::vars(y, lambda), seed_with_lambda);
    const ad::GraphInfo *with_lambda_boundary = find_kind(ad::inspect(hvp_with_lambda), ad::GraphNodeKind::FunctionCall);
    ad::Vars with_lambda_deps = hvp_with_lambda.variables();

    bool ok = true;
    ok &= check(no_lambda_boundary != nullptr, "lambda-not-in-wrt HVP should retain generic FunctionCall");
    ok &= check(with_lambda_boundary != nullptr, "lambda-in-wrt HVP should retain generic FunctionCall");
    ok &= check(with_lambda_boundary && contains_kind(*with_lambda_boundary, ad::GraphNodeKind::Slice), "lambda tangent should be gathered from seed when lambda is in wrt");
    ok &= check(no_lambda_deps.size() == 6, "lambda-not-in-wrt HVP should depend on y, lambda, and only y seed components");
    ok &= check(with_lambda_deps.size() == 8, "lambda-in-wrt HVP should additionally depend on lambda seed components");
    ok &= check(hvp_no_lambda.size() == y.size(), "lambda-not-in-wrt HVP size is wrong");
    ok &= check(hvp_with_lambda.size() == y.size(), "lambda-in-wrt HVP size is wrong");
    return ok;
}

bool test_callee_parameter_survives_reverse_forward_boundary() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Expr h = ad::parameter("h");
    ad::Function local({x}, h * ad::sigmoid(x));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Function outer({y}, local.call({y}));

    ad::Vec lambda = ad::vec_variable("lambda", outer.output_size());
    ad::Vec seed = ad::vec_variable("seed", outer.input_size());
    ad::Vec grad = outer.reverse(lambda);
    ad::Vec hvp = grad.forward_diff(outer.input_vars(), seed);
    ad::GraphInfo info = ad::inspect(hvp);
    ad::Vars deps = hvp.variables();

    bool ok = true;
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "parameterized callee HVP should retain generic transformed FunctionCall");
    ok &= check(hvp.parameters().size() == 1 && hvp.parameters()[0] == h.param(), "callee parameter should propagate through transformed function-call boundary");
    ok &= check(std::none_of(deps.values().begin(), deps.values().end(), [&](const ad::Var &var) {
        return var == x.variables()[0] || var == x.variables()[1];
    }), "parameterized callee local variables should not leak through HVP boundary");
    return ok;
}

bool test_nested_function_call_hvp_boundary() {
    ad::Vec x = ad::vec_variable("x", 2);
    ad::Function inner({x}, ad::sigmoid(x));

    ad::Vec m = ad::vec_variable("m", 2);
    ad::Function middle({m}, inner.call({m}));

    ad::Vec y = ad::vec_variable("y", 2);
    ad::Function outer({y}, middle.call({y}));

    ad::Vec lambda = ad::vec_variable("lambda", outer.output_size());
    ad::Vec seed = ad::vec_variable("seed", outer.input_size());
    ad::Vec hvp = outer.reverse(lambda).forward_diff(outer.input_vars(), seed);
    ad::GraphInfo info = ad::inspect(hvp);
    ad::Vars deps = hvp.variables();

    bool ok = true;
    ok &= check(contains_kind(info, ad::GraphNodeKind::FunctionCall), "nested function HVP should retain generic transformed FunctionCall boundary");
    ok &= check(ad::count_nodes(info, ad::GraphNodeKind::VectorUnary) == 0, "nested function HVP inspection should not inline nested callee bodies");
    ok &= check(std::none_of(deps.values().begin(), deps.values().end(), [&](const ad::Var &var) {
        return var == x.variables()[0] || var == m.variables()[0];
    }), "nested callee local variables should not leak through HVP boundary");
    return ok;
}

bool test_no_public_hessian_api() {
    const std::string ad_header = read_file(std::string(MOO_AD_SOURCE_DIR) + "/ad.h");
    const std::string function_header = read_file(std::string(MOO_AD_SOURCE_DIR) + "/function.h");
    return check(ad_header.find("hvp") == std::string::npos &&
                     ad_header.find("Hessian") == std::string::npos &&
                     function_header.find("hvp") == std::string::npos &&
                     function_header.find("Hessian") == std::string::npos,
                 "Chunk 6 should not add a public HVP/Hessian API");
}

} // namespace

int main() {
    bool ok = true;
    ok &= test_scalar_nonlinear_hvp_composition();
    ok &= test_vector_lagrangian_hvp_composition();
    ok &= test_linear_matvec_hvp_shape();
    ok &= test_collocation_like_local_defect_hvp_shape();
    ok &= test_function_call_boundary_survives_hvp_composition();
    ok &= test_lambda_tangent_behavior();
    ok &= test_callee_parameter_survives_reverse_forward_boundary();
    ok &= test_nested_function_call_hvp_boundary();
    ok &= test_no_public_hessian_api();
    return ok ? 0 : 1;
}
