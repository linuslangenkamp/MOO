// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2026 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <ad.h>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool near(double a, double b) {
    return std::abs(a - b) < 1e-12;
}

} // namespace

int main() {
    using namespace ad;

    GraphFunctionBuilder kb;
    auto x = kb.inputs("x", 4);
    auto u = kb.params("u", 1);

    auto first = kb.slice(x, 0, 2);
    auto second = kb.slice(x, 2, 2);
    DenseMatrix D(2, 2, {1.0, 2.0, 3.0, 4.0});
    auto linear = kb.dense_matvec(D, first);
    auto nonlinear = kb.vector({kb.at(second, 0) * kb.at(second, 0), kb.at(second, 1) + kb.at(u, 0)});
    auto residual = kb.add(linear, kb.scale(0.5, nonlinear));

    auto F = kb.function(x, residual, {u});
    if (!F.has_vector_structure() || F.vector_node_count() < 6) {
        std::cerr << "vector expression structure was not preserved in GraphFunction\n";
        return 1;
    }
    bool saw_dense_matvec = false;
    for (const auto &node : F.vector_nodes) {
        saw_dense_matvec = saw_dense_matvec || node.op == VectorOp::DenseMatVec;
    }
    if (!saw_dense_matvec) {
        std::cerr << "dense matvec vector node was not preserved in GraphFunction\n";
        return 1;
    }

    double xv[4] = {1.0, 2.0, 3.0, 4.0};
    double uv[1] = {5.0};
    double out[2] = {0.0, 0.0};
    VM(F).evaluate(EvalEnv().input("x", xv).param("u", uv), out);
    if (!near(out[0], 9.5) || !near(out[1], 15.5)) {
        std::cerr << "unexpected vector expression value: " << out[0] << ", " << out[1] << "\n";
        return 1;
    }

    auto JVP = forward_diff(F, "x", "v");
    if (!JVP.has_vector_structure()) {
        std::cerr << "forward-mode transform did not preserve vector structure\n";
        return 1;
    }
    bool jvp_saw_dense_matvec = false;
    for (const auto &node : JVP.vector_nodes) {
        jvp_saw_dense_matvec = jvp_saw_dense_matvec || node.op == VectorOp::DenseMatVec;
    }
    if (!jvp_saw_dense_matvec) {
        std::cerr << "forward-mode transform did not preserve dense matvec vector node\n";
        return 1;
    }
    double vv[4] = {0.25, -0.5, 1.5, 2.0};
    double jout[2] = {0.0, 0.0};
    VM(JVP).evaluate(EvalEnv().input("x", xv).param("u", uv).param("v", vv), jout);
    if (!near(jout[0], 3.75) || !near(jout[1], -0.25)) {
        std::cerr << "unexpected vector expression JVP: " << jout[0] << ", " << jout[1] << "\n";
        return 1;
    }

    auto c = to_c(F, "vector_expr_value");
    if (c.find("void vector_expr_value") == std::string::npos) {
        std::cerr << "generated C did not contain expected function name\n";
        return 1;
    }
    if (c.find("static const double mat") == std::string::npos || c.find("moo_ad_dense_matvec") == std::string::npos) {
        std::cerr << "generated C did not use vector-aware dense matvec lowering\n";
        return 1;
    }
    if (c.find("double vec0") != std::string::npos) {
        std::cerr << "generated C unnecessarily materialized the input vector in vector-aware lowering\n";
        return 1;
    }

    auto jvp_c = to_c(JVP, "vector_expr_jvp");
    if (jvp_c.find("void vector_expr_jvp") == std::string::npos || jvp_c.find("static const double mat") == std::string::npos ||
        jvp_c.find("moo_ad_dense_matvec") == std::string::npos) {
        std::cerr << "generated JVP C did not use vector-aware dense matvec lowering\n";
        return 1;
    }

    auto MixedVJP = reverse_diff(F, "lambda", "x");
    if (!MixedVJP.has_vector_structure()) {
        std::cerr << "mixed vector VJP did not preserve vector structure\n";
        return 1;
    }
    double mixed_lambda[2] = {5.0, 7.0};
    double mixed_grad[4] = {0.0, 0.0, 0.0, 0.0};
    VM(MixedVJP).evaluate(EvalEnv().input("x", xv).param("u", uv).param("lambda", mixed_lambda), mixed_grad);
    if (!near(mixed_grad[0], 26.0) || !near(mixed_grad[1], 38.0) || !near(mixed_grad[2], 15.0) || !near(mixed_grad[3], 3.5)) {
        std::cerr << "unexpected mixed vector VJP: " << mixed_grad[0] << ", " << mixed_grad[1] << ", " << mixed_grad[2] << ", " << mixed_grad[3] << "\n";
        return 1;
    }
    auto MixedHVP = forward_diff(MixedVJP, "x", "v");
    if (!MixedHVP.has_vector_structure()) {
        std::cerr << "mixed vector HVP did not preserve vector structure\n";
        return 1;
    }
    double mixed_hvp[4] = {0.0, 0.0, 0.0, 0.0};
    VM(MixedHVP).evaluate(EvalEnv().input("x", xv).param("u", uv).param("lambda", mixed_lambda).param("v", vv), mixed_hvp);
    if (!near(mixed_hvp[0], 0.0) || !near(mixed_hvp[1], 0.0) || !near(mixed_hvp[2], 7.5) || !near(mixed_hvp[3], 0.0)) {
        std::cerr << "unexpected mixed vector HVP: " << mixed_hvp[0] << ", " << mixed_hvp[1] << ", " << mixed_hvp[2] << ", " << mixed_hvp[3] << "\n";
        return 1;
    }

    GraphFunctionBuilder elem_builder;
    auto ex = elem_builder.inputs("x", 2);
    auto ey = elem_builder.params("y", 2);
    auto elem = elem_builder.add(elem_builder.mul(ex, ey), elem_builder.scale(0.5, elem_builder.pow_const(ex, 2.0)));
    auto EF = elem_builder.function(ex, elem, {ey});
    bool saw_mul = false;
    bool saw_pow = false;
    for (const auto &node : EF.vector_nodes) {
        saw_mul = saw_mul || node.op == VectorOp::Mul;
        saw_pow = saw_pow || node.op == VectorOp::PowConst;
    }
    if (!saw_mul || !saw_pow) {
        std::cerr << "elementwise vector nodes were not preserved\n";
        return 1;
    }
    double exv[2] = {2.0, 3.0};
    double eyv[2] = {5.0, 7.0};
    double eout[2] = {0.0, 0.0};
    VM(EF).evaluate(EvalEnv().input("x", exv).param("y", eyv), eout);
    if (!near(eout[0], 12.0) || !near(eout[1], 25.5)) {
        std::cerr << "unexpected elementwise vector value: " << eout[0] << ", " << eout[1] << "\n";
        return 1;
    }
    auto EJVP = forward_diff(EF, "x", "v");
    bool jvp_saw_mul = false;
    bool jvp_saw_pow = false;
    for (const auto &node : EJVP.vector_nodes) {
        jvp_saw_mul = jvp_saw_mul || node.op == VectorOp::Mul;
        jvp_saw_pow = jvp_saw_pow || node.op == VectorOp::PowConst;
    }
    if (!jvp_saw_mul || !jvp_saw_pow) {
        std::cerr << "elementwise vector JVP did not preserve nonlinear vector nodes\n";
        return 1;
    }
    double ev[2] = {0.25, -0.5};
    double ejvp[2] = {0.0, 0.0};
    VM(EJVP).evaluate(EvalEnv().input("x", exv).param("y", eyv).param("v", ev), ejvp);
    if (!near(ejvp[0], 1.75) || !near(ejvp[1], -5.0)) {
        std::cerr << "unexpected elementwise vector JVP: " << ejvp[0] << ", " << ejvp[1] << "\n";
        return 1;
    }
    auto EVJP = reverse_diff(EF, "lambda", "x");
    bool vjp_saw_mul = false;
    bool vjp_saw_pow = false;
    for (const auto &node : EVJP.vector_nodes) {
        vjp_saw_mul = vjp_saw_mul || node.op == VectorOp::Mul;
        vjp_saw_pow = vjp_saw_pow || node.op == VectorOp::PowConst;
    }
    if (!vjp_saw_mul || !vjp_saw_pow) {
        std::cerr << "elementwise vector VJP did not preserve nonlinear vector nodes\n";
        return 1;
    }
    double elambda[2] = {11.0, 13.0};
    double egrad[2] = {0.0, 0.0};
    VM(EVJP).evaluate(EvalEnv().input("x", exv).param("y", eyv).param("lambda", elambda), egrad);
    if (!near(egrad[0], 77.0) || !near(egrad[1], 130.0)) {
        std::cerr << "unexpected elementwise vector VJP: " << egrad[0] << ", " << egrad[1] << "\n";
        return 1;
    }
    auto EHVP = forward_diff(EVJP, "x", "v");
    double ehvp[2] = {0.0, 0.0};
    VM(EHVP).evaluate(EvalEnv().input("x", exv).param("y", eyv).param("lambda", elambda).param("v", ev), ehvp);
    if (!near(ehvp[0], 2.75) || !near(ehvp[1], -6.5)) {
        std::cerr << "unexpected elementwise vector HVP: " << ehvp[0] << ", " << ehvp[1] << "\n";
        return 1;
    }
    auto elem_c = to_c(EF, "elementwise_vector_value");
    if (elem_c.find("moo_ad_vec_mul") == std::string::npos || elem_c.find("moo_ad_vec_pow_const") == std::string::npos) {
        std::cerr << "generated elementwise vector C did not use vector helpers\n";
        return 1;
    }

    GraphFunctionBuilder unary_builder;
    auto ux = unary_builder.inputs("x", 2);
    auto unary_out = unary_builder.add(
        unary_builder.add(unary_builder.unary(Op::Sin, ux), unary_builder.unary(Op::Cos, ux)),
        unary_builder.add(unary_builder.add(unary_builder.unary(Op::Tan, ux), unary_builder.unary(Op::Exp, ux)), unary_builder.unary(Op::Log, ux)));
    auto UF = unary_builder.function(ux, unary_out);
    bool saw_unary = false;
    for (const auto &node : UF.vector_nodes) {
        saw_unary = saw_unary || node.op == VectorOp::Unary;
    }
    if (!saw_unary) {
        std::cerr << "vector unary node was not preserved\n";
        return 1;
    }
    double uxv[2] = {1.2, 1.7};
    double uout[2] = {0.0, 0.0};
    VM(UF).evaluate(EvalEnv().input("x", uxv), uout);
    auto unary_value = [](double x) { return std::sin(x) + std::cos(x) + std::tan(x) + std::exp(x) + std::log(x); };
    if (!near(uout[0], unary_value(uxv[0])) || !near(uout[1], unary_value(uxv[1]))) {
        std::cerr << "unexpected vector unary value: " << uout[0] << ", " << uout[1] << "\n";
        return 1;
    }
    auto UJVP = forward_diff(UF, "x", "v");
    bool unary_jvp_saw_unary = false;
    bool unary_jvp_saw_div = false;
    for (const auto &node : UJVP.vector_nodes) {
        unary_jvp_saw_unary = unary_jvp_saw_unary || node.op == VectorOp::Unary;
        unary_jvp_saw_div = unary_jvp_saw_div || node.op == VectorOp::Div;
    }
    if (!unary_jvp_saw_unary || !unary_jvp_saw_div) {
        std::cerr << "vector unary JVP did not preserve unary/div vector nodes\n";
        return 1;
    }
    double uvv[2] = {0.25, -0.5};
    double ujvp[2] = {0.0, 0.0};
    VM(UJVP).evaluate(EvalEnv().input("x", uxv).param("v", uvv), ujvp);
    auto unary_derivative = [](double x) { return std::cos(x) - std::sin(x) + 1.0 + std::tan(x) * std::tan(x) + std::exp(x) + 1.0 / x; };
    if (!near(ujvp[0], unary_derivative(uxv[0]) * uvv[0]) || !near(ujvp[1], unary_derivative(uxv[1]) * uvv[1])) {
        std::cerr << "unexpected vector unary JVP: " << ujvp[0] << ", " << ujvp[1] << "\n";
        return 1;
    }
    auto UVJP = reverse_diff(UF, "lambda", "x");
    double ulambda[2] = {3.0, 5.0};
    double ugrad[2] = {0.0, 0.0};
    VM(UVJP).evaluate(EvalEnv().input("x", uxv).param("lambda", ulambda), ugrad);
    if (!near(ugrad[0], ulambda[0] * unary_derivative(uxv[0])) || !near(ugrad[1], ulambda[1] * unary_derivative(uxv[1]))) {
        std::cerr << "unexpected vector unary VJP: " << ugrad[0] << ", " << ugrad[1] << "\n";
        return 1;
    }
    auto UHVP = forward_diff(UVJP, "x", "v");
    double uhvp[2] = {0.0, 0.0};
    VM(UHVP).evaluate(EvalEnv().input("x", uxv).param("lambda", ulambda).param("v", uvv), uhvp);
    auto unary_second = [](double x) {
        double tan_x = std::tan(x);
        return -std::sin(x) - std::cos(x) + 2.0 * tan_x * (1.0 + tan_x * tan_x) + std::exp(x) - 1.0 / (x * x);
    };
    if (!near(uhvp[0], ulambda[0] * unary_second(uxv[0]) * uvv[0]) || !near(uhvp[1], ulambda[1] * unary_second(uxv[1]) * uvv[1])) {
        std::cerr << "unexpected vector unary HVP: " << uhvp[0] << ", " << uhvp[1] << "\n";
        return 1;
    }
    auto unary_c = to_c(UF, "unary_vector_value");
    if (unary_c.find("moo_ad_vec_sin") == std::string::npos || unary_c.find("moo_ad_vec_cos") == std::string::npos ||
        unary_c.find("moo_ad_vec_tan") == std::string::npos || unary_c.find("moo_ad_vec_exp") == std::string::npos ||
        unary_c.find("moo_ad_vec_log") == std::string::npos) {
        std::cerr << "generated vector unary C did not use vector helpers\n";
        return 1;
    }
    auto ujac_pairs = jacobian_sparsity(UF, "x").to_pairs();
    auto ujac_fn = sparse_jacobian_function(UF, "x", ujac_pairs);
    std::vector<double> ujac_values(ujac_pairs.size(), 0.0);
    VM(ujac_fn).evaluate(EvalEnv().input("x", uxv), ujac_values.data());
    for (std::size_t i = 0; i < ujac_pairs.size(); ++i) {
        const auto [row, col] = ujac_pairs[i];
        const double expected = row == col ? unary_derivative(uxv[row]) : 0.0;
        if (!near(ujac_values[i], expected)) {
            std::cerr << "unexpected vector unary sparse Jacobian entry (" << row << ", " << col << "): " << ujac_values[i] << " expected " << expected << "\n";
            return 1;
        }
    }
    auto ujac_c = to_sparse_jacobian_c(UF, "x", ujac_pairs, "unary_vector_jac");
    if (ujac_c.find("void unary_vector_jac") == std::string::npos) {
        std::cerr << "generated vector unary sparse Jacobian C did not contain expected function\n";
        return 1;
    }

    GraphFunctionBuilder coll_builder;
    auto cx = coll_builder.inputs("x", 3);
    DenseMatrix CD(3, 3, {1.0, -2.0, 0.5, 0.0, 3.0, 4.0, -1.0, 0.0, 2.0});
    auto lin = coll_builder.dense_matvec(CD, cx);
    auto rhs = coll_builder.vector({coll_builder.at(cx, 0) * coll_builder.at(cx, 0), sin(coll_builder.at(cx, 1)), coll_builder.at(cx, 2)});
    auto CF = coll_builder.function(cx, coll_builder.sub(lin, coll_builder.scale(0.25, rhs)));
    auto CVJP = reverse_diff(CF, "lambda", "x");
    if (!CVJP.has_vector_structure()) {
        std::cerr << "collocation-like VJP did not preserve vector structure\n";
        return 1;
    }
    bool coll_vjp_saw_dense = false;
    for (const auto &node : CVJP.vector_nodes) {
        coll_vjp_saw_dense = coll_vjp_saw_dense || node.op == VectorOp::DenseMatVec;
    }
    if (!coll_vjp_saw_dense) {
        std::cerr << "collocation-like VJP did not keep dense transpose structured\n";
        return 1;
    }
    auto CHVP = forward_diff(CVJP, "x", "v");
    if (!CHVP.has_vector_structure()) {
        std::cerr << "collocation-like HVP did not preserve vector structure\n";
        return 1;
    }
    auto cjac_pairs = jacobian_sparsity(CF, "x").to_pairs();
    auto cjac_fn = sparse_jacobian_function(CF, "x", cjac_pairs);
    double cxv[3] = {1.0, 0.3, -2.0};
    std::vector<double> cjac_values(cjac_pairs.size(), 0.0);
    VM(cjac_fn).evaluate(EvalEnv().input("x", cxv), cjac_values.data());
    for (std::size_t i = 0; i < cjac_pairs.size(); ++i) {
        const auto [row, col] = cjac_pairs[i];
        double expected = CD(row, col);
        if (row == 0 && col == 0) {
            expected -= 0.5 * cxv[0];
        } else if (row == 1 && col == 1) {
            expected -= 0.25 * std::cos(cxv[1]);
        } else if (row == 2 && col == 2) {
            expected -= 0.25;
        }
        if (!near(cjac_values[i], expected)) {
            std::cerr << "unexpected mixed residual sparse Jacobian entry (" << row << ", " << col << "): " << cjac_values[i] << " expected " << expected << "\n";
            return 1;
        }
    }
    auto cjac_c = to_sparse_jacobian_c(CF, "x", cjac_pairs, "collocation_mixed_jac");
    if (cjac_c.find("void collocation_mixed_jac") == std::string::npos || cjac_c.find("constant_index") == std::string::npos ||
        cjac_c.find("constant_value") == std::string::npos || cjac_c.find("cos(") == std::string::npos) {
        std::cerr << "generated mixed residual sparse Jacobian C did not keep matrix constants and nonlinear terms together\n";
        return 1;
    }
    GraphFunctionBuilder vjp_builder;
    auto vx = vjp_builder.inputs("x", 2);
    DenseMatrix A(2, 2, {1.0, 2.0, 3.0, 4.0});
    auto AF = vjp_builder.function(vx, vjp_builder.dense_matvec(A, vx));
    auto VJP = reverse_diff(AF, "lambda", "x");
    if (!VJP.has_vector_structure()) {
        std::cerr << "reverse-mode transform did not preserve vector structure for dense matvec\n";
        return 1;
    }
    bool vjp_saw_dense_matvec = false;
    for (const auto &node : VJP.vector_nodes) {
        vjp_saw_dense_matvec = vjp_saw_dense_matvec || node.op == VectorOp::DenseMatVec;
    }
    if (!vjp_saw_dense_matvec) {
        std::cerr << "reverse-mode transform did not use transpose dense matvec metadata\n";
        return 1;
    }
    double vxv[2] = {0.0, 0.0};
    double lambda[2] = {5.0, 7.0};
    double gout[2] = {0.0, 0.0};
    VM(VJP).evaluate(EvalEnv().input("x", vxv).param("lambda", lambda), gout);
    if (!near(gout[0], 26.0) || !near(gout[1], 38.0)) {
        std::cerr << "unexpected dense matvec VJP: " << gout[0] << ", " << gout[1] << "\n";
        return 1;
    }
    auto vjp_c = to_c(VJP, "dense_matvec_vjp");
    if (vjp_c.find("void dense_matvec_vjp") == std::string::npos || vjp_c.find("static const double mat") == std::string::npos ||
        vjp_c.find("moo_ad_dense_matvec") == std::string::npos) {
        std::cerr << "generated VJP C did not use vector-aware transpose dense matvec lowering\n";
        return 1;
    }

    auto HVP = forward_diff(VJP, "x", "v");
    if (!HVP.has_vector_structure()) {
        std::cerr << "HVP transform did not preserve vector structure for linear dense matvec\n";
        return 1;
    }
    auto H = hessian_sparsity(HVP, "v");
    if (!H.empty()) {
        std::cerr << "linear dense matvec HVP produced nonzero Hessian sparsity\n";
        return 1;
    }
    double hv[2] = {11.0, 13.0};
    double hout[2] = {1.0, 1.0};
    VM(HVP).evaluate(EvalEnv().input("x", vxv).param("lambda", lambda).param("v", hv), hout);
    if (!near(hout[0], 0.0) || !near(hout[1], 0.0)) {
        std::cerr << "linear dense matvec HVP was not zero: " << hout[0] << ", " << hout[1] << "\n";
        return 1;
    }
    StagedVM staged_hvp(HVP, "v");
    auto prepared = staged_hvp.prepare(EvalEnv().input("x", vxv).param("lambda", lambda));
    hout[0] = 1.0;
    hout[1] = 1.0;
    prepared.apply(hv, hout);
    if (!near(hout[0], 0.0) || !near(hout[1], 0.0)) {
        std::cerr << "staged linear dense matvec HVP was not zero: " << hout[0] << ", " << hout[1] << "\n";
        return 1;
    }
    auto hvp_c = to_c(HVP, "dense_matvec_hvp");
    if (hvp_c.find("void dense_matvec_hvp") == std::string::npos) {
        std::cerr << "generated HVP C did not contain expected function name\n";
        return 1;
    }
    auto staged_c = to_staged_c(HVP, "dense_matvec_hvp_staged", "v");
    if (staged_c.find("dense_matvec_hvp_staged_prepare") == std::string::npos || staged_c.find("dense_matvec_hvp_staged_apply") == std::string::npos) {
        std::cerr << "generated staged HVP C did not contain expected functions\n";
        return 1;
    }

    GraphFunctionBuilder sparse_builder;
    auto sx = sparse_builder.inputs("x", 4);
    SparseMatrix S(3, 4, {0, 0, 1, 2}, {0, 2, 1, 3}, {2.0, -1.0, 3.0, 4.0});
    auto SF = sparse_builder.function(sx, sparse_builder.sparse_matvec(S, sx));
    if (!SF.has_vector_structure()) {
        std::cerr << "sparse matvec function did not preserve vector structure\n";
        return 1;
    }
    bool saw_sparse_matvec = false;
    for (const auto &node : SF.vector_nodes) {
        saw_sparse_matvec = saw_sparse_matvec || node.op == VectorOp::SparseMatVec;
    }
    if (!saw_sparse_matvec) {
        std::cerr << "sparse matvec vector node was not preserved\n";
        return 1;
    }
    double sxv[4] = {1.0, 2.0, 3.0, 4.0};
    double sout[3] = {0.0, 0.0, 0.0};
    VM(SF).evaluate(EvalEnv().input("x", sxv), sout);
    if (!near(sout[0], -1.0) || !near(sout[1], 6.0) || !near(sout[2], 16.0)) {
        std::cerr << "unexpected sparse matvec value: " << sout[0] << ", " << sout[1] << ", " << sout[2] << "\n";
        return 1;
    }
    auto SJVP = forward_diff(SF, "x", "v");
    bool sjvp_saw_sparse = false;
    for (const auto &node : SJVP.vector_nodes) {
        sjvp_saw_sparse = sjvp_saw_sparse || node.op == VectorOp::SparseMatVec;
    }
    if (!sjvp_saw_sparse) {
        std::cerr << "sparse matvec JVP did not preserve sparse vector node\n";
        return 1;
    }
    double sv[4] = {0.5, 1.0, -2.0, 3.0};
    double sjout[3] = {0.0, 0.0, 0.0};
    VM(SJVP).evaluate(EvalEnv().input("x", sxv).param("v", sv), sjout);
    if (!near(sjout[0], 3.0) || !near(sjout[1], 3.0) || !near(sjout[2], 12.0)) {
        std::cerr << "unexpected sparse matvec JVP: " << sjout[0] << ", " << sjout[1] << ", " << sjout[2] << "\n";
        return 1;
    }
    auto sparse_c = to_c(SF, "sparse_matvec_value");
    if (sparse_c.find("void sparse_matvec_value") == std::string::npos || sparse_c.find("sp_row") == std::string::npos || sparse_c.find("sp_col") == std::string::npos ||
        sparse_c.find("moo_ad_sparse_matvec") == std::string::npos) {
        std::cerr << "generated sparse matvec C did not use vector-aware sparse lowering\n";
        return 1;
    }

    auto SVJP = reverse_diff(SF, "lambda", "x");
    bool svjp_saw_sparse = false;
    for (const auto &node : SVJP.vector_nodes) {
        svjp_saw_sparse = svjp_saw_sparse || node.op == VectorOp::SparseMatVec;
    }
    if (!svjp_saw_sparse) {
        std::cerr << "sparse matvec VJP did not preserve transpose sparse vector node\n";
        return 1;
    }
    double slambda[3] = {5.0, 7.0, 11.0};
    double sgout[4] = {0.0, 0.0, 0.0, 0.0};
    VM(SVJP).evaluate(EvalEnv().input("x", sxv).param("lambda", slambda), sgout);
    if (!near(sgout[0], 10.0) || !near(sgout[1], 21.0) || !near(sgout[2], -5.0) || !near(sgout[3], 44.0)) {
        std::cerr << "unexpected sparse matvec VJP: " << sgout[0] << ", " << sgout[1] << ", " << sgout[2] << ", " << sgout[3] << "\n";
        return 1;
    }
    auto SHVP = forward_diff(SVJP, "x", "v");
    auto SH = hessian_sparsity(SHVP, "v");
    if (!SH.empty()) {
        std::cerr << "linear sparse matvec HVP produced nonzero Hessian sparsity\n";
        return 1;
    }
    GraphFunctionBuilder kron_builder;
    auto kx = kron_builder.inputs("x", 6);
    DenseMatrix K(2, 3, {1.0, 2.0, 0.0, -1.0, 0.0, 4.0});
    auto KF = kron_builder.function(kx, kron_builder.kron_eye_matvec(K, 2, kx));
    if (!KF.has_vector_structure()) {
        std::cerr << "kron-eye matvec function did not preserve vector structure\n";
        return 1;
    }
    bool saw_kron = false;
    for (const auto &node : KF.vector_nodes) {
        saw_kron = saw_kron || node.op == VectorOp::KronEyeMatVec;
    }
    if (!saw_kron) {
        std::cerr << "kron-eye matvec vector node was not preserved\n";
        return 1;
    }
    double kxv[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    double kout[4] = {0.0, 0.0, 0.0, 0.0};
    VM(KF).evaluate(EvalEnv().input("x", kxv), kout);
    if (!near(kout[0], 7.0) || !near(kout[1], 10.0) || !near(kout[2], 19.0) || !near(kout[3], 22.0)) {
        std::cerr << "unexpected kron-eye matvec value: " << kout[0] << ", " << kout[1] << ", " << kout[2] << ", " << kout[3] << "\n";
        return 1;
    }
    auto KJVP = forward_diff(KF, "x", "v");
    bool kjvp_saw_kron = false;
    for (const auto &node : KJVP.vector_nodes) {
        kjvp_saw_kron = kjvp_saw_kron || node.op == VectorOp::KronEyeMatVec;
    }
    if (!kjvp_saw_kron) {
        std::cerr << "kron-eye matvec JVP did not preserve vector node\n";
        return 1;
    }
    double kv[6] = {0.5, 1.0, -1.0, 2.0, 3.0, -2.0};
    double kjout[4] = {0.0, 0.0, 0.0, 0.0};
    VM(KJVP).evaluate(EvalEnv().input("x", kxv).param("v", kv), kjout);
    if (!near(kjout[0], -1.5) || !near(kjout[1], 5.0) || !near(kjout[2], 11.5) || !near(kjout[3], -9.0)) {
        std::cerr << "unexpected kron-eye matvec JVP: " << kjout[0] << ", " << kjout[1] << ", " << kjout[2] << ", " << kjout[3] << "\n";
        return 1;
    }
    auto kron_c = to_c(KF, "kron_eye_value");
    if (kron_c.find("void kron_eye_value") == std::string::npos || kron_c.find("kron_mat") == std::string::npos || kron_c.find("moo_ad_kron_eye_matvec") == std::string::npos) {
        std::cerr << "generated kron-eye C did not use vector-aware lowering\n";
        return 1;
    }
    auto kron_jac_c = to_sparse_jacobian_c(KF, "x", jacobian_sparsity(KF, "x").to_pairs(), "kron_eye_jac");
    if (kron_jac_c.find("void kron_eye_jac") == std::string::npos || kron_jac_c.find("constant_index") == std::string::npos ||
        kron_jac_c.find("constant_value") == std::string::npos) {
        std::cerr << "generated kron-eye sparse Jacobian C did not use vector-aware coefficient lowering\n";
        return 1;
    }
    auto KVJP = reverse_diff(KF, "lambda", "x");
    bool kvjp_saw_kron = false;
    for (const auto &node : KVJP.vector_nodes) {
        kvjp_saw_kron = kvjp_saw_kron || node.op == VectorOp::KronEyeMatVec;
    }
    if (!kvjp_saw_kron) {
        std::cerr << "kron-eye matvec VJP did not preserve transpose vector node\n";
        return 1;
    }
    double klambda[4] = {7.0, 11.0, 13.0, 17.0};
    double kgout[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    VM(KVJP).evaluate(EvalEnv().input("x", kxv).param("lambda", klambda), kgout);
    if (!near(kgout[0], -6.0) || !near(kgout[1], -6.0) || !near(kgout[2], 14.0) || !near(kgout[3], 22.0) || !near(kgout[4], 52.0) || !near(kgout[5], 68.0)) {
        std::cerr << "unexpected kron-eye matvec VJP: " << kgout[0] << ", " << kgout[1] << ", " << kgout[2] << ", " << kgout[3] << ", " << kgout[4] << ", " << kgout[5] << "\n";
        return 1;
    }
    return 0;
}
