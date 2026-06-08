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
