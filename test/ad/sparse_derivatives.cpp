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
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <ad.h>

#include <cmath>
#include <iostream>
#include <vector>

namespace {

bool close(double a, double b) {
    return std::abs(a - b) <= 1e-12;
}

} // namespace

int main() {
    using namespace ad;

    GraphFunctionBuilder g;
    auto x = g.inputs("x", 3);
    auto out = g.vector({
        g.at(x, 0) * g.at(x, 0) + g.at(x, 1),
        g.at(x, 1) * g.at(x, 2),
    });

    auto F = g.function(x, out);
    auto J_pat = jacobian_sparsity(F, "x").to_pairs();
    auto J = sparse_jacobian_function(F, "x", J_pat);

    double xv[3] = {2.0, 3.0, 5.0};
    double jac[4] = {};
    VM jvm(J);
    EvalEnv jenv;
    jenv.input("x", xv);
    jvm.evaluate(jenv, jac);

    std::vector<double> expected_jac = {4.0, 1.0, 5.0, 3.0};
    if (J_pat != std::vector<std::pair<int, int>>{{0, 0}, {0, 1}, {1, 1}, {1, 2}}) {
        std::cerr << "unexpected Jacobian sparsity\n";
        return 1;
    }
    for (int i = 0; i < 4; ++i) {
        if (!close(jac[i], expected_jac[i])) {
            std::cerr << "bad sparse Jacobian value " << i << ": " << jac[i] << "\n";
            return 1;
        }
    }

    auto Grad = reverse_diff(F, "lambda", "x");
    auto HVP = forward_diff(Grad, "x", "v");
    auto H_pat = hessian_sparsity(HVP, "v");
    auto H = sparse_hessian_function(HVP, "v", H_pat);

    double lambda[2] = {7.0, 11.0};
    double hes[2] = {};
    VM hvm(H);
    EvalEnv henv;
    henv.input("x", xv).param("lambda", lambda);
    hvm.evaluate(henv, hes);

    if (H_pat != std::vector<std::pair<int, int>>{{0, 0}, {2, 1}}) {
        std::cerr << "unexpected Hessian sparsity\n";
        return 1;
    }
    if (!close(hes[0], 14.0) || !close(hes[1], 11.0)) {
        std::cerr << "bad sparse Hessian values: " << hes[0] << ", " << hes[1] << "\n";
        return 1;
    }

    return 0;
}
