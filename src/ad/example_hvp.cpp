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

#include <advec.h>
#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    using namespace advec;

    Graph g;
    auto x = g.inputs("x", 2);

    std::vector<Expr> y(2);
    y[0] = sin(x[0]) * x[1];
    y[1] = pow_const(x[0], 2) + exp(x[1]);

    auto F = function_from(std::move(g), x, y);

    auto Grad = reverse_diff(F, "lambda", "x");
    auto JVP = forward_diff(F, "x", "v");
    auto HVP = forward_diff(Grad, "x", "v");

    VM vm(HVP);

    double xv[2] = {1.0, 2.0};
    double lambda[2] = {3.0, 4.0};
    double v[2] = {0.1, -0.2};
    double out[2] = {0.0, 0.0};

    EvalEnv env;
    env.input("x", xv).param("lambda", lambda).param("v", v);
    vm.evaluate(env, out);

    std::cout << "HVP VM: [" << out[0] << ", " << out[1] << "]\n";

    // Analytic check for f = [sin(x0)*x1, x0^2 + exp(x1)]
    double ref0 = (-lambda[0] * std::sin(xv[0]) * xv[1] + 2.0 * lambda[1]) * v[0] + (lambda[0] * std::cos(xv[0])) * v[1];
    double ref1 = (lambda[0] * std::cos(xv[0])) * v[0] + (lambda[1] * std::exp(xv[1])) * v[1];
    std::cout << "Reference: [" << ref0 << ", " << ref1 << "]\n";
    std::cout << "HVP - Ref: [" << out[0] - ref0 << ", " << out[1] - ref1 << "]\n";

    // Staged VM: cache all work independent of v, then apply many directions.
    StagedVM staged(HVP, "v");
    EvalEnv prep_env;
    prep_env.input("x", xv).param("lambda", lambda);
    auto prepared = staged.prepare(prep_env);

    prepared.apply(v, out);
    std::cout << "Staged HVP: [" << out[0] << ", " << out[1] << "]\n";
    std::cout << "Staged HVP - Ref: [" << out[0] - ref0 << ", " << out[1] - ref1 << "]\n";

    VM vm_grad(Grad);
    double gp[2], g0[2];
    const double h = 1e-8;
    double xp[2] = {xv[0] + h * v[0], xv[1] + h * v[1]};
    vm_grad.evaluate(EvalEnv{}.input("x", xp).param("lambda", lambda), gp);
    vm_grad.evaluate(EvalEnv{}.input("x", xv).param("lambda", lambda), g0);
    std::cout << "HVP FD: [" << (gp[0] - g0[0]) / (h) << ", " << (gp[1] - g0[1]) / (h) << "]\n";
    std::cout << "HVP FD - Ref: [" << (gp[0] - g0[0]) / (h)-ref0 << ", " << (gp[1] - g0[1]) / (h)-ref1 << "]\n";

    std::cout << "\nOne-shot generated C\n";
    std::cout << to_c(HVP, "hvp_eval") << "\n";

    std::cout << "\nStaged generated C\n";
    std::cout << to_staged_c(HVP, "hvp", "v") << "\n";

    auto J_pat = structural_sparsity(Grad, "lambda");
    std::cout << "\nTransposed Jacobian sparsity of F\n" << J_pat;

    J_pat = structural_sparsity(JVP, "x");
    std::cout << "\nJacobian sparsity of F\n" << J_pat;

    J_pat = structural_sparsity(F, "x");
    std::cout << "\nJacobian sparsity of F\n" << J_pat;

    auto H_lower = hessian_sparsity(HVP, "v");
    std::cout << "\nLagrangian Hessian sparsity (lower triangle)\n";
    for (auto &[row, col] : H_lower)
    {
        std::cout << "  H[" << row << "," << col << "]\n";
    }

    auto H_full = hessian_sparsity_full(HVP, "v");
    std::cout << "\nLagrangian Hessian sparsity (full symmetric)\n";
    for (auto &[row, col] : H_full)
    {
        std::cout << "  H[" << row << "," << col << "]\n";
    }
    return 0;
}
