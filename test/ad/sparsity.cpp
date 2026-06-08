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

#include <iostream>
#include <vector>

int main() {
    using namespace ad;

    GraphFunctionBuilder g;
    auto x = g.inputs("x", 3);
    auto out = g.vector({
        2.0 * g.at(x, 0) - g.at(x, 1) + 3.0,
        g.at(x, 2),
    });

    auto F = g.function(x, out);
    auto Grad = reverse_diff(F, "lambda", "x");
    auto HVP = forward_diff(Grad, "x", "v");

    auto J = structural_sparsity(F, "x");
    auto H = hessian_sparsity(HVP, "v");

    if (J.nnz() != 3) {
        std::cerr << "expected 3 Jacobian nonzeros, got " << J.nnz() << "\n";
        return 1;
    }
    if (!H.empty()) {
        std::cerr << "expected 0 Hessian nonzeros, got " << H.size() << "\n";
        return 1;
    }
    return 0;
}
