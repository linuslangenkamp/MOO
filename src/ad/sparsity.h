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

#ifndef MOO_AD_SPARSITY_H
#define MOO_AD_SPARSITY_H

#include "sparse_derivatives.h"

namespace ad {

struct SparsityEntry {
    int output;
    int var;
};

struct SparsityPattern {
    std::string variable_name;
    int output_size = 0;
    int variable_size = 0;
    std::vector<SparsityEntry> entries;

    std::size_t nnz() const;
    bool empty() const;
    std::vector<std::pair<int, int>> to_pairs() const;
    std::vector<std::pair<int, int>> contract_outputs() const;
    static std::vector<std::pair<int, int>> lower_triangular(std::vector<std::pair<int, int>> pairs);
    static std::vector<std::pair<int, int>> symmetrize(std::vector<std::pair<int, int>> pairs);
    std::string str() const;
};

std::ostream &operator<<(std::ostream &os, const SparsityPattern &p);
SparsityPattern structural_sparsity(const GraphFunction &fn, const std::string &wrt_name);
SparsityPattern jacobian_sparsity(const GraphFunction &F, const std::string &wrt = "x");
std::vector<std::pair<int, int>> hessian_sparsity(const GraphFunction &HVP, const std::string &v_name = "v");
std::vector<std::pair<int, int>> hessian_sparsity_full(const GraphFunction &HVP, const std::string &v_name = "v");
ExactDerivativePlan exact_derivative_plan(const GraphFunction &F,
                                          const std::string &wrt = "x",
                                          const std::string &direction_name = "v",
                                          const std::string &lambda_name = "lambda");

} // namespace ad

#endif // MOO_AD_SPARSITY_H
