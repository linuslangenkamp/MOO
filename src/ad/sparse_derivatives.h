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

#ifndef MOO_AD_SPARSE_DERIVATIVES_H
#define MOO_AD_SPARSE_DERIVATIVES_H

#include "diff.h"

namespace ad {

struct SparseDerivativePlan {
    std::vector<std::pair<int, int>> pattern;
    std::vector<int> column_color;
    int color_count = 0;
};

std::vector<int> greedy_column_coloring(int column_count, const std::vector<std::pair<int, int>> &pattern);
SparseDerivativePlan derivative_plan(int column_count, const std::vector<std::pair<int, int>> &pattern);

struct ExactDerivativePlan {
    GraphFunction jvp;
    GraphFunction grad;
    GraphFunction hvp;
    std::vector<std::pair<int, int>> jacobian_sparsity;
    std::vector<std::pair<int, int>> hessian_sparsity;
    std::vector<std::pair<int, int>> hessian_full_sparsity;
    std::vector<int> jacobian_colors;
    std::vector<int> hessian_colors;
    int jacobian_color_count = 0;
    int hessian_color_count = 0;
};

int graph_function_group_size(const GraphFunction &f, const std::string &name);

struct LinearForm {
    Expr constant;
    std::map<int, Expr> coeff;
};

bool has_coeff(const LinearForm &f);
LinearForm make_constant_form(Expr e);
LinearForm linear_add(OptimizingBuilder &b, const LinearForm &a, const LinearForm &c, double sign = 1.0);
LinearForm linear_scale(OptimizingBuilder &b, const LinearForm &a, Expr scale);
GraphFunction sparse_coefficients(const GraphFunction &fn, const std::string &seed_name, const std::vector<std::pair<int, int>> &entries);
GraphFunction sparse_jacobian_function(const GraphFunction &F,
                                       const std::string &wrt,
                                       const std::vector<std::pair<int, int>> &entries,
                                       const std::string &direction_name = "v");
GraphFunction sparse_hessian_function(const GraphFunction &HVP, const std::string &direction_name, const std::vector<std::pair<int, int>> &entries);

} // namespace ad

#endif // MOO_AD_SPARSE_DERIVATIVES_H
