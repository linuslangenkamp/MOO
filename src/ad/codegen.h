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

#ifndef MOO_AD_CODEGEN_H
#define MOO_AD_CODEGEN_H

#include "sparse_derivatives.h"
#include "sparsity.h"
#include "vm.h"

namespace ad {

struct CEmitter {
    static std::string number(double v);
    static std::string cname(const std::string &s);
    static std::string cache_slot(const Node &n);
    static std::string node_ref(const Node &n);
    static std::string expr_rhs(const Graph &g, NodeId id, const std::function<std::string(NodeId)> &ref);
    static bool emit_vector_function(const GraphFunction &f, const std::string &name, std::ostream &os);
    static void emit_function(const GraphFunction &fn, const std::string &name, std::ostream &os);
    static void emit_staged(const GraphFunction &fn, const std::string &basename, const std::string &direction_name, std::ostream &os);
};

std::string to_c(const GraphFunction &f, const std::string &name);
std::string to_staged_c(const GraphFunction &f, const std::string &basename, const std::string &direction_name = "v");
std::string to_sparse_coefficients_c(const GraphFunction &linear_f, const std::string &seed_name, const std::vector<std::pair<int, int>> &entries, const std::string &name);
std::string to_sparse_jacobian_c(const GraphFunction &F,
                                 const std::string &wrt,
                                 const std::vector<std::pair<int, int>> &entries,
                                 const std::string &name,
                                 const std::string &direction_name = "v");
std::string to_sparse_hessian_c(const GraphFunction &HVP, const std::string &direction_name, const std::vector<std::pair<int, int>> &entries, const std::string &name);
std::string derivative_callback_mode(const std::string &strategy, const std::vector<std::pair<int, int>> &pairs, const std::vector<int> &colors);
std::vector<std::string> derivative_section_keys(const std::string &prefix, const std::string &jac_mode, const std::string &hes_mode);
std::string render_jacobian_callback_body(const std::string &mode,
                                          const std::string &direct_call,
                                          const std::string &input_size,
                                          const std::string &output_size,
                                          const std::vector<std::pair<int, int>> &pairs,
                                          const std::vector<int> &colors,
                                          const std::string &colored_call);
std::string render_hessian_callback_body(const std::string &mode,
                                         const std::string &direct_body,
                                         const std::string &input_size,
                                         const std::string &tmp_size,
                                         const std::vector<std::pair<int, int>> &pairs,
                                         const std::vector<int> &colors,
                                         const std::string &prepare_body,
                                         const std::string &apply_call,
                                         const std::vector<int> &buf_indices = {});

struct ExactDerivativeCode {
    std::string value;
    std::string jvp;
    std::string hvp;
    std::string jacobian;
    std::string hessian;
    std::vector<std::pair<int, int>> jacobian_sparsity;
    std::vector<std::pair<int, int>> hessian_sparsity;
    std::vector<int> jacobian_colors;
    std::vector<int> hessian_colors;
    int jacobian_color_count = 0;
    int hessian_color_count = 0;
    std::size_t value_bytes = 0;
    std::size_t jvp_bytes = 0;
    std::size_t hvp_bytes = 0;
    std::size_t jacobian_bytes = 0;
    std::size_t hessian_bytes = 0;
};

ExactDerivativeCode emit_exact_derivative_code(const GraphFunction &F,
                                               const std::string &wrt,
                                               const std::string &direction_name,
                                               const std::string &lambda_name,
                                               const std::string &value_name,
                                               const std::string &jvp_name,
                                               const std::string &hvp_name);

} // namespace ad

#endif // MOO_AD_CODEGEN_H
