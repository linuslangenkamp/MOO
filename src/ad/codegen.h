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

} // namespace ad

#endif // MOO_AD_CODEGEN_H
