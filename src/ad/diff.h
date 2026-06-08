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

#ifndef MOO_AD_DIFF_H
#define MOO_AD_DIFF_H

#include "optimize.h"

namespace ad {

std::vector<Expr> clone_nodes(const Graph &src, OptimizingBuilder &dst, const std::optional<std::string> &input_as_param = std::nullopt);
GraphFunction forward_diff(const GraphFunction &f, const std::string &wrt_input_name = "x", const std::string &direction_name = "v");
GraphFunction reverse_diff(const GraphFunction &f, const std::string &lambda_name = "lambda", const std::string &wrt_input_name = "x");

} // namespace ad

#endif // MOO_AD_DIFF_H
