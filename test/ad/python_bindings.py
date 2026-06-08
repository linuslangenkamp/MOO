# SPDX-License-Identifier: LGPL-3.0-or-later
#
# This file is part of MOO - Modelica / Model Optimizer
# Copyright (C) 2026 University of Applied Sciences and Arts
# Bielefeld, Faculty of Engineering and Mathematics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from moo import ad


g = ad.GraphFunctionBuilder()
x = g.inputs("x", 3)
f = g.function(x, g.vector([2.0 * x[0] - x[1] + 3.0, x[2]]))
hvp = f.reverse_diff("lambda", "x").forward_diff("x", "v")

assert f.has_vector_structure
assert f.vector_node_count >= 2
assert f.jacobian_sparsity("x") == [(0, 0), (0, 1), (1, 2)]
assert hvp.hessian_sparsity("v") == []
assert f.evaluate(inputs={"x": [1.0, 2.0, 3.0]}) == [3.0, 3.0]
assert hvp.evaluate(
    inputs={"x": [1.0, 2.0, 3.0]},
    params={"lambda": [1.0, 1.0], "v": [0.0, 0.0, 1.0]},
) == [0.0, 0.0, 0.0]
prepared = hvp.prepare(inputs={"x": [1.0, 2.0, 3.0]}, params={"lambda": [1.0, 1.0]})
assert prepared.apply([0.0, 0.0, 1.0]) == [0.0, 0.0, 0.0]
assert "void linear_value" in f.to_c("linear_value")
assert "linear_hvp_prepare" in hvp.to_staged_c("linear_hvp", "v")
assert "void linear_jac" in f.to_sparse_jacobian_c("x", f.jacobian_sparsity("x"), "linear_jac")
assert "void linear_hes" in hvp.to_sparse_hessian_c("v", hvp.hessian_sparsity("v"), "linear_hes")
assert f.coloring(f.jacobian_sparsity("x"), 3)["color_count"] == 2

plan = f.exact_derivative_plan("x", "v", "lambda")
assert plan["jacobian_sparsity"] == [(0, 0), (0, 1), (1, 2)]
assert plan["hessian_sparsity"] == []
assert plan["jacobian_color_count"] == 2
assert plan["hessian_color_count"] == 0
assert plan["jvp"].evaluate(
    inputs={"x": [1.0, 2.0, 3.0]},
    params={"v": [1.0, 1.0, 1.0]},
) == [1.0, 1.0]
assert plan["hvp"].evaluate(
    inputs={"x": [1.0, 2.0, 3.0]},
    params={"lambda": [1.0, 1.0], "v": [1.0, 1.0, 1.0]},
) == [0.0, 0.0, 0.0]
