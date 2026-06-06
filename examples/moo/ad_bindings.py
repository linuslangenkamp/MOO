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
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from moo import ad


g = ad.GraphBuilder()
x = g.inputs("x", 3)
f = g.function(x, [2.0 * x[0] - x[1] + 3.0, x[2] * x[2] + 1.0])

grad = f.reverse_diff("lambda", "x")
hvp = grad.forward_diff("x", "v")

print(f.jacobian_sparsity("x"))
print(hvp.hessian_sparsity("v"))
print(f.evaluate(inputs={"x": [1.0, 2.0, 3.0]}))
print(grad.evaluate(inputs={"x": [1.0, 2.0, 3.0]}, params={"lambda": [1.0, 1.0]}))
print(hvp.evaluate(inputs={"x": [1.0, 2.0, 3.0]}, params={"lambda": [1.0, 1.0], "v": [0.0, 0.0, 1.0]}))
prepared = hvp.prepare(inputs={"x": [1.0, 2.0, 3.0]}, params={"lambda": [1.0, 1.0]})
print(prepared.apply([0.0, 0.0, 1.0]))
print(f.to_c("demo_value"))
