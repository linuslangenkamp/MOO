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

from pathlib import Path

from moo import gdop_model


def build_model():
    model = gdop_model("representative_codegen")

    gain = model.add_runtime_parameter("gain", 2.0)
    x = model.add_state("x", start=1.0, final=0.0)
    y = model.add_state("y", start=0.0)
    u = model.add_control("u")
    v = model.add_control("v")
    p = model.add_parameter("p", lb=0.0, ub=4.0, guess=1.0)

    model.set_time_fixed(tf=1.0)
    model.mesh(intervals=20, nodes=3)

    model.set_dynamics(x, u)
    model.set_dynamics(y, p + v)

    model.add_path(v, eq=0.0)
    model.add_boundary(y - gain, eq=0.0)

    model.add_lagrange((u + 1.0) * (u + 1.0) + v * v + (p - gain) * (p - gain))
    model.add_mayer((x * x) + ((y - gain) * (y - gain)) + ((p - gain) * (p - gain)))
    model.solver(tolerance=1e-12, derivative_test=True)
    return model


if __name__ == "__main__":
    out = Path("build/moo/representative")
    result = build_model().run(out, solver="Ipopt")
    raise SystemExit(result.returncode)
