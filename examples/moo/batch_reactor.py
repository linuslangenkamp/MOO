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
    model = gdop_model("BatchReactor")
    x1 = model.add_state("x1", start=1.0)
    x2 = model.add_state("x2", start=0.0)
    u = model.add_control("u", lb=0.0, ub=5.0)

    model.set_time_fixed(t0 = 0.0, tf = 1.0)
    model.mesh(intervals=40, nodes=3)

    model.set_dynamics(x1, -(u + u**2 / 2) * x1)
    model.set_dynamics(x2, u * x1)

    model.add_mayer(-x2)

    model.mesh_refinement(0, 0)

    model.solver(tolerance=1e-12)

    return model


if __name__ == "__main__":
    out = Path("build/moo/batchReactor")
    result = build_model().run(out)
    result.result.plot.all(show=True)
    raise SystemExit(result.returncode)
