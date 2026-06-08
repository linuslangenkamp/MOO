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

from pathlib import Path

from moo import gdop_model, exp


def build_model():
    model = gdop_model("OilShalePyrolysis")

    x1 = model.add_state("kerogen", start=1.0)
    x2 = model.add_state("pyrolytic bitumen", start=0.0)
    x3 = model.add_state("oil & gas", start=0.0)
    x4 = model.add_state("organic carbon", start=0.0)

    T = model.add_control(
        "temperature",
        lb=698.15,
        ub=748.15,
        nominal=700,
        guess=(748.15 + 698.15) / 2,
    )

    model.set_time_fixed(t0=0.0, tf=8.0)
    model.mesh(intervals=16, nodes=3)

    k1 = exp(8.86 - (20300 / 1.9872) / T)
    k2 = exp(24.25 - (37400 / 1.9872) / T)
    k3 = exp(23.67 - (33800 / 1.9872) / T)
    k4 = exp(18.75 - (28200 / 1.9872) / T)
    k5 = exp(20.70 - (31000 / 1.9872) / T)

    model.set_dynamics(x1, -k1 * x1 - (k3 + k4 + k5) * x1 * x2)
    model.set_dynamics(x2, k1 * x1 - k2 * x2 + k3 * x1 * x2)
    model.set_dynamics(x3, k2 * x2 + k4 * x1 * x2)
    model.set_dynamics(x4, k5 * x1 * x2)

    model.add_mayer(-x2)

    model.solver(tolerance=1e-8)

    model.mesh_refinement(2, 2)

    return model


if __name__ == "__main__":
    out = Path("build/moo/oilShalePyrolysis")
    result = build_model().run(out)
    result.result.plot.all(save=out / "solution.png")
    raise SystemExit(result.returncode)
