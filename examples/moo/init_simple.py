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

from moo import init_model


def build_model():
    model = init_model("init_simple_codegen")
    y = model.add_variable("y", guess=0.0)
    p = model.add_parameter("p", lb=0.0, ub=5.0, base=1.0)
    dp = model.delta(p)

    model.add_f(y + p - 3.0)
    model.add_g(y, lb=0.0)
    model.set_objective((y - 1.0) * (y - 1.0) + dp * dp)
    model.solver(tolerance=1e-12, derivative_test=True, qp=True)
    return model


if __name__ == "__main__":
    out = Path("build/moo/init_simple")
    model = build_model()
    c_path, h_path = model.generate(out)
    exe_path = model.compile(out)
    result = model.optimize(out, solver="Ipopt")
    raise SystemExit(result.returncode)
