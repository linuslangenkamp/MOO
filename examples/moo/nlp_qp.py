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

from moo import nlp_model


def build_model():
    model = nlp_model("nlp_qp_codegen")
    x = model.add_variable("x", lb=-10.0, ub=10.0, guess=1.0)
    y = model.add_variable("y", lb=-10.0, ub=10.0, guess=1.0)

    model.minimize((x - 1.0) ** 2 + (y - 2.0) ** 2)
    model.add_constraint(x + y, eq=2.0, name="sum")
    model.solver(backend="Uno", tolerance=1e-12, derivative_test=True)
    return model


if __name__ == "__main__":
    out = Path("build/moo/nlp_qp")
    model = build_model()
    c_path, h_path = model.generate(out)
    exe_path = model.compile(out)
    result = model.optimize(out, solver="Ipopt")
    print(result.result.variables, flush=True)
    raise SystemExit(result.returncode)
