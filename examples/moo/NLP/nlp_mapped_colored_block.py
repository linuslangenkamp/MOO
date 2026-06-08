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

from moo import nlp_model


def build_model(name: str = "nlp_mapped_colored_block", blocks: int = 200, block_size: int = 64):
    model = nlp_model(name)
    x = model.add_variables("x", blocks + block_size, lb=-10.0, ub=10.0, guess=0.1, nominal=1.0)

    model.minimize_sum(blocks, lambda i: (x[i] - 1.0) ** 2, name="tracking")
    model.add_constraints_map(
        blocks,
        lambda i: [
            (x[i + j] + x[i + j + 1]) ** 2
            for j in range(block_size - 1)
        ],
        lb=-1.0,
        ub=25.0,
        name="band_block",
    )
    model.solver(tolerance=1e-10, derivative_test=False)
    return model


if __name__ == "__main__":
    out = Path("build/moo/nlp_mapped_colored_block")
    model = build_model().codegen("colored")
    c_path, _ = model.generate(out)
    exe_path = model.compile(out)
    print(c_path)
    print(exe_path)
    print((out / "codegen_report.txt").read_text(), end="")
