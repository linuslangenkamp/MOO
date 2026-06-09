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

from moo import blockvec, matrix, nlp_model, radauIIA, vec


def f(x, u):
    return vec([
        -(u[0] + 0.5 * u[0] * u[0]) * x[0],
        u[0] * x[0],
    ])


def M(x, u):
    return -x[1]


def build_model(name: str = "BatchReactor", intervals: int = 25, stages: int = 3):
    r = radauIIA(stages)
    model = nlp_model(name)

    num_states = 2
    num_controls = 1

    state_values = intervals * stages + 1
    control_values = intervals * stages

    x = [
        model.add_variables("x0", state_values, guess=1.0),
        model.add_variables("x1", state_values, guess=0.0),
    ]

    u = [ model.add_variables("u", control_values, lb=0.0, ub=5.0, guess=2) ]

    h = 1 / intervals
    D_otimes_I = matrix(r.D).otimes_eye(num_states)

    initial_states = [1.0, 0.0]

    for i in range(num_states):
        x[i].fix(0, initial_states[i])

    def state_index(interval, stage):
        return interval * stages + stage + 1

    def control_index(interval, stage):
        return interval * stages + stage

    def state(interval, stage):
        idx = state_index(interval, stage)
        return vec([x[i][idx] for i in range(num_states)])

    def left_state(interval):
        if interval == 0:
            return vec(initial_states)
        return state(interval - 1, stages - 1)

    def control(interval, stage):
        idx = control_index(interval, stage)
        return vec([u[i][idx] for i in range(num_controls)])

    def collocation_defects(interval):
        X = blockvec([left_state(interval)] + [state(interval, stage) for stage in range(stages)])
        derivative = D_otimes_I @ X
        return blockvec([derivative.block(stage + 1, num_states) - h * f(state(interval, stage), control(interval, stage)) for stage in range(stages)])

    model.add_constraints(range(intervals), lambda interval: collocation_defects(interval), eq=0.0, name="collocation")

    model.minimize(M(state(intervals - 1, stages - 1), control(intervals - 1, stages - 1)), name="mayer")

    model.codegen("direct")
    model.solver(tolerance=1e-10)
    return model

if __name__ == "__main__":
    out = Path("build/moo/NLP/batch_reactor")
    intervals = 100
    stages = 3
    run = build_model(intervals=intervals, stages=stages).run(out)
    raise SystemExit(run.returncode)
