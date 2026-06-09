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
        -x[0]**3 + u[0],
    ])


def L(x, u):
    return 0.5 * (x[0]**2 + u[0]**2)


def build_model(name: str = "SpectralHypersensitiveFreeTime", intervals: int = 1, stages: int = 100):
    r = radauIIA(stages)
    model = nlp_model(name)

    num_states = 1
    num_controls = 1

    state_values = intervals * stages + 1
    control_values = intervals * stages

    x = [
        model.add_variables("x", state_values, guess=1.0),
    ]

    u = [
        model.add_variables("u", control_values, guess=0.0),
    ]

    tf = model.add_variables("tf", 1, lb=0.5, guess=50.0, nominal=50.0)

    D_otimes_I = matrix(r.D).otimes_eye(num_states)

    initial_states = [1.0]

    for i in range(num_states):
        x[i].fix(0, initial_states[i])

    def final_time():
        return tf[0]

    def step_size():
        return final_time() / intervals

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
        h = step_size()
        X = blockvec([left_state(interval)] + [state(interval, stage) for stage in range(stages)])
        derivative = D_otimes_I @ X
        return blockvec([
            derivative.block(stage + 1, num_states) - vec([h * f(state(interval, stage), control(interval, stage))[0]])
            for stage in range(stages)
        ])

    def lagrange(interval):
        h = step_size()
        return sum(h * r.b[stage] * L(state(interval, stage), control(interval, stage)) for stage in range(stages))

    model.add_constraints(range(intervals), lambda interval: collocation_defects(interval), eq=0.0, name="collocation")
    model.add_constraint(state(intervals - 1, stages - 1)[0], eq=1.5, name="final_state")

    model.minimize_sum(range(intervals), lambda interval: lagrange(interval), name="tracking")

    model.codegen("direct")
    model.solver(tolerance=1e-12)
    return model


if __name__ == "__main__":
    out = Path("build/moo/NLP/spectral_hypersensitive_free_time")
    intervals = 1
    stages = 100
    run = build_model(intervals=intervals, stages=stages).run(out)
    print(f"tf = {run.result.variables['tf0']:.12g}")
    raise SystemExit(run.returncode)
