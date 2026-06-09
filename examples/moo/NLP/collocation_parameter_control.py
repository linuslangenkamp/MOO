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
from math import tanh as math_tanh

from moo import blockvec, exp, matrix, nlp_model, radauIIA, vec


STATE_DIM = 2
HIDDEN_DIM = 4
FIRST_WEIGHT_OFFSET = 0
FIRST_BIAS_OFFSET = FIRST_WEIGHT_OFFSET + HIDDEN_DIM * STATE_DIM
SECOND_WEIGHT_OFFSET = FIRST_BIAS_OFFSET + HIDDEN_DIM
SECOND_BIAS_OFFSET = SECOND_WEIGHT_OFFSET + HIDDEN_DIM * HIDDEN_DIM
OUTPUT_WEIGHT_OFFSET = SECOND_BIAS_OFFSET + HIDDEN_DIM
OUTPUT_BIAS_OFFSET = OUTPUT_WEIGHT_OFFSET + HIDDEN_DIM
PARAMETER_DIM = OUTPUT_BIAS_OFFSET + 1


def activation(z):
    return 2.0 / (1.0 + exp(-2.0 * z)) - 1.0


def control_nn(x, p):
    hidden_1 = [
        activation(
            p[FIRST_WEIGHT_OFFSET + j * STATE_DIM] * x[0]
            + p[FIRST_WEIGHT_OFFSET + j * STATE_DIM + 1] * x[1]
            + p[FIRST_BIAS_OFFSET + j]
        )
        for j in range(HIDDEN_DIM)
    ]
    hidden_2 = [
        activation(
            sum(p[SECOND_WEIGHT_OFFSET + j * HIDDEN_DIM + k] * hidden_1[k] for k in range(HIDDEN_DIM))
            + p[SECOND_BIAS_OFFSET + j]
        )
        for j in range(HIDDEN_DIM)
    ]
    return p[OUTPUT_BIAS_OFFSET] + sum(p[OUTPUT_WEIGHT_OFFSET + j] * hidden_2[j] for j in range(HIDDEN_DIM))


def f(x, p):
    u = control_nn(x, p)
    return vec([
        -(u + 0.5 * u * u) * x[0],
        u * x[0],
    ])


def M(x, p):
    return -x[1]

def build_model(name: str = "BatchReactor", intervals: int = 25, stages: int = 3):
    r = radauIIA(stages)
    model = nlp_model(name)

    state_values = intervals * stages + 1

    x = [
        model.add_variables("x0", state_values, guess=1.0),
        model.add_variables("x1", state_values, guess=0.0),
    ]

    p = model.add_variables("p", PARAMETER_DIM, guess=0.0)
    directions = [(1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)]
    for j, (w0, w1) in enumerate(directions[:HIDDEN_DIM]):
        p.set_guess(FIRST_WEIGHT_OFFSET + j * STATE_DIM, w0)
        p.set_guess(FIRST_WEIGHT_OFFSET + j * STATE_DIM + 1, w1)
    for j in range(HIDDEN_DIM):
        p.set_guess(SECOND_WEIGHT_OFFSET + j * HIDDEN_DIM + j, 1.0)
        p.set_guess(OUTPUT_WEIGHT_OFFSET + j, 0.05 if j % 2 == 0 else -0.05)
    p.set_guess(OUTPUT_BIAS_OFFSET, 2.0)

    h = 1 / intervals
    D_otimes_I = matrix(r.D).otimes_eye(STATE_DIM)

    initial_states = [1.0, 0.0]

    for i in range(STATE_DIM):
        x[i].fix(0, initial_states[i])

    def state_index(interval, stage):
        return interval * stages + stage + 1

    def state(interval, stage):
        idx = state_index(interval, stage)
        return vec([x[i][idx] for i in range(STATE_DIM)])

    def left_state(interval):
        if interval == 0:
            return vec(initial_states)
        return state(interval - 1, stages - 1)

    def collocation_defects(interval):
        stage_states = [state(interval, stage) for stage in range(stages)]
        X = blockvec([left_state(interval)] + stage_states)
        derivative = D_otimes_I @ X
        return blockvec([
            derivative.block(stage + 1, STATE_DIM) - h * f(stage_states[stage], p)
            for stage in range(stages)
        ])

    model.add_constraints(range(intervals), lambda interval: collocation_defects(interval), eq=0.0, name="collocation")
    model.add_constraints(
        range(intervals),
        lambda interval: vec([control_nn(state(interval, stage), p) for stage in range(stages)]),
        lb=0.0,
        ub=5.0,
        name="control_bounds",
    )
    model.minimize(M(state(intervals - 1, stages - 1), p), name="mayer")

    model.codegen("direct")
    model.solver(tolerance=1e-8)
    return model


def plot_solution(run, out: Path, intervals: int, stages: int):
    import matplotlib
    import matplotlib.pyplot as plt

    r = radauIIA(stages)
    values = run.result.variables
    h = 1.0 / intervals
    params = [values[f"p{i}"] for i in range(PARAMETER_DIM)]

    def u_of(x0, x1):
        hidden_1 = [
            math_tanh(
                params[FIRST_WEIGHT_OFFSET + j * STATE_DIM] * x0
                + params[FIRST_WEIGHT_OFFSET + j * STATE_DIM + 1] * x1
                + params[FIRST_BIAS_OFFSET + j]
            )
            for j in range(HIDDEN_DIM)
        ]
        hidden_2 = [
            math_tanh(
                sum(params[SECOND_WEIGHT_OFFSET + j * HIDDEN_DIM + k] * hidden_1[k] for k in range(HIDDEN_DIM))
                + params[SECOND_BIAS_OFFSET + j]
            )
            for j in range(HIDDEN_DIM)
        ]
        return params[OUTPUT_BIAS_OFFSET] + sum(params[OUTPUT_WEIGHT_OFFSET + j] * hidden_2[j] for j in range(HIDDEN_DIM))

    t = [0.0]
    x0 = [1.0]
    x1 = [0.0]
    tu = []
    uu = []

    for interval in range(intervals):
        for stage, c in enumerate(r.C):
            idx = interval * stages + stage + 1
            time = (interval + c) * h
            x0_value = values[f"x0{idx}"]
            x1_value = values[f"x1{idx}"]
            t.append(time)
            x0.append(x0_value)
            x1.append(x1_value)
            tu.append(time)
            uu.append(u_of(x0_value, x1_value))

    fig, (ax_x, ax_u) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    ax_x.plot(t, x0, ".-", label="x0")
    ax_x.plot(t, x1, ".-", label="x1")
    ax_x.set_ylabel("state")
    ax_x.grid(True, alpha=0.3)
    ax_x.legend()

    ax_u.plot(tu, uu, ".-", label="u(x, p)")
    ax_u.set_xlabel("time")
    ax_u.set_ylabel("control")
    ax_u.grid(True, alpha=0.3)
    ax_u.legend(title=f"{HIDDEN_DIM}-hidden-unit feedback NN")

    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    path = out / "solution.png"
    fig.savefig(path, dpi=160)
    if matplotlib.get_backend().lower() not in {"agg", "pdf", "svg", "ps", "template"}:
        plt.show()
    return path


def print_overview(intervals: int, stages: int):
    print("Batch reactor NN feedback overview")
    print(f"  collocation: {intervals} intervals x {stages} Radau IIA stages")
    print(f"  states: {STATE_DIM}")
    print(f"  controller: {STATE_DIM} -> {HIDDEN_DIM} tanh -> {HIDDEN_DIM} tanh -> 1 direct output")
    print(f"  controller parameters: {PARAMETER_DIM}")
    print("  control range: 0 <= u(x, p) <= 5 via path constraints")
    print(f"  objective: maximize final x1")


if __name__ == "__main__":
    out = Path("build/moo/NLP/batch_reactor_parameters")
    intervals = 25
    stages = 5
    run = build_model(intervals=intervals, stages=stages).run(out)
    print_overview(intervals, stages)
    print(plot_solution(run, out, intervals, stages))
    raise SystemExit(run.returncode)
