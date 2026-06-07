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

from moo import matrix, nlp_model, radauIIA, vec


STAGES = 3
STATE_DIM = 2


def f(x, u):
    return vec([
        -(u[0] + 0.5 * u[0] * u[0]) * x[0],
        u[0] * x[0],
    ])

def M(x, u):
    return -x[1]


def build_model(name: str, intervals: int = 50, segments: int = 1):
    """
    intervals : number of collocation intervals
    segments : number of intervals sharing one control value (segments=1 -> one u per interval)
    intervals must be divisible by segments.
    """
    assert intervals % segments == 0, f"intervals={intervals} must be divisible by segments={segments}"

    r = radauIIA(STAGES)
    model = nlp_model(name)

    n_controls = intervals // segments

    x1 = model.add_variables("x1", intervals * STAGES, guess=1.0, lb=0.0, ub=1.0)
    x2 = model.add_variables("x2", intervals * STAGES, guess=0.0, lb=0.0, ub=1.0)
    u  = model.add_variables("u", n_controls, guess=2.0, lb=0.0, ub=5.0)

    h = 1.0 / intervals
    D_otimes_I = matrix(r.D).otimes_eye(STATE_DIM)
    x_initial = vec([1.0, 0.0])

    def node(interval, node):
        return interval * STAGES + node

    def state(interval, stage):
        return vec([x1[node(interval, stage)], x2[node(interval, stage)]])

    def controls(interval):
        return vec([u[interval // segments]])

    def state_block(previous, interval):
        return vec([previous] + [state(interval, stage) for stage in range(STAGES)])

    def collocation_state(block, stage):
        return block.block(stage + 1, STATE_DIM)

    def collocation_derivative(block, stage):
        return (D_otimes_I @ block).block(stage + 1, STATE_DIM)

    def defect(interval, stage, block):
        x_stage = collocation_state(block, stage)
        return collocation_derivative(block, stage) - h * f(x_stage, controls(interval))

    model.add_constraints(
        range(0, 1),
        lambda interval: [
            value
            for stage in range(STAGES)
            for value in defect(interval, stage, state_block(x_initial, interval))
        ],
        eq=0.0,
        name="initial_collocation",
    )

    model.add_constraints(
        range(1, intervals),
        lambda interval: [
            value
            for stage in range(STAGES)
            for value in defect(interval, stage, state_block(state(interval - 1, STAGES - 1), interval))
        ],
        eq=0.0,
        name="collocation",
    )

    model.minimize(M(state(intervals - 1, STAGES - 1), controls(intervals - 1)), name="final_x2")

    model.codegen("loop-direct")
    model.solver(tolerance=1e-12)

    return model


def plot_solution(run, out: Path, intervals: int, segments: int):
    import matplotlib
    import matplotlib.pyplot as plt

    r = radauIIA(STAGES)
    values = run.result.variables
    h = 1.0 / intervals

    t = [0.0]
    x1 = [1.0]
    x2 = [0.0]
    for interval in range(intervals):
        for node, c in enumerate(r.C):
            idx = interval * STAGES + node
            t.append((interval + c) * h)
            x1.append(values[f"x1{idx}"])
            x2.append(values[f"x2{idx}"])

    tu = [i * h for i in range(intervals + 1)]
    uu = [values[f"u{interval // segments}"] for interval in range(intervals)]
    uu.append(uu[-1])

    fig, (ax_x, ax_u) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    ax_x.plot(t, x1, ".-", label="x1")
    ax_x.plot(t, x2, ".-", label="x2")
    ax_x.set_ylabel("state")
    ax_x.grid(True, alpha=0.3)
    ax_x.legend()

    ax_u.step(tu, uu, where="post", label=f"u  (segments={segments})")
    ax_u.set_xlabel("time")
    ax_u.set_ylabel("control")
    ax_u.grid(True, alpha=0.3)
    ax_u.legend()

    fig.tight_layout()
    out.mkdir(parents=True, exist_ok=True)
    path = out / "solution.png"
    fig.savefig(path, dpi=160)
    if matplotlib.get_backend().lower() not in {"agg", "pdf", "svg", "ps", "template"}:
        plt.show()
    return path


if __name__ == "__main__":
    out = Path("build/moo/batch_reactor_piecewise_ctrls")
    intervals = 100
    segments = 2

    model = build_model("BatchReactorPiecewise", intervals=intervals, segments=segments)
    c_path, h_path = model.generate(out)
    exe_path = model.compile(out)
    run = model.optimize(out)

    print(plot_solution(run, out, intervals, segments))
    raise SystemExit(run.returncode)
