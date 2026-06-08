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


STAGES = 100
STATE_DIM = 1


def f(x, u):
    return vec([
        -x[0]**3 + u[0]
    ])


def L(x, u):
    return 0.5 * (x[0]**2 + u[0]**2)


def build_model(name: str, intervals: int = 1, tf: float = 50.0):
    r = radauIIA(STAGES)
    model = nlp_model(name)
    x = model.add_variables("x", intervals * STAGES, guess=1.0)
    u = model.add_variables("u", intervals * STAGES, guess=0.0)

    h = tf / intervals
    D_otimes_I = matrix(r.D).otimes_eye(STATE_DIM)
    x_initial = vec([1.0])

    def node(interval, j):
        return interval * STAGES + j

    def state(interval, stage):
        return vec([x[node(interval, stage)]])

    def state_nodes(previous, interval):
        return blockvec([previous] + [state(interval, stage) for stage in range(STAGES)])

    def collocation_defects(previous, interval):
        nodes = state_nodes(previous, interval)
        derivative = D_otimes_I @ nodes
        return blockvec([
            derivative.block(stage + 1, STATE_DIM) - h * f(nodes[stage + 1], vec([u[node(interval, stage)]]))
            for stage in range(STAGES)
        ])

    def lagrange(interval):
        return sum(h * r.b[stage] * L(state(interval, stage), vec([u[node(interval, stage)]])) for stage in range(STAGES))

    model.add_constraints(
        range(0, 1),
        lambda interval: collocation_defects(x_initial, interval),
        eq=0.0,
        name="initial_collocation",
    )
    model.add_constraints(
        range(1, intervals),
        lambda interval: collocation_defects(state(interval - 1, STAGES - 1), interval),
        eq=0.0,
        name="collocation",
    )

    model.add_constraint(x[-1], eq=1.5)

    model.minimize_sum(range(intervals), lambda interval: lagrange(interval), name="tracking")
    model.codegen("direct")
    model.solver(tolerance=1e-12)
    return model


def plot_solution(run, out: Path, intervals: int, tf: float):
    import matplotlib
    import matplotlib.pyplot as plt

    r = radauIIA(STAGES)
    values = run.result.variables
    h = tf / intervals

    t = [0.0]
    x = [1.0]
    tu = []
    uu = []
    for interval in range(intervals):
        for j, c in enumerate(r.C):
            idx = interval * STAGES + j
            time = (interval + c) * h
            t.append(time)
            x.append(values[f"x{idx}"])
            tu.append(time)
            uu.append(values[f"u{idx}"])

    fig, (ax_x, ax_u) = plt.subplots(2, 1, sharex=True, figsize=(9, 6))
    ax_x.plot(t, x, ".-", label="x")
    ax_x.set_ylabel("state")
    ax_x.grid(True, alpha=0.3)
    ax_x.legend()

    ax_u.plot(tu, uu, ".-", label="u at collocation nodes")
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
    out = Path("build/moo/spectral_hypersensitive")
    intervals = 1
    tf = 50.0
    run = build_model("SpectralHypersensitive", intervals=intervals, tf=tf).run(out)
    print(plot_solution(run, out, intervals, tf))
    raise SystemExit(run.returncode)
