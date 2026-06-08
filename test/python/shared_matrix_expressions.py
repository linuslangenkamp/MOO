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
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

from pathlib import Path
import shutil

import moo.ad as ad
from moo import Expr, gdop_model, init_model, matrix, nlp_model, sparse_matrix, vec
from moo.graph_expression import GraphExpressionEmitter


def assert_backend_dense_matvec_jvp() -> None:
    builder = ad.GraphFunctionBuilder()
    x = builder.inputs("x", 2)
    u = builder.params("u", 1)
    D = [[1.0, 2.0], [3.0, 4.0]]

    dx = builder.dense_matvec(D, x)
    rhs = builder.vector([x[1], -x[0] + u[0] * u[0]])
    residual = builder.vector_sub(dx, builder.vector_scale(0.5, rhs))
    fn = builder.function(x, residual, [u])

    value = fn.evaluate(inputs={"x": [2.0, 3.0]}, params={"u": [4.0]})
    assert value == [6.5, 11.0]

    jvp = fn.forward_diff("x", "v")
    direction = {"u": [4.0], "v": [0.25, -0.5]}
    got = jvp.evaluate(inputs={"x": [2.0, 3.0]}, params=direction)
    expected = [
        1.0 * 0.25 + 2.0 * -0.5 - 0.5 * -0.5,
        3.0 * 0.25 + 4.0 * -0.5 - 0.5 * (-0.25),
    ]
    assert all(abs(a - b) < 1e-12 for a, b in zip(got, expected)), (got, expected)


def assert_backend_vector_composition() -> None:
    builder = ad.GraphFunctionBuilder()
    x = builder.inputs("x", 4)
    first = builder.slice(x, 0, 2)
    second = builder.slice(x, 2, 2)
    y = builder.vector_add(builder.vector_scale(2.0, first), second)
    fn = builder.function(x, y)

    assert fn.evaluate(inputs={"x": [1.0, 2.0, 3.0, 4.0]}) == [5.0, 8.0]
    jvp = fn.forward_diff("x", "v")
    got = jvp.evaluate(inputs={"x": [1.0, 2.0, 3.0, 4.0]}, params={"v": [0.5, 1.0, 1.5, 2.0]})
    assert got == [2.5, 4.0]


def assert_graph_expression_emitter_unwraps_backend_handles() -> None:
    builder = ad.GraphFunctionBuilder()
    x = builder.inputs("x", 2)
    emitter = GraphExpressionEmitter(builder, {"x": x})

    scalar = Expr(graph_scalar=x[0] + 2.0)
    emitted_scalar = emitter.output(scalar)
    assert emitted_scalar.source == "native_backend_scalar"
    assert emitted_scalar.value is not None
    assert emitted_scalar.vector is None

    vector_expr = vec(x)
    shifted = vector_expr + vector_expr
    emitted_vector = emitter.output(shifted)
    assert emitted_vector.source == "native_backend_vector"
    assert emitted_vector.value is None
    assert emitted_vector.vector is not None

    matvec = matrix([[1.0, 2.0], [3.0, 4.0]]) @ vec(x)
    emitted_matvec = emitter.output(matvec)
    assert emitted_matvec.source == "native_backend_vector"
    assert emitted_matvec.vector is not None

    fn = builder.function(x, emitted_vector.vector)
    assert fn.evaluate(inputs={"x": [3.0, 4.0]}) == [6.0, 8.0]


def assert_backend_sparse_matvec_jvp() -> None:
    builder = ad.GraphFunctionBuilder()
    x = builder.inputs("x", 4)
    y = builder.sparse_matvec([0, 0, 1, 2], [0, 2, 1, 3], [2.0, -1.0, 3.0, 4.0], (3, 4), x)
    fn = builder.function(x, y)

    assert fn.evaluate(inputs={"x": [1.0, 2.0, 3.0, 4.0]}) == [-1.0, 6.0, 16.0]
    jvp = fn.forward_diff("x", "v")
    got = jvp.evaluate(inputs={"x": [1.0, 2.0, 3.0, 4.0]}, params={"v": [0.5, 1.0, -2.0, 3.0]})
    assert got == [3.0, 3.0, 12.0]


def assert_backend_kron_eye_matvec_jvp() -> None:
    builder = ad.GraphFunctionBuilder()
    x = builder.inputs("x", 6)
    y = builder.kron_eye_matvec([[1.0, 2.0, 0.0], [-1.0, 0.0, 4.0]], 2, x)
    fn = builder.function(x, y)

    assert fn.evaluate(inputs={"x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]}) == [7.0, 10.0, 19.0, 22.0]
    jvp = fn.forward_diff("x", "v")
    got = jvp.evaluate(inputs={"x": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]}, params={"v": [0.5, 1.0, -1.0, 2.0, 3.0, -2.0]})
    assert got == [-1.5, 5.0, 11.5, -9.0]


def assert_dense_matvec_structure() -> None:
    M = matrix([[1.0, 2.0], [3.0, 4.0]])
    x = vec([Expr.variable("x", 0), Expr.variable("x", 1)])
    y = M @ x
    builder = ad.GraphFunctionBuilder()
    xb = builder.inputs("x", 2)
    emitted = GraphExpressionEmitter(builder, {"x": xb}).output(y)
    assert emitted.source == "native_backend_vector"
    fn = builder.function(xb, emitted.vector)
    assert fn.has_vector_structure


def assert_sparse_matvec_structure() -> None:
    S = sparse_matrix([0, 0, 1], [0, 2, 1], [2.0, -1.0, 3.0], (2, 3))
    x = vec([Expr.variable("x", 0), Expr.variable("x", 1), Expr.variable("x", 2)])
    y = S @ x
    builder = ad.GraphFunctionBuilder()
    xb = builder.inputs("x", 3)
    emitted = GraphExpressionEmitter(builder, {"x": xb}).output(y)
    assert emitted.source == "native_backend_vector"
    fn = builder.function(xb, emitted.vector)
    assert fn.has_vector_structure


def assert_kron_eye_matvec_structure() -> None:
    K = matrix([[1.0, 2.0, 0.0], [-1.0, 0.0, 4.0]]).otimes_eye(2)
    x = vec([Expr.variable("x", i) for i in range(6)])
    y = K @ x
    assert len(y) == 4
    builder = ad.GraphFunctionBuilder()
    xb = builder.inputs("x", 6)
    emitted = GraphExpressionEmitter(builder, {"x": xb}).output(y)
    assert emitted.source == "native_backend_vector"
    fn = builder.function(xb, emitted.vector)
    assert fn.has_vector_structure


def assert_python_expressions_are_backend_buildable() -> None:
    x0 = Expr.variable("x", 0)
    x1 = Expr.variable("x", 1)
    product = (x0 + 2.0) * x1
    summed = x0 + x1 + 3.0
    y = matrix([[1.0, 2.0]]) @ vec([x0, x1])
    builder = ad.GraphFunctionBuilder()
    xb = builder.inputs("x", 2)
    emitter = GraphExpressionEmitter(builder, {"x": xb})
    assert emitter.output(product).source == "native_backend_scalar"
    assert emitter.output(summed).source == "native_backend_scalar"
    assert emitter.output(y).source == "native_backend_vector"


def assert_dense_matvec_matches_scalar_fallback(out: Path) -> None:
    out = out / "dense_matvec_equivalence"
    M = matrix([[1.0, 2.0], [3.0, 4.0]])

    structured = nlp_model("DenseMatvecStructured")
    xs = structured.add_variables("x", 2, guess=1.0)
    ys = M @ xs
    structured.minimize(0.0)
    structured.add_constraints(ys, eq=0.0, name="matvec")
    structured.generate(out / "structured")

    scalar = nlp_model("DenseMatvecScalar")
    xq = scalar.add_variables("x", 2, guess=1.0)
    scalar.minimize(0.0)
    scalar.add_constraint(1.0 * xq[0] + 2.0 * xq[1], eq=0.0, name="matvec_0")
    scalar.add_constraint(3.0 * xq[0] + 4.0 * xq[1], eq=0.0, name="matvec_1")
    scalar.generate(out / "scalar")

    structured_report = (out / "structured" / "codegen_report.txt").read_text(encoding="utf-8")
    scalar_report = (out / "scalar" / "codegen_report.txt").read_text(encoding="utf-8")
    assert "graph_source=native_expr" in structured_report, structured_report
    assert "ad_graph_skipped" not in structured_report, structured_report
    assert "structured_dense_matvec" not in structured_report, structured_report
    for expected in ("jacobian_nnz=4", "hessian_nnz=0"):
        assert expected in structured_report, structured_report
        assert expected in scalar_report, scalar_report


def assert_sparse_matvec_matches_scalar_fallback(out: Path) -> None:
    out = out / "sparse_matvec_equivalence"
    S = sparse_matrix([0, 0, 1], [0, 2, 1], [2.0, -1.0, 3.0], (2, 3))

    structured = nlp_model("SparseMatvecStructured")
    xs = structured.add_variables("x", 3, guess=1.0)
    ys = S @ xs
    structured.minimize(0.0)
    structured.add_constraints(ys, eq=0.0, name="spmatvec")
    structured.generate(out / "structured")

    scalar = nlp_model("SparseMatvecScalar")
    xq = scalar.add_variables("x", 3, guess=1.0)
    scalar.minimize(0.0)
    scalar.add_constraint(2.0 * xq[0] - xq[2], eq=0.0, name="spmatvec_0")
    scalar.add_constraint(3.0 * xq[1], eq=0.0, name="spmatvec_1")
    scalar.generate(out / "scalar")

    structured_report = (out / "structured" / "codegen_report.txt").read_text(encoding="utf-8")
    scalar_report = (out / "scalar" / "codegen_report.txt").read_text(encoding="utf-8")
    assert "graph_source=native_expr" in structured_report, structured_report
    for expected in ("jacobian_nnz=3", "hessian_nnz=0"):
        assert expected in structured_report, structured_report
        assert expected in scalar_report, scalar_report


def assert_dense_matvec_solve_matches_scalar(out: Path) -> None:
    out = out / "dense_matvec_equivalence"

    def build(name: str, structured: bool):
        model = nlp_model(name)
        x = model.add_variables("x", 2, lb=-10.0, ub=10.0, guess=1.0)
        model.minimize(0.0)
        if structured:
            y = matrix([[1.0, 2.0], [3.0, 4.0]]) @ x
            model.add_constraint(y[0], eq=0.0, name="g0")
            model.add_constraint(y[1], eq=0.0, name="g1")
        else:
            model.add_constraint(1.0 * x[0] + 2.0 * x[1], eq=0.0, name="g0")
            model.add_constraint(3.0 * x[0] + 4.0 * x[1], eq=0.0, name="g1")
        return model

    structured = build("DenseMatvecSolveStructured", True).run(out / "dense_matvec_solve_structured", capture=True)
    scalar = build("DenseMatvecSolveScalar", False).run(out / "dense_matvec_solve_scalar", capture=True)
    assert structured.returncode == 0, structured.process.stdout
    assert scalar.returncode == 0, scalar.process.stdout
    structured_report = (out / "dense_matvec_solve_structured" / "codegen_report.txt").read_text(encoding="utf-8")
    assert "graph_source=native_expr" in structured_report, structured_report
    for name in ("x0", "x1"):
        assert abs(structured.result.variables[name] - scalar.result.variables[name]) < 1e-10


def main() -> None:
    out = Path("/tmp/moo_shared_matrix_expression_test")
    if out.exists():
        shutil.rmtree(out)
    assert_backend_dense_matvec_jvp()
    assert_backend_vector_composition()
    assert_graph_expression_emitter_unwraps_backend_handles()
    assert_backend_sparse_matvec_jvp()
    assert_backend_kron_eye_matvec_jvp()
    assert_dense_matvec_structure()
    assert_sparse_matvec_structure()
    assert_kron_eye_matvec_structure()
    assert_python_expressions_are_backend_buildable()
    assert_dense_matvec_matches_scalar_fallback(out)
    assert_sparse_matvec_matches_scalar_fallback(out)
    assert_dense_matvec_solve_matches_scalar(out)
    M = matrix([[1.0, 2.0], [3.0, 4.0]])

    nlp = nlp_model("SharedMatrixNLP")
    x = nlp.add_variables("x", 2, guess=1.0)
    y = M @ x
    nlp.minimize(y[0] * y[0] + y[1] * y[1])
    nlp.generate(out / "nlp")

    init = init_model("SharedMatrixInit")
    a = init.add_variable("a", guess=1.0)
    b = init.add_variable("b", guess=1.0)
    r = M @ vec([a, b])
    init.set_objective(r[0] * r[0] + r[1] * r[1])
    init.generate(out / "init")

    gdop = gdop_model("SharedMatrixGDOP")
    s0 = gdop.add_state("s0", start=1.0, final=0.0)
    s1 = gdop.add_state("s1", start=0.0, final=0.0)
    u = gdop.add_control("u")
    r = M @ vec([s0, s1])
    gdop.set_time_fixed(0.0, 1.0)
    gdop.mesh(intervals=2, nodes=2)
    gdop.set_dynamics(s0, r[0] + u)
    gdop.set_dynamics(s1, r[1])
    gdop.add_lagrange(u * u)
    gdop.generate(out / "gdop")

    for report in out.glob("*/codegen_report.txt"):
        text = report.read_text(encoding="utf-8")
        assert "graph_source=native_expr" in text or "structured_loop_blocks" in text, report


if __name__ == "__main__":
    main()
