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

from __future__ import annotations

from dataclasses import dataclass

from . import ad


@dataclass
class EmittedFunction:
    value: str
    jvp: str
    hvp: str
    jac_sparsity: list[tuple[int, int]]
    hes_sparsity: list[tuple[int, int]]


def _safe_env(builder: ad.GraphBuilder, variables: dict[str, object]) -> dict[str, object]:
    env: dict[str, object] = {
        "__builtins__": {},
        "sin": ad.sin,
        "cos": ad.cos,
        "tan": ad.tan,
        "exp": ad.exp,
        "log": ad.log,
        "pow_const": ad.pow_const,
    }
    env.update(variables)
    return env


def _as_expr(builder: ad.GraphBuilder, value: object):
    if isinstance(value, (int, float)):
        return builder.constant(float(value))
    return value


def emit_function(
    input_name: str,
    input_size: int,
    outputs: list[str | None],
    params: list[tuple[str, int]],
    value_name: str,
    jvp_name: str,
    hvp_name: str,
) -> EmittedFunction:
    builder = ad.GraphBuilder()
    x = builder.inputs(input_name, input_size)
    param_vectors = {name: builder.params(name, size) for name, size in params}
    env = _safe_env(builder, {input_name: x, **param_vectors})

    out = []
    for code in outputs:
        if code is None:
            continue
        out.append(_as_expr(builder, eval(code, env)))
    if not out:
        out.append(builder.constant(0.0))

    fn = builder.function(x, out, list(param_vectors.values()))
    jvp = fn.forward_diff(input_name, "v")
    grad = fn.reverse_diff("lambda", input_name)
    hvp = grad.forward_diff(input_name, "v")
    return EmittedFunction(
        value=fn.to_c(value_name),
        jvp=jvp.to_c(jvp_name),
        hvp=hvp.to_staged_c(hvp_name, "v"),
        jac_sparsity=fn.jacobian_sparsity(input_name),
        hes_sparsity=hvp.hessian_sparsity("v"),
    )


def emit_value_function(
    input_name: str,
    input_size: int,
    outputs: list[str | None],
    params: list[tuple[str, int]],
    value_name: str,
) -> tuple[str, list[tuple[int, int]]]:
    builder = ad.GraphBuilder()
    x = builder.inputs(input_name, input_size)
    param_vectors = {name: builder.params(name, size) for name, size in params}
    env = _safe_env(builder, {input_name: x, **param_vectors})

    out = [_as_expr(builder, eval(code, env)) for code in outputs if code is not None]
    if not out:
        out.append(builder.constant(0.0))

    fn = builder.function(x, out, list(param_vectors.values()))
    return fn.to_c(value_name), fn.jacobian_sparsity(input_name)
