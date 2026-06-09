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

from __future__ import annotations

from dataclasses import dataclass

from . import ad
from .graph_expression import GraphExpressionEmitter


@dataclass
class EmittedFunction:
    value: str
    jvp: str
    hvp: str
    jac: str
    hes: str
    jac_sparsity: list[tuple[int, int]]
    hes_sparsity: list[tuple[int, int]]
    jac_colors: list[int]
    hes_colors: list[int]
    report: dict[str, object]


@dataclass
class BuiltGraphFunction:
    function: object
    report: dict[str, object]


def dedupe_ad_vector_helpers(sections: list[str]) -> list[str]:
    marker = "#ifndef MOO_AD_VECTOR_HELPERS_DEFINED"
    end_marker = "#endif\n\n"
    seen = False
    out = []
    for section in sections:
        start = section.find(marker)
        if start < 0:
            out.append(section)
            continue
        end = section.find(end_marker, start)
        if end < 0:
            out.append(section)
            continue
        end += len(end_marker)
        if seen:
            out.append(section[:start] + section[end:])
        else:
            seen = True
            out.append(section)
    return out


def build_graph_function(
    input_name: str,
    input_size: int,
    outputs: list[object],
    params: list[tuple[str, int]],
) -> BuiltGraphFunction:
    builder = ad.GraphFunctionBuilder()
    x = builder.inputs(input_name, input_size)
    param_vectors = {name: builder.params(name, size) for name, size in params}
    return build_graph_function_from_builder(builder, input_name, x, outputs, param_vectors)


def build_graph_function_from_builder(
    builder: object,
    input_name: str,
    input_vector: object,
    outputs: list[object],
    param_vectors: dict[str, object],
) -> BuiltGraphFunction:
    compiler = GraphExpressionEmitter(builder, {input_name: input_vector, **param_vectors})

    out = []
    vector_out = None
    sources: set[str] = set()
    for output in outputs:
        compiled = compiler.output(output)
        sources.add(compiled.source)
        if compiled.vector is not None and compiled.value is None and vector_out is None and not out:
            vector_out = compiled.vector
        elif compiled.vector is not None:
            for value in compiled.vector.values:
                out.append(value)
            vector_out = None
        elif compiled.value is not None:
            out.append(compiled.value)
    if not out:
        if vector_out is None:
            out.append(builder.constant(0.0))
            if not sources:
                sources.add("native_constant")

    native_sources = {"native_expr", "native_vector_expr", "native_backend_scalar", "native_backend_vector", "native_constant", "empty"}
    source_summary = "native_expr" if sources and sources <= native_sources else "unknown"
    fn = builder.function(input_vector, vector_out if vector_out is not None else builder.vector(out), list(param_vectors.values()))
    return BuiltGraphFunction(
        function=fn,
        report={
            "graph_source": source_summary,
            "graph_sources": ",".join(sorted(sources)),
            "outputs": len(vector_out) if vector_out is not None else len(out),
            "has_vector_structure": bool(fn.has_vector_structure),
            "vector_node_count": int(fn.vector_node_count),
            "value_lowering": "vector" if bool(fn.has_vector_structure) else "scalar",
        },
    )


def _emit_built_function(built: BuiltGraphFunction, input_name: str, value_name: str, jvp_name: str, hvp_name: str) -> EmittedFunction:
    fn = built.function
    plan = fn.exact_derivative_plan(input_name, "v", "lambda")
    jvp = plan["jvp"]
    hvp = plan["hvp"]
    jac_sparsity = list(plan["jacobian_sparsity"])
    hes_sparsity = list(plan["hessian_sparsity"])
    jac_colors = list(plan["jacobian_colors"])
    hes_colors = list(plan["hessian_colors"])
    report = dict(built.report)
    report.update(
        {
            "jacobian_nnz": len(jac_sparsity),
            "jacobian_colors": int(plan["jacobian_color_count"]),
            "hessian_nnz": len(hes_sparsity),
            "hessian_colors": int(plan["hessian_color_count"]),
        }
    )
    value, jvp_code, hvp_code, jac_code, hes_code = dedupe_ad_vector_helpers([
        fn.to_c(value_name),
        jvp.to_c(jvp_name),
        hvp.to_c(hvp_name),
        fn.to_sparse_jacobian_c(input_name, jac_sparsity, f"{jvp_name}_sparse"),
        hvp.to_sparse_hessian_c("v", hes_sparsity, f"{hvp_name}_sparse"),
    ])
    return EmittedFunction(
        value=value,
        jvp=jvp_code,
        hvp=hvp_code,
        jac=jac_code,
        hes=hes_code,
        jac_sparsity=jac_sparsity,
        hes_sparsity=hes_sparsity,
        jac_colors=jac_colors,
        hes_colors=hes_colors,
        report=report,
    )


def emit_function_from_builder(
    builder: object,
    input_name: str,
    input_vector: object,
    outputs: list[object],
    param_vectors: dict[str, object],
    value_name: str,
    jvp_name: str,
    hvp_name: str,
) -> EmittedFunction:
    built = build_graph_function_from_builder(builder, input_name, input_vector, outputs, param_vectors)
    return _emit_built_function(built, input_name, value_name, jvp_name, hvp_name)


def emit_function(
    input_name: str,
    input_size: int,
    outputs: list[object],
    params: list[tuple[str, int]],
    value_name: str,
    jvp_name: str,
    hvp_name: str,
) -> EmittedFunction:
    built = build_graph_function(input_name, input_size, outputs, params)
    return _emit_built_function(built, input_name, value_name, jvp_name, hvp_name)


def emit_value_function(
    input_name: str,
    input_size: int,
    outputs: list[object],
    params: list[tuple[str, int]],
    value_name: str,
) -> tuple[str, list[tuple[int, int]]]:
    fn = build_graph_function(input_name, input_size, outputs, params).function
    return fn.to_c(value_name), fn.jacobian_sparsity(input_name)
