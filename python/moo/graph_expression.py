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
from .expressions import GraphBuildContext


@dataclass
class GraphExpression:
    value: object | None
    vector: object | None
    source: str


class GraphExpressionEmitter:
    """Graph emission helper for Python modeling expressions.

    Expression nodes emit native AD graph handles through their own ``to_graph``
    methods. This class owns only the shared builder state and small helpers
    needed while Python is still a lightweight modeling layer.
    """

    def __init__(self, builder: ad.GraphFunctionBuilder, variables: dict[str, object]):
        self.builder = builder
        self.ad = ad
        self.variables = variables
        self.context = GraphBuildContext(builder, variables)

    def output(self, output: object) -> GraphExpression:
        if output is None:
            return GraphExpression(None, None, "empty")
        if hasattr(output, "build_graph_vector"):
            return GraphExpression(None, output.build_graph_vector(self.context), "native_backend_vector")
        if hasattr(output, "build_graph"):
            return GraphExpression(output.build_graph(self.context), None, "native_backend_scalar")
        if isinstance(output, (int, float)):
            return GraphExpression(self.builder.constant(float(output)), None, "native_constant")
        if isinstance(output, str):
            raise TypeError("String AD emission has been removed; pass Expr objects instead")
        raise TypeError(f"Cannot emit AD for output of type {type(output)!r}")

    def balanced_sum(self, terms: list[object]):
        if not terms:
            return self.builder.constant(0.0)
        layer = terms
        while len(layer) > 1:
            nxt = []
            for i in range(0, len(layer), 2):
                nxt.append(layer[i] + layer[i + 1] if i + 1 < len(layer) else layer[i])
            layer = nxt
        return layer[0]
