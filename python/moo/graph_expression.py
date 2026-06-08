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

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from . import ad


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
        self._structured_cache: dict[object, object] = {}
        self._vector_cache: dict[object, object] = {}

    def output(self, output: object) -> GraphExpression:
        if output is None:
            return GraphExpression(None, None, "empty")
        vector_node = getattr(output, "vector_node", None)
        if vector_node is not None:
            return GraphExpression(None, self.vector_node(vector_node), "native_vector_expr")
        node = getattr(output, "lfg_node", None)
        if node is not None:
            return GraphExpression(self.node(node), None, "native_expr")
        if isinstance(output, (int, float)):
            return GraphExpression(self.builder.constant(float(output)), None, "native_constant")
        if isinstance(output, str):
            raise TypeError("String AD emission has been removed; pass Expr objects instead")
        raise TypeError(f"Cannot emit AD for output of type {type(output)!r}")

    def node(self, node: object):
        if node is None:
            return None
        if hasattr(node, "to_graph"):
            return node.to_graph(self)
        raise TypeError(f"Unsupported AD expression node type {type(node)!r}")

    def vector_node(self, node: object):
        cached = self._vector_cache.get(node)
        if cached is not None:
            return cached
        if hasattr(node, "to_graph_vector"):
            value = node.to_graph_vector(self)
            self._vector_cache[node] = value
            return value
        raise TypeError(f"Unsupported AD vector expression node type {type(node)!r}")

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

    def dense_matvec_row(self, node: object):
        cache_key = ("dense_matvec", int(node.origin))
        if cache_key not in self._structured_cache:
            rhs = [term.to_graph(self) for term in node.rhs]
            self._structured_cache[cache_key] = self.builder.dense_matvec_values([list(row) for row in node.matrix_values], rhs)
        return self._structured_cache[cache_key][int(node.row_index)]

    def sparse_matvec_row(self, node: object):
        cache_key = ("sparse_matvec", int(node.origin))
        if cache_key not in self._structured_cache:
            rhs = [term.to_graph(self) for term in node.rhs]
            self._structured_cache[cache_key] = self.builder.sparse_matvec_values(list(node.rows), list(node.cols), list(node.vals), tuple(node.shape), rhs)
        return self._structured_cache[cache_key][int(node.row_index)]

    def kron_eye_matvec_row(self, node: object):
        cache_key = ("kron_eye_matvec", int(node.origin))
        if cache_key not in self._structured_cache:
            rhs = [term.to_graph(self) for term in node.rhs]
            self._structured_cache[cache_key] = self.builder.kron_eye_matvec_values([list(row) for row in node.base_values], int(node.eye_size), rhs)
        return self._structured_cache[cache_key][int(node.row_index)]


VariableMapper = Callable[[str, int], object | None]


def remap_expression_node(node: object | None, variable_mapper: VariableMapper) -> object | None:
    if node is None:
        return None
    if hasattr(node, "remap"):
        return node.remap(variable_mapper)
    raise TypeError(f"Unsupported expression node type {type(node)!r}")
