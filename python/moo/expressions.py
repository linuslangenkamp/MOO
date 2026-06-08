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
from collections.abc import Iterable
from typing import Callable


_MISSING = object()


class ExprNode:
    def to_text(self) -> str | None:
        raise NotImplementedError

    def to_graph(self, compiler):
        raise NotImplementedError

    def remap(self, variable_mapper: Callable[[str, int], object | None]):
        raise NotImplementedError


class VectorNode:
    def to_graph_vector(self, compiler):
        raise NotImplementedError

    def remap(self, variable_mapper: Callable[[str, int], object | None]):
        raise NotImplementedError


@dataclass(frozen=True)
class ValuesVectorNode(VectorNode):
    values: tuple[ExprNode, ...]

    def to_graph_vector(self, compiler):
        return compiler.builder.vector([compiler.node(value) for value in self.values])

    def remap(self, variable_mapper):
        return ValuesVectorNode(tuple(_remap_node(value, variable_mapper) for value in self.values))


@dataclass(frozen=True)
class BinaryVectorNode(VectorNode):
    op: str
    lhs: VectorNode
    rhs: VectorNode

    def to_graph_vector(self, compiler):
        lhs = compiler.vector_node(self.lhs)
        rhs = compiler.vector_node(self.rhs)
        if self.op == "+":
            return compiler.builder.vector_add(lhs, rhs)
        if self.op == "-":
            return compiler.builder.vector_sub(lhs, rhs)
        raise ValueError(f"Unsupported vector binary op {self.op!r}")

    def remap(self, variable_mapper):
        return BinaryVectorNode(self.op, self.lhs.remap(variable_mapper), self.rhs.remap(variable_mapper))


@dataclass(frozen=True)
class ScaleVectorNode(VectorNode):
    factor: float
    rhs: VectorNode

    def to_graph_vector(self, compiler):
        return compiler.builder.vector_scale(float(self.factor), compiler.vector_node(self.rhs))

    def remap(self, variable_mapper):
        return ScaleVectorNode(float(self.factor), self.rhs.remap(variable_mapper))


@dataclass(frozen=True)
class SliceVectorNode(VectorNode):
    rhs: VectorNode
    start: int
    length: int
    stride: int = 1

    def to_graph_vector(self, compiler):
        return compiler.builder.slice(compiler.vector_node(self.rhs), int(self.start), int(self.length), int(self.stride))

    def remap(self, variable_mapper):
        return SliceVectorNode(self.rhs.remap(variable_mapper), int(self.start), int(self.length), int(self.stride))


@dataclass(frozen=True)
class ConcatVectorNode(VectorNode):
    lhs: VectorNode
    rhs: VectorNode

    def to_graph_vector(self, compiler):
        return compiler.builder.concat(compiler.vector_node(self.lhs), compiler.vector_node(self.rhs))

    def remap(self, variable_mapper):
        return ConcatVectorNode(self.lhs.remap(variable_mapper), self.rhs.remap(variable_mapper))


@dataclass(frozen=True)
class DenseMatVecVectorNode(VectorNode):
    matrix_values: tuple[tuple[float, ...], ...]
    rhs: tuple[ExprNode, ...]

    def to_graph_vector(self, compiler):
        rhs = [compiler.node(term) for term in self.rhs]
        return compiler.builder.dense_matvec_values([list(row) for row in self.matrix_values], rhs)

    def remap(self, variable_mapper):
        return DenseMatVecVectorNode(
            self.matrix_values,
            tuple(_remap_node(term, variable_mapper) for term in self.rhs),
        )


@dataclass(frozen=True)
class SparseMatVecVectorNode(VectorNode):
    rows: tuple[int, ...]
    cols: tuple[int, ...]
    vals: tuple[float, ...]
    shape: tuple[int, int]
    rhs: tuple[ExprNode, ...]

    def to_graph_vector(self, compiler):
        rhs = [compiler.node(term) for term in self.rhs]
        return compiler.builder.sparse_matvec_values(list(self.rows), list(self.cols), list(self.vals), tuple(self.shape), rhs)

    def remap(self, variable_mapper):
        return SparseMatVecVectorNode(
            self.rows,
            self.cols,
            self.vals,
            self.shape,
            tuple(_remap_node(term, variable_mapper) for term in self.rhs),
        )


@dataclass(frozen=True)
class KronEyeMatVecVectorNode(VectorNode):
    base_values: tuple[tuple[float, ...], ...]
    eye_size: int
    rhs: tuple[ExprNode, ...]

    def to_graph_vector(self, compiler):
        rhs = [compiler.node(term) for term in self.rhs]
        return compiler.builder.kron_eye_matvec_values([list(row) for row in self.base_values], int(self.eye_size), rhs)

    def remap(self, variable_mapper):
        return KronEyeMatVecVectorNode(
            self.base_values,
            int(self.eye_size),
            tuple(_remap_node(term, variable_mapper) for term in self.rhs),
        )


@dataclass(frozen=True)
class VectorAtNode(ExprNode):
    vector: VectorNode
    index: int

    def to_text(self) -> str | None:
        return None

    def to_graph(self, compiler):
        return compiler.builder.at(compiler.vector_node(self.vector), int(self.index))

    def remap(self, variable_mapper):
        return VectorAtNode(self.vector.remap(variable_mapper), int(self.index))


@dataclass(frozen=True)
class ConstNode(ExprNode):
    value: float

    def to_text(self) -> str:
        return repr(float(self.value))

    def to_graph(self, compiler):
        return compiler.builder.constant(float(self.value))

    def remap(self, variable_mapper):
        return self


@dataclass(frozen=True)
class VarNode(ExprNode):
    group: str
    index: int

    def to_text(self) -> str:
        return f"{self.group}[{int(self.index)}]"

    def to_graph(self, compiler):
        return compiler.variables[self.group][int(self.index)]

    def remap(self, variable_mapper):
        mapped = variable_mapper(self.group, int(self.index))
        return self if mapped is None else mapped


@dataclass(frozen=True)
class SumNode(ExprNode):
    terms: tuple[ExprNode, ...]

    def to_text(self) -> str | None:
        parts = [term.to_text() for term in self.terms]
        if any(part is None for part in parts):
            return None
        return _balanced_join([part for part in parts if part is not None], "+")

    def to_graph(self, compiler):
        return compiler.balanced_sum([term.to_graph(compiler) for term in self.terms])

    def remap(self, variable_mapper):
        return SumNode(tuple(_remap_node(term, variable_mapper) for term in self.terms))


@dataclass(frozen=True)
class BinNode(ExprNode):
    op: str
    lhs: ExprNode
    rhs: ExprNode

    def to_text(self) -> str | None:
        lhs = self.lhs.to_text()
        rhs = self.rhs.to_text()
        if lhs is None or rhs is None:
            return None
        return f"({lhs} {self.op} {rhs})"

    def to_graph(self, compiler):
        lhs = self.lhs.to_graph(compiler)
        rhs = self.rhs.to_graph(compiler)
        if self.op == "+":
            return lhs + rhs
        if self.op == "-":
            return lhs - rhs
        if self.op == "*":
            return lhs * rhs
        if self.op == "/":
            return lhs / rhs
        raise ValueError(f"Unsupported AD binary op {self.op!r}")

    def remap(self, variable_mapper):
        return BinNode(self.op, _remap_node(self.lhs, variable_mapper), _remap_node(self.rhs, variable_mapper))


@dataclass(frozen=True)
class UnaryNode(ExprNode):
    name: str
    arg: ExprNode

    def to_text(self) -> str | None:
        arg = self.arg.to_text()
        if arg is None:
            return None
        return f"{self.name}({arg})" if self.name != "neg" else f"(-{arg})"

    def to_graph(self, compiler):
        arg = self.arg.to_graph(compiler)
        if self.name == "neg":
            return -arg
        try:
            fn = getattr(compiler.ad, self.name)
        except AttributeError as exc:
            raise ValueError(f"Unsupported AD unary op {self.name!r}") from exc
        return fn(arg)

    def remap(self, variable_mapper):
        return UnaryNode(self.name, _remap_node(self.arg, variable_mapper))


@dataclass(frozen=True)
class PowConstNode(ExprNode):
    arg: ExprNode
    power: float

    def to_text(self) -> str | None:
        arg = self.arg.to_text()
        if arg is None:
            return None
        return f"pow_const({arg}, {repr(float(self.power))})"

    def to_graph(self, compiler):
        return compiler.ad.pow_const(self.arg.to_graph(compiler), float(self.power))

    def remap(self, variable_mapper):
        return PowConstNode(_remap_node(self.arg, variable_mapper), self.power)


def _remap_node(node: object | None, variable_mapper):
    if node is None:
        return None
    if isinstance(node, ExprNode):
        return node.remap(variable_mapper)
    raise TypeError(f"Unsupported expression node type {type(node)!r}")


class Expr:
    def __init__(self, lfg_node: object | None = None, mr_node: object = _MISSING):
        self.lfg_node = lfg_node
        self.mr_node = lfg_node if mr_node is _MISSING else mr_node

    @property
    def lfg(self) -> str | None:
        return _node_to_text(self.lfg_node)

    @property
    def mr(self) -> str | None:
        return _node_to_text(self.mr_node)

    @staticmethod
    def const(value: float | int) -> "Expr":
        node = ConstNode(float(value))
        return Expr(node, node)

    def is_zero(self) -> bool:
        return (
            self.lfg_node is not None
            and self.mr_node is not None
            and isinstance(self.lfg_node, ConstNode)
            and isinstance(self.mr_node, ConstNode)
            and float(self.lfg_node.value) == 0.0
            and float(self.mr_node.value) == 0.0
        )

    @staticmethod
    def variable(group: str, idx: int, mr_node: object = _MISSING) -> "Expr":
        node = VarNode(group, int(idx))
        return Expr(node, node if mr_node is _MISSING else mr_node)

    def _binary(self, other: object, op: str) -> "Expr":
        rhs = as_expr(other)
        if op == "+":
            return SumExpr.from_terms([self, rhs])
        lfg_node = BinNode(op, self.lfg_node, rhs.lfg_node) if self.lfg_node is not None and rhs.lfg_node is not None else None
        mr_node = BinNode(op, self.mr_node, rhs.mr_node) if self.mr_node is not None and rhs.mr_node is not None else None
        return Expr(lfg_node, mr_node)

    def __add__(self, other: object) -> "Expr":
        return self._binary(other, "+")

    def __radd__(self, other: object) -> "Expr":
        return as_expr(other)._binary(self, "+")

    def __sub__(self, other: object) -> "Expr":
        return self._binary(other, "-")

    def __rsub__(self, other: object) -> "Expr":
        return as_expr(other)._binary(self, "-")

    def __mul__(self, other: object) -> "Expr":
        return self._binary(other, "*")

    def __rmul__(self, other: object) -> "Expr":
        return as_expr(other)._binary(self, "*")

    def __truediv__(self, other: object) -> "Expr":
        return self._binary(other, "/")

    def __rtruediv__(self, other: object) -> "Expr":
        return as_expr(other)._binary(self, "/")

    def __neg__(self) -> "Expr":
        return Expr(
            UnaryNode("neg", self.lfg_node) if self.lfg_node is not None else None,
            UnaryNode("neg", self.mr_node) if self.mr_node is not None else None,
        )

    def __pow__(self, power: float | int) -> "Expr":
        return pow_const(self, power)


def _balanced_join(parts: list[str], op: str) -> str:
    if not parts:
        return "0.0"
    layer = parts
    while len(layer) > 1:
        nxt = []
        for i in range(0, len(layer), 2):
            if i + 1 < len(layer):
                nxt.append(f"({layer[i]} {op} {layer[i + 1]})")
            else:
                nxt.append(layer[i])
        layer = nxt
    return layer[0]


def _node_to_text(node: object | None) -> str | None:
    if node is None:
        return None
    if isinstance(node, ExprNode):
        return node.to_text()
    return None


class SumExpr(Expr):
    def __init__(self, terms: list[Expr]):
        self.terms = terms

    @staticmethod
    def from_terms(terms: list[Expr]) -> Expr:
        flat: list[Expr] = []
        for term in terms:
            expr = as_expr(term)
            if isinstance(expr, SumExpr):
                flat.extend(expr.terms)
            elif expr.is_zero():
                continue
            else:
                flat.append(expr)
        if not flat:
            return Expr.const(0.0)
        if len(flat) == 1:
            return flat[0]
        return SumExpr(flat)

    @property
    def lfg(self) -> str | None:
        parts = [term.lfg for term in self.terms]
        if any(part is None for part in parts):
            return None
        return _balanced_join([part for part in parts if part is not None], "+")

    @property
    def mr(self) -> str | None:
        parts = [term.mr for term in self.terms]
        if any(part is None for part in parts):
            return None
        return _balanced_join([part for part in parts if part is not None], "+")

    @property
    def lfg_node(self):
        nodes = [term.lfg_node for term in self.terms]
        if any(node is None for node in nodes):
            return None
        return SumNode(tuple(nodes))

    @property
    def mr_node(self):
        nodes = [term.mr_node for term in self.terms]
        if any(node is None for node in nodes):
            return None
        return SumNode(tuple(nodes))


def as_expr(value: object) -> Expr:
    if isinstance(value, Expr):
        return value
    if isinstance(value, (int, float)):
        return Expr.const(value)
    raise TypeError(f"Cannot convert {type(value)!r} to Expr")


def _unary(name: str, value: object) -> Expr:
    expr = as_expr(value)
    return Expr(
        UnaryNode(name, expr.lfg_node) if expr.lfg_node is not None else None,
        UnaryNode(name, expr.mr_node) if expr.mr_node is not None else None,
    )


def sin(value: object) -> Expr:
    return _unary("sin", value)


def cos(value: object) -> Expr:
    return _unary("cos", value)


def tan(value: object) -> Expr:
    return _unary("tan", value)


def exp(value: object) -> Expr:
    return _unary("exp", value)


def log(value: object) -> Expr:
    return _unary("log", value)


def pow_const(value: object, power: float | int) -> Expr:
    expr = as_expr(value)
    return Expr(
        PowConstNode(expr.lfg_node, float(power)) if expr.lfg_node is not None else None,
        PowConstNode(expr.mr_node, float(power)) if expr.mr_node is not None else None,
    )


def sum_expr(values: Iterable[object]) -> Expr:
    return SumExpr.from_terms([as_expr(value) for value in values])


def _vector_at_expr(vector_node: VectorNode | None, row_index: int, fallback: Expr) -> Expr:
    if vector_node is None:
        return fallback
    node = VectorAtNode(vector_node, int(row_index))
    return Expr(node, node)


def _dense_matvec_element_expr(matrix_values: list[list[float]], row_index: int, rhs: list[Expr], vector_node: VectorNode | None) -> Expr:
    fallback = sum_expr(weight * value for weight, value in zip(matrix_values[row_index], rhs))
    return _vector_at_expr(vector_node, row_index, fallback)


def _sparse_matvec_element_expr(rows: list[int], cols: list[int], vals: list[float], row_index: int, rhs: list[Expr], vector_node: VectorNode | None) -> Expr:
    buckets: list[tuple[float, Expr]] = [(val, rhs[col]) for row, col, val in zip(rows, cols, vals) if int(row) == int(row_index)]
    fallback = sum_expr(val * value for val, value in buckets)
    return _vector_at_expr(vector_node, row_index, fallback)


def _kron_eye_matvec_element_expr(base_values: list[list[float]], eye_size: int, row_index: int, rhs: list[Expr], vector_node: VectorNode | None) -> Expr:
    eye = int(eye_size)
    base_row = int(row_index) // eye
    inner = int(row_index) % eye
    fallback = sum_expr(weight * rhs[col * eye + inner] for col, weight in enumerate(base_values[base_row]))
    return _vector_at_expr(vector_node, row_index, fallback)


def _is_vector_like(value: object) -> bool:
    return isinstance(value, (VectorExpr, BlockVectorExpr)) or bool(getattr(value, "__moo_vector__", False))


def _as_vector_node(value: object) -> VectorNode | None:
    if isinstance(value, VectorExpr):
        return value.vector_node
    if isinstance(value, BlockVectorExpr):
        return value.flatten().vector_node
    if getattr(value, "__moo_vector__", False):
        return vec(value).vector_node
    if isinstance(value, (list, tuple)):
        return vec(value).vector_node
    return None


def _expr_node_for_vector_value(value: Expr) -> ExprNode | None:
    return value.lfg_node if value.lfg_node is not None else value.mr_node


def _values_vector_node(values: list[Expr]) -> VectorNode | None:
    nodes = [_expr_node_for_vector_value(value) for value in values]
    if any(node is None for node in nodes):
        return None
    return ValuesVectorNode(tuple(node for node in nodes if node is not None))


class VectorExpr:
    def __init__(self, values: Iterable[object], vector_node: VectorNode | None = None):
        self.values = [as_expr(value) for value in values]
        self.vector_node = vector_node if vector_node is not None else _values_vector_node(self.values)

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            start, stop, stride = idx.indices(len(self.values))
            values = self.values[idx]
            length = len(values)
            vector_node = SliceVectorNode(self.vector_node, start, length, stride) if self.vector_node is not None else None
            return VectorExpr(values, vector_node)
        return self.values[idx]

    def __iter__(self):
        return iter(self.values)

    def block(self, index: int, size: int) -> "VectorExpr":
        if size <= 0:
            raise ValueError("block size must be positive")
        start = int(index) * size
        stop = start + size
        if start < 0 or stop > len(self.values):
            raise IndexError(index)
        vector_node = SliceVectorNode(self.vector_node, start, size, 1) if self.vector_node is not None else None
        return VectorExpr(self.values[start:stop], vector_node)

    def blocks(self, size: int) -> "BlockVectorExpr":
        if size <= 0:
            raise ValueError("block size must be positive")
        if len(self.values) % size != 0:
            raise ValueError("vector length must be divisible by block size")
        return BlockVectorExpr([self.block(i, size) for i in range(len(self.values) // size)])

    def _binary(self, other: object, op: str, reverse: bool = False) -> "VectorExpr":
        if _is_vector_like(other):
            rhs_values = _as_vector_values(other)
            if len(rhs_values) != len(self.values):
                raise ValueError("vector sizes must match")
            lhs_node = self.vector_node
            rhs_node = _as_vector_node(other)
            vector_node = None
            if lhs_node is not None and rhs_node is not None and op in {"+", "-"}:
                vector_node = BinaryVectorNode(op, rhs_node, lhs_node) if reverse else BinaryVectorNode(op, lhs_node, rhs_node)
            if reverse:
                return VectorExpr([as_expr(rhs)._binary(lhs, op) for lhs, rhs in zip(self.values, rhs_values)], vector_node)
            return VectorExpr([lhs._binary(rhs, op) for lhs, rhs in zip(self.values, rhs_values)], vector_node)
        if reverse:
            return VectorExpr([as_expr(other)._binary(lhs, op) for lhs in self.values])
        return VectorExpr([lhs._binary(other, op) for lhs in self.values])

    def __add__(self, other: object) -> "VectorExpr":
        return self._binary(other, "+")

    def __radd__(self, other: object) -> "VectorExpr":
        return self._binary(other, "+", reverse=True)

    def __sub__(self, other: object) -> "VectorExpr":
        return self._binary(other, "-")

    def __rsub__(self, other: object) -> "VectorExpr":
        return self._binary(other, "-", reverse=True)

    def __mul__(self, other: object) -> "VectorExpr":
        if _is_vector_like(other):
            return self._binary(other, "*")
        vector_node = ScaleVectorNode(float(other), self.vector_node) if self.vector_node is not None and isinstance(other, (int, float)) else None
        return VectorExpr([value * other for value in self.values], vector_node)

    def __rmul__(self, other: object) -> "VectorExpr":
        return self.__mul__(other)

    def __truediv__(self, other: object) -> "VectorExpr":
        vector_node = ScaleVectorNode(1.0 / float(other), self.vector_node) if self.vector_node is not None and isinstance(other, (int, float)) else None
        return VectorExpr([value / other for value in self.values], vector_node)

    def __neg__(self) -> "VectorExpr":
        vector_node = ScaleVectorNode(-1.0, self.vector_node) if self.vector_node is not None else None
        return VectorExpr([-value for value in self.values], vector_node)


class BlockVectorExpr:
    def __init__(self, blocks: Iterable[object]):
        values = []
        block_size: int | None = None
        for block in blocks:
            vector_block = vec(block)
            if block_size is None:
                block_size = len(vector_block)
            elif len(vector_block) != block_size:
                raise ValueError("all blocks must have equal length")
            values.append(vector_block)
        self.values = values
        self.block_size = block_size or 0

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return BlockVectorExpr(self.values[idx])
        return self.values[idx]

    def __iter__(self):
        return iter(self.values)

    def flatten(self) -> VectorExpr:
        values = [value for block in self.values for value in block]
        node = None
        block_nodes = [block.vector_node for block in self.values]
        if all(block_node is not None for block_node in block_nodes):
            node = block_nodes[0]
            for block_node in block_nodes[1:]:
                node = ConcatVectorNode(node, block_node)
        return VectorExpr(values, node)

    def _binary(self, other: object, op: str, reverse: bool = False) -> "BlockVectorExpr":
        if isinstance(other, BlockVectorExpr):
            if len(other) != len(self) or other.block_size != self.block_size:
                raise ValueError("block vector sizes must match")
            if reverse:
                return BlockVectorExpr([rhs._binary(lhs, op) for lhs, rhs in zip(self.values, other.values)])
            return BlockVectorExpr([lhs._binary(rhs, op) for lhs, rhs in zip(self.values, other.values)])
        if _is_vector_like(other) or isinstance(other, (list, tuple)):
            return self._binary(vec(other).blocks(self.block_size), op, reverse=reverse)
        if reverse:
            return BlockVectorExpr([
                VectorExpr([as_expr(other)._binary(value, op) for value in block])
                for block in self.values
            ])
        return BlockVectorExpr([block._binary(other, op) for block in self.values])

    def __add__(self, other: object) -> "BlockVectorExpr":
        return self._binary(other, "+")

    def __radd__(self, other: object) -> "BlockVectorExpr":
        return self._binary(other, "+", reverse=True)

    def __sub__(self, other: object) -> "BlockVectorExpr":
        return self._binary(other, "-")

    def __rsub__(self, other: object) -> "BlockVectorExpr":
        return self._binary(other, "-", reverse=True)

    def __mul__(self, other: object) -> "BlockVectorExpr":
        return BlockVectorExpr([block * other for block in self.values])

    def __rmul__(self, other: object) -> "BlockVectorExpr":
        return self.__mul__(other)

    def __truediv__(self, other: object) -> "BlockVectorExpr":
        return BlockVectorExpr([block / other for block in self.values])

    def __neg__(self) -> "BlockVectorExpr":
        return BlockVectorExpr([-block for block in self.values])


class ExprMatrix:
    def __init__(self, values: object, vector: bool = False):
        raw = _tolist(values)
        if vector:
            self.values = [float(value) for value in raw]
            self.is_vector = True
            return
        self.values = [[float(value) for value in row] for row in raw]
        self.is_vector = False
        if self.values:
            width = len(self.values[0])
            if any(len(row) != width for row in self.values):
                raise ValueError("matrix rows must have equal length")

    def __matmul__(self, other: object):
        if isinstance(other, BlockVectorExpr):
            if self.is_vector:
                if len(self.values) != len(other):
                    raise ValueError("dot dimensions do not match")
                return VectorExpr([
                    sum_expr(weight * block[row] for weight, block in zip(self.values, other))
                    for row in range(other.block_size)
                ])
            if self.values and len(self.values[0]) != len(other):
                raise ValueError("matrix/block-vector dimensions do not match")
            return BlockVectorExpr([
                VectorExpr([
                    sum_expr(weight * block[row] for weight, block in zip(matrix_row, other))
                    for row in range(other.block_size)
                ])
                for matrix_row in self.values
            ])
        rhs = _as_vector_values(other)
        if self.is_vector:
            if len(self.values) != len(rhs):
                raise ValueError("dot dimensions do not match")
            return sum_expr(weight * value for weight, value in zip(self.values, rhs))
        if self.values and len(self.values[0]) != len(rhs):
            raise ValueError("matrix/vector dimensions do not match")
        frozen_matrix = tuple(tuple(float(value) for value in row) for row in self.values)
        rhs_nodes = [_expr_node_for_vector_value(value) for value in rhs]
        vector_node = DenseMatVecVectorNode(frozen_matrix, tuple(rhs_nodes)) if all(node is not None for node in rhs_nodes) else None
        return VectorExpr(
            [_dense_matvec_element_expr(self.values, row_index, rhs, vector_node) for row_index, _row in enumerate(self.values)],
            vector_node,
        )

    def otimes_eye(self, size: int) -> "KroneckerEyeMatrix":
        if self.is_vector:
            raise ValueError("otimes_eye is only defined for matrices")
        if size <= 0:
            raise ValueError("identity size must be positive")
        return KroneckerEyeMatrix(self, size)

    def kron_eye(self, size: int) -> "KroneckerEyeMatrix":
        return self.otimes_eye(size)


class SparseMatrixExpr:
    def __init__(self, rows: Iterable[int], cols: Iterable[int], vals: Iterable[float], shape: tuple[int, int]):
        self.rows = [int(row) for row in rows]
        self.cols = [int(col) for col in cols]
        self.vals = [float(val) for val in vals]
        self.shape = (int(shape[0]), int(shape[1]))
        self.is_vector = False
        if not (len(self.rows) == len(self.cols) == len(self.vals)):
            raise ValueError("sparse rows, cols, and vals must have equal length")

    def __matmul__(self, other: object) -> VectorExpr:
        rhs = _as_vector_values(other)
        if len(rhs) != self.shape[1]:
            raise ValueError("sparse matrix/vector dimensions do not match")
        for row, col in zip(self.rows, self.cols):
            if row < 0 or row >= self.shape[0] or col < 0 or col >= self.shape[1]:
                raise ValueError("sparse matrix entry out of shape bounds")
        rhs_nodes = [_expr_node_for_vector_value(value) for value in rhs]
        vector_node = SparseMatVecVectorNode(
            tuple(self.rows),
            tuple(self.cols),
            tuple(self.vals),
            self.shape,
            tuple(rhs_nodes),
        ) if all(node is not None for node in rhs_nodes) else None
        return VectorExpr(
            [_sparse_matvec_element_expr(self.rows, self.cols, self.vals, row_index, rhs, vector_node) for row_index in range(self.shape[0])],
            vector_node,
        )


class KroneckerEyeMatrix:
    def __init__(self, base: ExprMatrix, eye_size: int):
        if base.is_vector:
            raise ValueError("KroneckerEyeMatrix requires a matrix base")
        self.base = base
        self.eye_size = eye_size
        self.is_vector = False
        self.values = None

    def __matmul__(self, other: object):
        rhs = vec(other)
        expected = len(self.base.values[0]) * self.eye_size if self.base.values else 0
        if len(rhs) != expected:
            raise ValueError("kron-eye matrix/vector dimensions do not match")
        frozen_base = tuple(tuple(float(value) for value in row) for row in self.base.values)
        rhs_nodes = [_expr_node_for_vector_value(value) for value in rhs.values]
        vector_node = KronEyeMatVecVectorNode(frozen_base, self.eye_size, tuple(rhs_nodes)) if all(node is not None for node in rhs_nodes) else None
        return VectorExpr(
            [
                _kron_eye_matvec_element_expr(self.base.values, self.eye_size, row_index, rhs.values, vector_node)
                for row_index in range(len(self.base.values) * self.eye_size)
            ],
            vector_node,
        )

    def kron_eye(self, size: int) -> "KroneckerEyeMatrix":
        if size != self.eye_size:
            raise ValueError("nested kron_eye is not supported")
        return self


def _tolist(values: object):
    if hasattr(values, "tolist"):
        return values.tolist()
    return values


def _as_vector_values(values: object) -> list[Expr]:
    if isinstance(values, VectorExpr):
        return list(values.values)
    if isinstance(values, BlockVectorExpr):
        return _as_vector_values(values.flatten())
    if getattr(values, "__moo_vector__", False):
        return [values[i] for i in range(len(values))]
    raw = _tolist(values)
    if isinstance(raw, (list, tuple)):
        flattened = []
        for value in raw:
            if _is_vector_like(value):
                flattened.extend(_as_vector_values(value))
            elif isinstance(value, (list, tuple)):
                flattened.extend(_as_vector_values(value))
            else:
                flattened.append(as_expr(value))
        return flattened
    raise TypeError(f"Cannot convert {type(values)!r} to expression vector")


def vec(values: object) -> VectorExpr:
    if isinstance(values, VectorExpr):
        return values
    return VectorExpr(_as_vector_values(values))


def blockvec(values: object) -> BlockVectorExpr:
    return BlockVectorExpr(_tolist(values))


def vector(values: object) -> ExprMatrix:
    return ExprMatrix(values, vector=True)


def matrix(values: object) -> ExprMatrix:
    return ExprMatrix(values, vector=False)


def sparse_matrix(rows: Iterable[int], cols: Iterable[int], vals: Iterable[float], shape: tuple[int, int]) -> SparseMatrixExpr:
    return SparseMatrixExpr(rows, cols, vals, shape)


def dot(lhs: object, rhs: object) -> Expr:
    return vector(lhs) @ vec(rhs)
