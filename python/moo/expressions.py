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

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Callable


GraphScalarBuilder = Callable[[object], object]
GraphVectorBuilder = Callable[[object], object]


class GraphBuildContext:
    def __init__(self, builder: object, variables: dict[str, object]):
        self.builder = builder
        self.variables = variables

    def scalar(self, group: str, index: int):
        try:
            return self.variables[group][int(index)]
        except KeyError as exc:
            raise KeyError(f"unknown graph variable group {group!r}") from exc

    def vector(self, values: list[object]):
        return self.builder.vector(values)

    def constant(self, value: float):
        return self.builder.constant(float(value))


class RemapContext:
    def __init__(self, parent: GraphBuildContext, mapper: Callable[[str, int], tuple[str, int] | None]):
        self.parent = parent
        self.builder = parent.builder
        self.mapper = mapper

    def scalar(self, group: str, index: int):
        mapped = self.mapper(group, int(index))
        if mapped is None:
            return self.parent.scalar(group, int(index))
        return self.parent.scalar(mapped[0], mapped[1])

    def vector(self, values: list[object]):
        return self.builder.vector(values)

    def constant(self, value: float):
        return self.builder.constant(float(value))


@dataclass(frozen=True)
class Symbol:
    group: str
    index: int


class Expr:
    def __init__(
        self,
        build: GraphScalarBuilder | None = None,
        *,
        mr_build: GraphScalarBuilder | None = None,
        symbol: Symbol | None = None,
        mr_symbol: Symbol | None = None,
        const_value: float | None = None,
        graph_scalar: object | None = None,
    ):
        self._build = build
        self._mr_build = build if mr_build is None else mr_build
        self.symbol = symbol
        self.mr_symbol = symbol if mr_symbol is None else mr_symbol
        self.const_value = const_value
        self.graph_scalar = graph_scalar

    @staticmethod
    def const(value: float | int) -> "Expr":
        value = float(value)
        return Expr(lambda ctx, value=value: ctx.constant(value), const_value=value)

    @staticmethod
    def variable(group: str, idx: int, *, mr_group: str | None = None, mr_idx: int | None = None, graph_scalar: object | None = None) -> "Expr":
        idx = int(idx)
        symbol = Symbol(group, idx)
        mr_symbol = Symbol(mr_group, int(mr_idx if mr_idx is not None else idx)) if mr_group is not None else symbol
        return Expr(
            lambda ctx, group=group, idx=idx: ctx.scalar(group, idx),
            mr_build=lambda ctx, group=mr_symbol.group, idx=mr_symbol.index: ctx.scalar(group, idx),
            symbol=symbol,
            mr_symbol=mr_symbol,
            graph_scalar=graph_scalar,
        )

    def build_graph(self, ctx: GraphBuildContext):
        if self.graph_scalar is not None:
            return self.graph_scalar
        if self._build is None:
            raise ValueError("expression has no graph representation")
        return self._build(ctx)

    def build_mr_graph(self, ctx: GraphBuildContext):
        if self._mr_build is None:
            raise ValueError("expression has no boundary graph representation")
        return self._mr_build(ctx)

    def remap_mr(self, mapper: Callable[[str, int], tuple[str, int] | None]) -> "Expr":
        return Expr(lambda ctx: self.build_mr_graph(RemapContext(ctx, mapper)))

    def is_zero(self) -> bool:
        return self.const_value == 0.0

    def _binary(self, other: object, op: str) -> "Expr":
        rhs = as_expr(other)

        def apply(lhs_value, rhs_value):
            if op == "+":
                return lhs_value + rhs_value
            if op == "-":
                return lhs_value - rhs_value
            if op == "*":
                return lhs_value * rhs_value
            if op == "/":
                return lhs_value / rhs_value
            raise ValueError(f"unsupported binary operation {op!r}")

        graph_scalar = None
        if self.graph_scalar is not None and rhs.graph_scalar is not None:
            graph_scalar = apply(self.graph_scalar, rhs.graph_scalar)

        const_value = None
        if self.const_value is not None and rhs.const_value is not None:
            const_value = float(apply(self.const_value, rhs.const_value))

        return Expr(
            lambda ctx: apply(self.build_graph(ctx), rhs.build_graph(ctx)),
            mr_build=lambda ctx: apply(self.build_mr_graph(ctx), rhs.build_mr_graph(ctx)),
            const_value=const_value,
            graph_scalar=graph_scalar,
        )

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
        graph_scalar = -self.graph_scalar if self.graph_scalar is not None else None
        const_value = -self.const_value if self.const_value is not None else None
        return Expr(
            lambda ctx: -self.build_graph(ctx),
            mr_build=lambda ctx: -self.build_mr_graph(ctx),
            const_value=const_value,
            graph_scalar=graph_scalar,
        )

    def __pow__(self, power: float | int) -> "Expr":
        return pow_const(self, power)


class SumExpr(Expr):
    def __init__(self, terms: list[Expr]):
        self.terms = terms
        super().__init__(
            lambda ctx: _balanced_backend_sum(ctx, [term.build_graph(ctx) for term in self.terms]),
            mr_build=lambda ctx: _balanced_backend_sum(ctx, [term.build_mr_graph(ctx) for term in self.terms]),
            const_value=sum(term.const_value for term in self.terms) if all(term.const_value is not None for term in self.terms) else None,
            graph_scalar=_sum_graph_scalars(self.terms),
        )

    @staticmethod
    def from_terms(terms: Iterable[object]) -> Expr:
        flat: list[Expr] = []
        for term in terms:
            expr = as_expr(term)
            if isinstance(expr, SumExpr):
                flat.extend(expr.terms)
            elif not expr.is_zero():
                flat.append(expr)
        if not flat:
            return Expr.const(0.0)
        if len(flat) == 1:
            return flat[0]
        return SumExpr(flat)


def _sum_graph_scalars(terms: list[Expr]):
    if not terms or any(term.graph_scalar is None for term in terms):
        return None
    out = terms[0].graph_scalar
    for term in terms[1:]:
        out = out + term.graph_scalar
    return out


def _balanced_backend_sum(ctx: GraphBuildContext, terms: list[object]):
    if not terms:
        return ctx.constant(0.0)
    layer = terms
    while len(layer) > 1:
        nxt = []
        for i in range(0, len(layer), 2):
            nxt.append(layer[i] + layer[i + 1] if i + 1 < len(layer) else layer[i])
        layer = nxt
    return layer[0]


def as_expr(value: object) -> Expr:
    if isinstance(value, Expr):
        return value
    if _is_graph_scalar(value):
        return Expr(graph_scalar=value)
    if isinstance(value, (int, float)):
        return Expr.const(value)
    raise TypeError(f"Cannot convert {type(value)!r} to Expr")


def _unary(name: str, value: object) -> Expr:
    expr = as_expr(value)
    graph_scalar = None
    if expr.graph_scalar is not None:
        from . import ad

        graph_scalar = getattr(ad, name)(expr.graph_scalar)
    return Expr(
        lambda ctx: getattr(_ad_module(), name)(expr.build_graph(ctx)),
        mr_build=lambda ctx: getattr(_ad_module(), name)(expr.build_mr_graph(ctx)),
        graph_scalar=graph_scalar,
    )


def _ad_module():
    from . import ad

    return ad


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


def pow_const(value: object, power: float | int) -> Expr | "VectorExpr":
    power = float(power)
    if _is_vector_like(value):
        vector_value = vec(value)
        graph_vector = vector_value.graph_vector.pow_const(power) if vector_value.graph_vector is not None else None
        return VectorExpr(
            [pow_const(item, power) for item in vector_value],
            graph_vector=graph_vector,
            vector_build=lambda ctx, power=power: ctx.builder.vector_pow_const(vector_value.build_graph_vector(ctx), power),
        )
    expr = as_expr(value)
    graph_scalar = None
    if expr.graph_scalar is not None:
        graph_scalar = _ad_module().pow_const(expr.graph_scalar, power)
    const_value = expr.const_value ** power if expr.const_value is not None else None
    return Expr(
        lambda ctx: _ad_module().pow_const(expr.build_graph(ctx), power),
        mr_build=lambda ctx: _ad_module().pow_const(expr.build_mr_graph(ctx), power),
        const_value=const_value,
        graph_scalar=graph_scalar,
    )


def sum_expr(values: Iterable[object]) -> Expr:
    return SumExpr.from_terms(values)


def _is_graph_scalar(value: object) -> bool:
    try:
        from . import ad
    except ImportError:
        return False
    return isinstance(value, ad.GraphScalar)


def _is_graph_vector(value: object) -> bool:
    try:
        from . import ad
    except ImportError:
        return False
    return isinstance(value, ad.GraphVector)


def _as_graph_vector(value: object):
    if isinstance(value, VectorExpr):
        return value.graph_vector
    if isinstance(value, BlockVectorExpr):
        return value.flatten().graph_vector
    if _is_graph_vector(value):
        return value
    return None


class VectorExpr:
    def __init__(self, values: Iterable[object], graph_vector: object | None = None, vector_build: GraphVectorBuilder | None = None):
        self.values = [as_expr(value) for value in values]
        self.graph_vector = graph_vector
        self._vector_build = vector_build
        self._build_cache: dict[int, object] = {}

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            values = self.values[idx]
            vector_build = None
            if self._vector_build is not None:
                start, _stop, stride = idx.indices(len(self.values))
                vector_build = lambda ctx, start=start, length=len(values), stride=stride: ctx.builder.slice(self.build_graph_vector(ctx), start, length, stride)
            return VectorExpr(values, vector_build=vector_build)
        return self.values[idx]

    def __iter__(self):
        return iter(self.values)

    def build_graph_vector(self, ctx: GraphBuildContext):
        if self.graph_vector is not None:
            return self.graph_vector
        cache_key = id(ctx)
        if cache_key in self._build_cache:
            return self._build_cache[cache_key]
        if self._vector_build is not None:
            result = self._vector_build(ctx)
        else:
            result = ctx.vector([value.build_graph(ctx) for value in self.values])
        self._build_cache[cache_key] = result
        return result

    def block(self, index: int, size: int) -> "VectorExpr":
        if size <= 0:
            raise ValueError("block size must be positive")
        start = int(index) * size
        stop = start + size
        if start < 0 or stop > len(self.values):
            raise IndexError(index)
        return self[start:stop]

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
            lhs_graph = self.graph_vector
            rhs_graph = _as_graph_vector(other)
            graph_vector = None
            if lhs_graph is not None and rhs_graph is not None and op in {"+", "-", "*"}:
                if reverse and op == "+":
                    graph_vector = rhs_graph + lhs_graph
                elif reverse and op == "-":
                    graph_vector = rhs_graph - lhs_graph
                elif op == "+":
                    graph_vector = lhs_graph + rhs_graph
                elif op == "*":
                    graph_vector = rhs_graph * lhs_graph if reverse else lhs_graph * rhs_graph
                else:
                    graph_vector = lhs_graph - rhs_graph

            def vector_build(ctx):
                lhs = self.build_graph_vector(ctx)
                rhs = vec(other).build_graph_vector(ctx)
                if reverse:
                    lhs, rhs = rhs, lhs
                if op == "+":
                    return ctx.builder.vector_add(lhs, rhs)
                if op == "-":
                    return ctx.builder.vector_sub(lhs, rhs)
                if op == "*":
                    return ctx.builder.vector_mul(lhs, rhs)
                return ctx.vector([a._binary(b, op).build_graph(ctx) for a, b in zip(vec(other).values, self.values)] if reverse else [a._binary(b, op).build_graph(ctx) for a, b in zip(self.values, rhs_values)])

            values = [as_expr(rhs)._binary(lhs, op) for lhs, rhs in zip(self.values, rhs_values)] if reverse else [lhs._binary(rhs, op) for lhs, rhs in zip(self.values, rhs_values)]
            return VectorExpr(values, graph_vector=graph_vector, vector_build=vector_build if op in {"+", "-", "*"} else None)
        values = [as_expr(other)._binary(lhs, op) for lhs in self.values] if reverse else [lhs._binary(other, op) for lhs in self.values]
        return VectorExpr(values)

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
        factor = float(other)
        graph_vector = self.graph_vector * factor if self.graph_vector is not None else None
        return VectorExpr([value * factor for value in self.values], graph_vector=graph_vector, vector_build=lambda ctx, factor=factor: ctx.builder.vector_scale(factor, self.build_graph_vector(ctx)))

    def __rmul__(self, other: object) -> "VectorExpr":
        return self.__mul__(other)

    def __truediv__(self, other: object) -> "VectorExpr":
        return self.__mul__(1.0 / float(other))

    def __neg__(self) -> "VectorExpr":
        return self.__mul__(-1.0)

    def __pow__(self, power: float | int) -> "VectorExpr":
        return pow_const(self, power)


class BlockVectorExpr:
    def __init__(self, blocks: Iterable[object]):
        self.values = [vec(block) for block in blocks]
        self.block_size = len(self.values[0]) if self.values else 0
        if any(len(block) != self.block_size for block in self.values):
            raise ValueError("all blocks must have equal length")

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
        return VectorExpr(values, vector_build=lambda ctx: _concat_vectors(ctx, [block.build_graph_vector(ctx) for block in self.values]))

    def _binary(self, other: object, op: str, reverse: bool = False) -> "BlockVectorExpr":
        if isinstance(other, BlockVectorExpr):
            if len(other) != len(self) or other.block_size != self.block_size:
                raise ValueError("block vector sizes must match")
            return BlockVectorExpr([rhs._binary(lhs, op) for lhs, rhs in zip(self.values, other.values)] if reverse else [lhs._binary(rhs, op) for lhs, rhs in zip(self.values, other.values)])
        if _is_vector_like(other) or isinstance(other, (list, tuple)):
            return self._binary(vec(other).blocks(self.block_size), op, reverse=reverse)
        return BlockVectorExpr([VectorExpr([as_expr(other)._binary(value, op) for value in block]) for block in self.values] if reverse else [block._binary(other, op) for block in self.values])

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


def _concat_vectors(ctx: GraphBuildContext, vectors: list[object]):
    if not vectors:
        return ctx.vector([])
    out = vectors[0]
    for vector_value in vectors[1:]:
        out = ctx.builder.concat(out, vector_value)
    return out


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
                return VectorExpr([sum_expr(weight * block[row] for weight, block in zip(self.values, other)) for row in range(other.block_size)])
            if self.values and len(self.values[0]) != len(other):
                raise ValueError("matrix/block-vector dimensions do not match")
            return BlockVectorExpr([VectorExpr([sum_expr(weight * block[row] for weight, block in zip(matrix_row, other)) for row in range(other.block_size)]) for matrix_row in self.values])
        rhs = vec(other)
        if self.is_vector:
            if len(self.values) != len(rhs):
                raise ValueError("dot dimensions do not match")
            return sum_expr(weight * value for weight, value in zip(self.values, rhs))
        if self.values and len(self.values[0]) != len(rhs):
            raise ValueError("matrix/vector dimensions do not match")
        graph_rhs = _as_graph_vector(rhs)
        graph_vector = None
        if graph_rhs is not None:
            # GraphVector matvec is available only through the active builder; keep
            # the callable path as the canonical representation.
            graph_vector = None
        return VectorExpr(
            [sum_expr(weight * value for weight, value in zip(row, rhs)) for row in self.values],
            graph_vector=graph_vector,
            vector_build=lambda ctx, matrix_values=self.values, rhs=rhs: ctx.builder.dense_matvec(matrix_values, rhs.build_graph_vector(ctx)),
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
        rhs = vec(other)
        if len(rhs) != self.shape[1]:
            raise ValueError("sparse matrix/vector dimensions do not match")
        for row, col in zip(self.rows, self.cols):
            if row < 0 or row >= self.shape[0] or col < 0 or col >= self.shape[1]:
                raise ValueError("sparse matrix entry out of shape bounds")
        return VectorExpr(
            [_sparse_row_expr(self.rows, self.cols, self.vals, row, rhs.values) for row in range(self.shape[0])],
            vector_build=lambda ctx, self=self, rhs=rhs: ctx.builder.sparse_matvec(self.rows, self.cols, self.vals, self.shape, rhs.build_graph_vector(ctx)),
        )


class KroneckerEyeMatrix:
    def __init__(self, base: ExprMatrix, eye_size: int):
        if base.is_vector:
            raise ValueError("KroneckerEyeMatrix requires a matrix base")
        self.base = base
        self.eye_size = int(eye_size)
        self.is_vector = False
        self.values = None

    def __matmul__(self, other: object) -> VectorExpr:
        rhs = vec(other)
        expected = len(self.base.values[0]) * self.eye_size if self.base.values else 0
        if len(rhs) != expected:
            raise ValueError("kron-eye matrix/vector dimensions do not match")
        rows = len(self.base.values) * self.eye_size
        return VectorExpr(
            [_kron_row_expr(self.base.values, self.eye_size, row, rhs.values) for row in range(rows)],
            vector_build=lambda ctx, self=self, rhs=rhs: ctx.builder.kron_eye_matvec(self.base.values, self.eye_size, rhs.build_graph_vector(ctx)),
        )

    def kron_eye(self, size: int) -> "KroneckerEyeMatrix":
        if int(size) != self.eye_size:
            raise ValueError("nested kron_eye is not supported")
        return self


def _sparse_row_expr(rows: list[int], cols: list[int], vals: list[float], row_index: int, rhs: list[Expr]) -> Expr:
    return sum_expr(val * rhs[col] for row, col, val in zip(rows, cols, vals) if int(row) == int(row_index))


def _kron_row_expr(base_values: list[list[float]], eye_size: int, row_index: int, rhs: list[Expr]) -> Expr:
    base_row = int(row_index) // int(eye_size)
    inner = int(row_index) % int(eye_size)
    return sum_expr(weight * rhs[col * int(eye_size) + inner] for col, weight in enumerate(base_values[base_row]))


def _tolist(values: object):
    if hasattr(values, "tolist"):
        return values.tolist()
    return values


def _is_vector_like(value: object) -> bool:
    return isinstance(value, (VectorExpr, BlockVectorExpr)) or bool(getattr(value, "__moo_vector__", False)) or _is_graph_vector(value)


def _as_vector_values(values: object) -> list[Expr]:
    if isinstance(values, VectorExpr):
        return list(values.values)
    if isinstance(values, BlockVectorExpr):
        return _as_vector_values(values.flatten())
    if _is_graph_vector(values):
        return [as_expr(value) for value in values.values]
    if getattr(values, "__moo_vector__", False):
        return [values[i] for i in range(len(values))]
    raw = _tolist(values)
    if isinstance(raw, (list, tuple)):
        flattened: list[Expr] = []
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
    if _is_graph_vector(values):
        return VectorExpr(values.values, graph_vector=values)
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
