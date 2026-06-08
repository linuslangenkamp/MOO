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

import math
import os
import re
import shutil
import subprocess
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from .callback_codegen import (
    render_hessian_callback_body,
    render_jacobian_callback_body,
    render_local_colored_hes_lines,
    render_local_colored_jac_lines,
    render_local_direct_constraint_jac_lines,
    render_local_direct_hes_lines,
    render_local_direct_objective_jac_lines,
    static_int_array,
)
from .common import select_derivative_callback_mode, parse_sparsity_pairs
from . import paths
from . import ad
from .ad_codegen import emit_function_from_builder
from .expressions import (
    BlockVectorExpr,
    Expr,
    ExprMatrix,
    KroneckerEyeMatrix,
    VectorExpr,
    as_expr,
    blockvec,
    dot,
    matrix,
    sparse_matrix,
    vec,
    vector,
)
from .local_function import InputGroup, LocalGraphFunction
from .model import BaseModel, _arr, _c_bool, _num
from .results import NLPResult, OptimizationRun, read_results


@dataclass
class NLPVar:
    name: str
    lb: float = -math.inf
    ub: float = math.inf
    guess: float = 0.0
    nominal: float = 1.0


@dataclass
class NLPConstraint:
    name: str
    expr: Expr
    lb: float
    ub: float
    nominal: float = 1.0


@dataclass
class MappedObjectiveBlock:
    name: str
    indices: list[int]
    local_size: int
    local_globals: list[list[int]]
    expr: Expr
    body: Callable[[LoopIndex], object] | None = None

    @property
    def count(self) -> int:
        return len(self.indices)


@dataclass
class MappedConstraintBlock:
    name: str
    indices: list[int]
    local_size: int
    local_globals: list[list[int]]
    exprs: list[object]
    lb: float
    ub: float
    nominal: float = 1.0
    body: Callable[[LoopIndex], object] | None = None

    @property
    def expr(self) -> Expr:
        return self.exprs[0]

    @property
    def output_size(self) -> int:
        total = 0
        for expr in self.exprs:
            total += len(expr) if isinstance(expr, (VectorExpr, BlockVectorExpr)) else 1
        return total

    @property
    def count(self) -> int:
        return len(self.indices)


@dataclass
class LocalBlockEmission:
    value: str
    jvp: str
    hvp: str
    jac: str
    hes: str
    jac_sparsity: list[tuple[int, int]]
    hes_sparsity: list[tuple[int, int]]
    jac_colors: list[int]
    hes_colors: list[int]
    jac_mode: str
    hes_mode: str
    report: dict[str, object]


class LoopIndex:
    def __init__(self, offset: int = 0, scale: int = 1, divisor: int = 1):
        self.offset = offset
        self.scale = scale
        self.divisor = divisor

    def eval(self, rep: int) -> int:
        return (self.scale * rep + self.offset) // self.divisor

    def __add__(self, other: int) -> "LoopIndex":
        # only valid when divisor == 1, i.e. no active //
        if self.divisor != 1:
            raise TypeError("cannot add after //")
        return LoopIndex(self.offset + int(other), self.scale)

    def __radd__(self, other: int) -> "LoopIndex":
        return self.__add__(other)

    def __sub__(self, other: int) -> "LoopIndex":
        if self.divisor != 1:
            raise TypeError("cannot subtract after //")
        return LoopIndex(self.offset - int(other), self.scale)

    def __mul__(self, other: int) -> "LoopIndex":
        if self.divisor != 1:
            raise TypeError("cannot multiply after //")
        factor = int(other)
        return LoopIndex(self.offset * factor, self.scale * factor)

    def __rmul__(self, other: int) -> "LoopIndex":
        return self.__mul__(other)

    def __floordiv__(self, other: int) -> "LoopIndex":
        if self.divisor != 1:
            raise TypeError("cannot chain //")
        if self.offset != 0:
            raise TypeError("// only supported for pure scale*rep, i.e. offset must be 0")
        return LoopIndex(0, self.scale, int(other))

class NLPVector:
    __moo_vector__ = True

    def __init__(self, model: "NLPModel", name: str, indices: Sequence[int]):
        self.model = model
        self.name = name
        self.indices = tuple(int(index) for index in indices)
        self.size = len(self.indices)
        self.start = self.indices[0] if self.indices else 0

    def __len__(self) -> int:
        return self.size

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return NLPVector(self.model, self.name, self.indices[idx])
        if isinstance(idx, LoopIndex):
            return self.model._loop_var(self, idx)
        idx = int(idx)
        if idx < 0:
            idx += self.size
        if idx < 0 or idx >= self.size:
            raise IndexError(idx)
        if self.model._loop_globals is not None:
            return self.model._loop_var(self, idx)
        global_idx = self.indices[idx]
        return Expr.variable("x", global_idx)

    def __sub__(self, other: object) -> "VectorExpr":
        if isinstance(other, (NLPVector, VectorExpr)):
            rhs = list(other)
            if len(rhs) != self.size:
                raise ValueError("vector sizes must match")
            return VectorExpr([self[i] - rhs[i] for i in range(self.size)])
        return VectorExpr([self[i] - other for i in range(self.size)])

    def __rsub__(self, other: object) -> "VectorExpr":
        if isinstance(other, (NLPVector, VectorExpr)):
            lhs = list(other)
            if len(lhs) != self.size:
                raise ValueError("vector sizes must match")
            return VectorExpr([lhs[i] - self[i] for i in range(self.size)])
        return VectorExpr([other - self[i] for i in range(self.size)])

    def __add__(self, other: object) -> "VectorExpr":
        if isinstance(other, (NLPVector, VectorExpr)):
            rhs = list(other)
            if len(rhs) != self.size:
                raise ValueError("vector sizes must match")
            return VectorExpr([self[i] + rhs[i] for i in range(self.size)])
        return VectorExpr([self[i] + other for i in range(self.size)])

    def __radd__(self, other: object) -> "VectorExpr":
        return self.__add__(other)

    def __mul__(self, other: object) -> "VectorExpr":
        return VectorExpr([self[i] * other for i in range(self.size)])

    def __rmul__(self, other: object) -> "VectorExpr":
        return self.__mul__(other)

    def __truediv__(self, other: object) -> "VectorExpr":
        return VectorExpr([self[i] / other for i in range(self.size)])

    def fix(self, idx, value: float) -> None:
        self.set_bounds(idx, lb=value, ub=value)
        self.set_guess(idx, value)

    def set_bounds(self, idx, lb: float = -math.inf, ub: float = math.inf) -> None:
        for global_idx in self._selected_indices(idx):
            self.model.variables[global_idx].lb = lb
            self.model.variables[global_idx].ub = ub

    def set_guess(self, idx, value: float) -> None:
        for global_idx in self._selected_indices(idx):
            self.model.variables[global_idx].guess = value

    def set_nominal(self, idx, value: float) -> None:
        for global_idx in self._selected_indices(idx):
            self.model.variables[global_idx].nominal = value

    def _selected_indices(self, idx) -> list[int]:
        if isinstance(idx, slice):
            return list(self.indices[idx])
        if isinstance(idx, Iterable) and not isinstance(idx, (str, bytes)):
            return [self.indices[int(i)] for i in idx]
        i = int(idx)
        if i < 0:
            i += self.size
        if i < 0 or i >= self.size:
            raise IndexError(idx)
        return [self.indices[i]]

class NLPModel(BaseModel):
    def __init__(self, name: str):
        super().__init__(name)
        self.variables: list[NLPVar] = []
        self.constraints: list[NLPConstraint] = []
        self.mapped_objectives: list[MappedObjectiveBlock] = []
        self.mapped_constraints: list[MappedConstraintBlock] = []
        self.objective: Expr | None = None
        self.obj_nominal = 1.0
        self._loop_indices: list[int] | None = None
        self._loop_globals: dict[tuple[int, ...], int] | None = None
        self._loop_global_lists: list[list[int]] | None = None
        self._loop_backend_values: list[object] | None = None
        self._loop_replay: bool = False

    @property
    def x_size(self) -> int:
        return len(self.variables)

    def add_variable(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> Expr:
        idx = len(self.variables)
        self.variables.append(NLPVar(name or f"x{idx}", lb, ub, guess, nominal))
        return Expr.variable("x", idx)

    def add_variables(self, name: str, size: int, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> NLPVector:
        start = len(self.variables)
        for i in range(size):
            self.variables.append(NLPVar(f"{name}{i}", lb, ub, guess, nominal))
        return NLPVector(self, name, range(start, start + size))

    def add_runtime_parameter(self, name: str, value: float) -> Expr:
        idx = len(self.runtime_parameters)
        self.runtime_parameters.append((name, value))
        return Expr.variable("rp", idx)

    def minimize(self, expr: object, nominal: float = 1.0, name: str | None = None) -> None:
        expr_value = as_expr(expr)
        self.objective = expr_value if self.objective is None else self.objective + expr_value
        self.obj_nominal = nominal

    def minimize_sum(self, count: int | Iterable[int], body: Callable[[LoopIndex], object], name: str | None = None) -> None:
        indices = self._normalize_loop_indices(count)
        idx = LoopIndex()
        self._start_loop_context(indices)
        expr = as_expr(body(idx))
        local_globals = self._finish_loop_context()
        self.mapped_objectives.append(MappedObjectiveBlock(name or f"obj_map{len(self.mapped_objectives)}", indices, len(local_globals), local_globals, expr, body))

    def add_constraint(self, expr: object, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        idx = len(self.constraints)
        self.constraints.append(NLPConstraint(name or f"g{idx}", as_expr(expr), low, high, nominal))

    def add_constraints(self, count_or_expr: int | Iterable[int] | VectorExpr, body: Callable[[LoopIndex], object] | None = None, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        if isinstance(count_or_expr, VectorExpr):
            for i, expr in enumerate(count_or_expr):
                self.add_constraint(expr, lb=low, ub=high, name=f"{name or 'g'}_{i}", nominal=nominal)
            return
        if body is None:
            raise ValueError("add_constraints(indices, body, ...) requires a body callable")
        self.add_constraints_map(count_or_expr, body, lb=low, ub=high, name=name, nominal=nominal)

    def add_constraints_map(self, count: int | Iterable[int], body: Callable[[LoopIndex], object], lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        indices = self._normalize_loop_indices(count)
        idx = LoopIndex()
        self._start_loop_context(indices)
        raw = body(idx)
        exprs = self._mapped_constraint_outputs(raw)
        local_globals = self._finish_loop_context()
        block_name = name or f"g_map{len(self.mapped_constraints)}"
        self.mapped_constraints.append(MappedConstraintBlock(block_name, indices, len(local_globals), local_globals, exprs, low, high, nominal, body))

    def _mapped_constraint_outputs(self, raw: object) -> list[object]:
        if isinstance(raw, VectorExpr):
            return [raw]
        if isinstance(raw, BlockVectorExpr):
            return [raw.flatten()]
        if isinstance(raw, (list, tuple)):
            raw_vec = vec(raw)
            if isinstance(raw_vec, VectorExpr) and any(isinstance(item, (VectorExpr, BlockVectorExpr, list, tuple)) for item in raw):
                return [raw_vec]
            return [as_expr(expr) for expr in raw]
        return [as_expr(raw)]

    def sumsqr(self, values: NLPVector | VectorExpr) -> Expr:
        from .expressions import sum_expr
        return sum_expr(value * value for value in values)

    def _loop_var(self, vector: NLPVector, idx: int | LoopIndex) -> Expr:
        if self._loop_globals is None:
            if isinstance(idx, LoopIndex):
                raise RuntimeError("loop index used outside a mapped block")
            globals_for_idx = vector._selected_indices(idx)
            global_idx = globals_for_idx[0]
            return Expr.variable("x", global_idx)
        if isinstance(idx, LoopIndex):
            globals_for_reps = [vector._selected_indices(idx.eval(iter_value))[0] for iter_value in self._loop_indices or []]
        else:
            globals_for_reps = [vector._selected_indices(idx)[0] for _ in self._loop_indices or []]
        key = tuple(globals_for_reps)
        local = self._loop_globals.get(key)
        if local is None:
            if self._loop_replay:
                raise RuntimeError("mapped block replay used a variable pattern that was not present during initial tracing")
            local = len(self._loop_globals)
            self._loop_globals[key] = local
            assert self._loop_global_lists is not None
            self._loop_global_lists.append(globals_for_reps)
        if self._loop_backend_values is not None:
            return Expr.variable("xl", local, graph_scalar=self._loop_backend_values[local])
        return Expr.variable("xl", local)

    def _normalize_loop_indices(self, indices: int | Iterable[int]) -> list[int]:
        if isinstance(indices, int):
            if indices < 0:
                raise ValueError("mapped block count must be non-negative")
            return list(range(indices))
        values = [int(index) for index in indices]
        return values

    def _start_loop_context(self, indices: list[int]) -> None:
        self._loop_indices = indices
        self._loop_globals = {}
        self._loop_global_lists = []
        self._loop_backend_values = None
        self._loop_replay = False

    def _start_loop_replay_context(self, indices: list[int], local_globals: list[list[int]], backend_values: list[object]) -> None:
        self._loop_indices = indices
        self._loop_globals = {tuple(values): local for local, values in enumerate(local_globals)}
        self._loop_global_lists = [list(values) for values in local_globals]
        self._loop_backend_values = backend_values
        self._loop_replay = True

    def _finish_loop_context(self) -> list[list[int]]:
        if self._loop_globals is None or self._loop_global_lists is None:
            raise RuntimeError("internal loop context is not active")
        local_globals = self._loop_global_lists
        self._clear_loop_context()
        return local_globals

    def _clear_loop_context(self) -> None:
        self._loop_indices = None
        self._loop_globals = None
        self._loop_global_lists = None
        self._loop_backend_values = None
        self._loop_replay = False

    def generate(self, out_dir: str | os.PathLike[str], repo_root: str | os.PathLike[str] | None = None, standalone_main: bool = True) -> tuple[Path, Path]:
        self._validate()
        out_path = Path(out_dir)
        if not out_path.is_absolute():
            out_path = out_path.resolve()
        out_path.mkdir(parents=True, exist_ok=True)
        emitted = self._emit_ad()
        c_path = out_path / f"{self.name}.c"
        h_path = out_path / f"{self.name}.h"
        c_path.write_text(self._render_c(emitted, standalone_main), encoding="utf-8")
        h_path.write_text(self._render_h(), encoding="utf-8")
        self._write_codegen_report(out_path, emitted, c_path)
        return c_path, h_path

    def compile(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, generate: bool = False) -> Path:
        out_path = Path(out_dir)
        abs_out_path = out_path if out_path.is_absolute() else out_path.resolve()
        c_path = abs_out_path / f"{self.name}.c"
        if generate or not c_path.exists():
            c_path, _ = self.generate(abs_out_path, repo_root, standalone_main=True)
        exe = abs_out_path / self.name
        cc = shutil.which("cc") or shutil.which("gcc")
        if cc is None:
            raise RuntimeError("No C compiler found in PATH")
        include = paths.include_dir()
        lib_dir = paths.library_dir(build_dir)
        subprocess.run([
            cc, "-std=c99", "-O3",
            f"-I{include}",
            f"-I{abs_out_path}",
            str(c_path),
            f"-L{lib_dir}",
            "-lmoo",
            f"-Wl,-rpath,{lib_dir}",
            "-lm",
            "-o", str(exe),
        ], check=True)
        return exe

    def optimize(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, solver: str | None = None, solver_args: list[str] | None = None, capture: bool = False, generate: bool = False, run_cwd: str | os.PathLike[str] | None = None) -> OptimizationRun:
        exe = self.compile(out_dir, build_dir=build_dir, repo_root=repo_root, generate=generate)
        cwd = Path(run_cwd) if run_cwd is not None else paths.library_dir(build_dir)
        if not cwd.is_absolute():
            cwd = cwd.resolve()
        cwd.mkdir(parents=True, exist_ok=True)
        process = subprocess.run([str(exe)] + self._solver_args(solver, solver_args), cwd=cwd, text=True, capture_output=capture, check=False)
        raw_results = read_results(cwd)
        result_view = NLPResult(
            raw_results,
            [var.name for var in self.variables],
            self._constraint_names(),
        )
        return OptimizationRun(process, raw_results, result_view)

    def run(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, solver: str | None = None, solver_args: list[str] | None = None, capture: bool = False, generate: bool = True, run_cwd: str | os.PathLike[str] | None = None) -> OptimizationRun:
        return self.optimize(out_dir, build_dir, repo_root, solver, solver_args, capture, generate, run_cwd)

    def _solver_args(self, solver: str | None, solver_args: list[str] | None) -> list[str]:
        return self.solver_settings.cli_args(solver, solver_args)

    def _validate(self) -> None:
        if not self.variables:
            raise ValueError("NLPModel requires at least one variable")
        if self.objective is None and not self.mapped_objectives:
            raise ValueError("NLPModel requires an objective; call minimize(...) or minimize_sum(...)")

    def _emit_ad(self) -> dict[str, str | list[tuple[int, int]]]:
        if self.mapped_objectives or self.mapped_constraints:
            return self._emit_structured_ad()
        outputs = [self.objective if self.objective is not None else 0.0]
        outputs.extend(c.expr for c in self.constraints)
        emitted = LocalGraphFunction(
            name="moo_nlp",
            input=InputGroup("x", self.x_size),
            outputs=outputs,
            params=[InputGroup("rp", len(self.runtime_parameters))],
        ).emit()
        jac_mode = select_derivative_callback_mode(self.derivative_strategy, emitted.jac_sparsity, emitted.jac_colors)
        hes_mode = select_derivative_callback_mode(self.derivative_strategy, emitted.hes_sparsity, emitted.hes_colors)
        return {
            "VALUE": emitted.value,
            "JVP": emitted.jvp,
            "HVP": emitted.hvp,
            "JAC": emitted.jac,
            "HES": emitted.hes,
            "JAC_SPARSITY": emitted.jac_sparsity,
            "HES_SPARSITY": emitted.hes_sparsity,
            "JAC_COLORS": emitted.jac_colors,
            "HES_COLORS": emitted.hes_colors,
            "JAC_MODE": jac_mode,
            "HES_MODE": hes_mode,
            "REPORT": emitted.report,
        }

    @property
    def g_size(self) -> int:
        return len(self.constraints) + sum(block.count * block.output_size for block in self.mapped_constraints)

    def _constraint_names(self) -> list[str]:
        names = [constraint.name for constraint in self.constraints]
        for block in self.mapped_constraints:
            if block.output_size == 1:
                names.extend(f"{block.name}_{i}" for i in range(block.count))
            else:
                names.extend(f"{block.name}_{i}_{j}" for i in range(block.count) for j in range(block.output_size))
        return names

    def _emit_local_block(self, prefix: str, block: MappedObjectiveBlock | MappedConstraintBlock) -> LocalBlockEmission:
        if block.body is not None:
            builder = ad.GraphFunctionBuilder()
            xl = builder.inputs("xl", block.local_size)
            rp = builder.params("rp", len(self.runtime_parameters))
            self._start_loop_replay_context(block.indices, block.local_globals, [xl[i] for i in range(block.local_size)])
            try:
                raw = block.body(LoopIndex())
                outputs = [as_expr(raw)] if isinstance(block, MappedObjectiveBlock) else self._mapped_constraint_outputs(raw)
                replayed_globals = self._finish_loop_context()
            except Exception:
                self._clear_loop_context()
                raise
            if replayed_globals != block.local_globals:
                raise RuntimeError("mapped block replay changed local variable ordering")
            emitted = emit_function_from_builder(
                builder,
                "xl",
                xl,
                outputs,
                {"rp": rp},
                f"{prefix}_value",
                f"{prefix}_jvp",
                f"{prefix}_hvp",
            )
        else:
            outputs = [block.expr] if isinstance(block, MappedObjectiveBlock) else list(block.exprs)
            emitted = LocalGraphFunction(
                name=prefix,
                input=InputGroup("xl", block.local_size),
                outputs=outputs,
                params=[InputGroup("rp", len(self.runtime_parameters))],
            ).emit()
        jac_mode = select_derivative_callback_mode(self.derivative_strategy, emitted.jac_sparsity, emitted.jac_colors)
        hes_mode = select_derivative_callback_mode(self.derivative_strategy, emitted.hes_sparsity, emitted.hes_colors)
        return LocalBlockEmission(
            value=emitted.value,
            jvp=emitted.jvp,
            hvp=emitted.hvp,
            jac=emitted.jac,
            hes=emitted.hes,
            jac_sparsity=emitted.jac_sparsity,
            hes_sparsity=emitted.hes_sparsity,
            jac_colors=emitted.jac_colors,
            hes_colors=emitted.hes_colors,
            jac_mode=jac_mode,
            hes_mode=hes_mode,
            report=emitted.report,
        )

    def _emit_structured_ad(self) -> dict[str, object]:
        scalar_emitted = None
        scalar_jac_mode = "direct"
        scalar_hes_mode = "direct"
        if self.objective is not None or self.constraints:
            scalar_outputs = [self.objective if self.objective is not None else 0.0]
            scalar_outputs.extend(c.expr for c in self.constraints)
            emitted = LocalGraphFunction(
                name="moo_nlp_scalar",
                input=InputGroup("x", self.x_size),
                outputs=scalar_outputs,
                params=[InputGroup("rp", len(self.runtime_parameters))],
            ).emit()
            scalar_jac_mode = select_derivative_callback_mode(self.derivative_strategy, emitted.jac_sparsity, emitted.jac_colors)
            scalar_hes_mode = select_derivative_callback_mode(self.derivative_strategy, emitted.hes_sparsity, emitted.hes_colors)
            scalar_emitted = emitted
        local_objectives = [
            (block, self._emit_local_block(f"moo_nlp_obj_map_{idx}", block))
            for idx, block in enumerate(self.mapped_objectives)
        ]
        local_constraints = [
            (block, self._emit_local_block(f"moo_nlp_g_map_{idx}", block))
            for idx, block in enumerate(self.mapped_constraints)
        ]

        obj_cols: dict[int, int] = {}
        g_jac: list[tuple[int, int]] = []
        g_buf: list[int] = []
        hes_pairs: dict[tuple[int, int], int] = {}

        def obj_buf_for(col: int) -> int:
            if col not in obj_cols:
                obj_cols[col] = len(obj_cols)
            return obj_cols[col]

        def hes_buf_for(row: int, col: int) -> int:
            if row < col:
                row, col = col, row
            key = (row, col)
            if key not in hes_pairs:
                hes_pairs[key] = len(hes_pairs)
            return hes_pairs[key]

        if scalar_emitted is not None:
            for row, col in scalar_emitted.jac_sparsity:
                if row == 0:
                    obj_buf_for(col)
            for row, col in scalar_emitted.hes_sparsity:
                hes_buf_for(row, col)

        for block, local in local_objectives:
            for rep in range(block.count):
                for row, col in local.jac_sparsity:
                    if row == 0:
                        obj_buf_for(block.local_globals[col][rep])
                for row, col in local.hes_sparsity:
                    hes_buf_for(block.local_globals[row][rep], block.local_globals[col][rep])

        g_buf_next = len(obj_cols)
        scalar_jac_buf: list[int] = []
        if scalar_emitted is not None:
            for row, col in scalar_emitted.jac_sparsity:
                if row == 0:
                    scalar_jac_buf.append(obj_buf_for(col))
                else:
                    g_jac.append((row - 1, col))
                    scalar_jac_buf.append(g_buf_next)
                    g_buf.append(g_buf_next)
                    g_buf_next += 1

        g_base = len(self.constraints)
        for block, local in local_constraints:
            for rep in range(block.count):
                for row, col in local.jac_sparsity:
                    global_row = g_base + rep * block.output_size + row
                    g_jac.append((global_row, block.local_globals[col][rep]))
                    g_buf.append(g_buf_next)
                    g_buf_next += 1
                for row, col in local.hes_sparsity:
                    hes_buf_for(block.local_globals[row][rep], block.local_globals[col][rep])
            g_base += block.count * block.output_size

        obj_jac = [(0, col) for col, _ in sorted(obj_cols.items(), key=lambda item: item[1])]
        hes = [pair for pair, _ in sorted(hes_pairs.items(), key=lambda item: item[1])]
        all_local = [local for _, local in local_objectives] + [local for _, local in local_constraints]
        jac_modes = {local.jac_mode for local in all_local}
        hes_modes = {local.hes_mode for local in all_local}
        report = {
            "graph_source": "structured_loop_blocks",
            "scalar_constraints": len(self.constraints),
            "scalar_objective": self.objective is not None,
            "mapped_objective_blocks": len(local_objectives),
            "mapped_constraint_blocks": len(local_constraints),
            "mapped_repetitions": sum(block.count for block, _ in local_objectives) + sum(block.count for block, _ in local_constraints),
            "mapped_constraint_outputs": sum(block.count * block.output_size for block, _ in local_constraints),
            "jacobian_nnz": len(obj_jac) + len(g_jac),
            "hessian_nnz": len(hes),
        }
        if scalar_emitted is not None:
            report["scalar_jacobian_nnz"] = len(scalar_emitted.jac_sparsity)
            report["scalar_hessian_nnz"] = len(scalar_emitted.hes_sparsity)
            report["scalar_jacobian_mode"] = scalar_jac_mode
            report["scalar_hessian_mode"] = scalar_hes_mode
        for idx, (block, local) in enumerate(local_objectives):
            prefix = f"local_objective_{idx}"
            report[f"{prefix}_repetitions"] = block.count
            report[f"{prefix}_input_size"] = block.local_size
            report[f"{prefix}_output_size"] = 1
            report[f"{prefix}_jacobian_nnz"] = len(local.jac_sparsity)
            report[f"{prefix}_hessian_nnz"] = len(local.hes_sparsity)
            report[f"{prefix}_jacobian_colors"] = local.report.get("jacobian_colors", 0)
            report[f"{prefix}_hessian_colors"] = local.report.get("hessian_colors", 0)
            report[f"{prefix}_jacobian_mode"] = local.jac_mode
            report[f"{prefix}_hessian_mode"] = local.hes_mode
        for idx, (block, local) in enumerate(local_constraints):
            prefix = f"local_constraint_{idx}"
            report[f"{prefix}_repetitions"] = block.count
            report[f"{prefix}_input_size"] = block.local_size
            report[f"{prefix}_output_size"] = block.output_size
            report[f"{prefix}_jacobian_nnz"] = len(local.jac_sparsity)
            report[f"{prefix}_hessian_nnz"] = len(local.hes_sparsity)
            report[f"{prefix}_jacobian_colors"] = local.report.get("jacobian_colors", 0)
            report[f"{prefix}_hessian_colors"] = local.report.get("hessian_colors", 0)
            report[f"{prefix}_jacobian_mode"] = local.jac_mode
            report[f"{prefix}_hessian_mode"] = local.hes_mode

        return {
            "STRUCTURED": True,
            "SCALAR": scalar_emitted,
            "SCALAR_JAC_MODE": scalar_jac_mode,
            "SCALAR_HES_MODE": scalar_hes_mode,
            "SCALAR_JAC_BUF": scalar_jac_buf,
            "LOCAL_OBJECTIVES": local_objectives,
            "LOCAL_CONSTRAINTS": local_constraints,
            "LOCAL_JAC_MODE": next(iter(jac_modes)) if len(jac_modes) == 1 else "mixed",
            "LOCAL_HES_MODE": next(iter(hes_modes)) if len(hes_modes) == 1 else "mixed",
            "OBJ_JAC": obj_jac,
            "OBJ_BUF": list(range(len(obj_jac))),
            "G_JAC": g_jac,
            "G_BUF": g_buf,
            "HES_SPARSITY": hes,
            "REPORT": report,
        }

    def _render_c(self, emitted: dict[str, str], standalone_main: bool) -> str:
        if emitted.get("STRUCTURED"):
            return self._render_structured_c(emitted, standalone_main)
        jac = parse_sparsity_pairs(emitted.get("JAC_SPARSITY", ""))
        hes = parse_sparsity_pairs(emitted.get("HES_SPARSITY", ""))
        obj_jac = [(0, col) for row, col in jac if row == 0]
        g_jac = [(row - 1, col) for row, col in jac if row >= 1]
        obj_buf = [idx for idx, (row, _) in enumerate(jac) if row == 0]
        g_buf = [idx for idx, (row, _) in enumerate(jac) if row >= 1]

        def coo(name: str, pairs: list[tuple[int, int]], bufs: list[int] | None = None) -> str:
            return f"""static coo_t {name} = {{
    .row = {_arr([r for r, _ in pairs])},
    .col = {_arr([c for _, c in pairs])},
    .buf_index = {_arr(list(range(len(pairs))) if bufs is None else bufs)},
    .nnz = {len(pairs)}
}};
"""

        def bounds(name: str, values) -> str:
            if not values:
                return f"static bounds_t {name}[1];\n"
            return f"static bounds_t {name}[{len(values)}] = {{ {', '.join(f'{{ {_num(v.lb)}, {_num(v.ub)} }}' for v in values)} }};\n"

        main = f"""
int main(int argc, char** argv) {{
    return main_{self.name}(argc, argv);
}}
""" if standalone_main else ""
        jac_mode = str(emitted.get("JAC_MODE", "direct"))
        hes_mode = str(emitted.get("HES_MODE", "direct"))
        sections = []
        sections.append("JAC" if jac_mode == "direct" else "JVP")
        sections.append("HES" if hes_mode == "direct" else "HVP")
        derivative_sections = "\n".join(
            emitted.get(key, "")
            for key in sections
        )
        jac_body = render_jacobian_callback_body(
            jac_mode,
            "moo_nlp_jvp_sparse(x, rp, out);",
            "X_SIZE",
            "OUT_SIZE",
            jac,
            emitted.get("JAC_COLORS", []) if isinstance(emitted.get("JAC_COLORS", []), list) else [],
            "moo_nlp_jvp(x, rp, v, tmp);",
        )
        hes_body = render_hessian_callback_body(
            hes_mode,
            """    f64 seed[OUT_SIZE] = {0};
    seed[0] = obj_factor;
    for (int i = 0; i < G_SIZE; ++i) { seed[1 + i] = lambda[i]; }
    moo_nlp_hvp_sparse(x, seed, rp, out);""",
            "X_SIZE",
            "X_SIZE",
            hes,
            emitted.get("HES_COLORS", []) if isinstance(emitted.get("HES_COLORS", []), list) else [],
            """    f64 seed[OUT_SIZE] = {0};
    moo_nlp_hvp_cache_t cache;
    seed[0] = obj_factor;
    for (int i = 0; i < G_SIZE; ++i) {{ seed[1 + i] = lambda[i]; }}
    moo_nlp_hvp_prepare(x, seed, rp, &cache);""",
            "moo_nlp_hvp_apply(&cache, v, tmp);",
        )

        return f"""#include <float.h>
#include <stdbool.h>
#include <math.h>

#include <interfaces/nlp/main_nlp.h>
#include "{self.name}.h"

#define X_SIZE {self.x_size}
#define RP_SIZE {len(self.runtime_parameters)}
#define G_SIZE {len(self.constraints)}
#define OUT_SIZE (1 + G_SIZE)

{emitted.get("VALUE", "")}
{derivative_sections}
static f64 globl_rp[RP_SIZE] = {{ {', '.join(_num(v) for _, v in self.runtime_parameters)} }};
static f64 globl_x0[X_SIZE] = {{ {', '.join(_num(v.guess) for v in self.variables)} }};
{bounds("globl_x_bounds", self.variables)}
{bounds("globl_g_bounds", self.constraints)}
static f64 globl_x_nominal[X_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.variables)} }};
static f64 globl_g_nominal[G_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.constraints)} }};
{coo("globl_obj_jac", obj_jac, obj_buf)}
{coo("globl_g_jac", g_jac, g_buf)}
{coo("globl_hes", hes)}

static void eval_all(const f64* x, const f64* rp, f64* out, void* user_data) {{
    moo_nlp_value(x, rp, out);
}}

static void jac_all(const f64* x, const f64* rp, f64* out, void* user_data) {{
{jac_body}
}}

static void hes_all(const f64* x, const f64* rp, const f64* lambda, f64 obj_factor, f64* out, void* user_data) {{
{hes_body}
}}

static c_nlp_callbacks_t globl_callbacks = {{ eval_all, jac_all, hes_all }};
static solver_ctx_t globl_solver_ctx = {{ .derivative_test = {_c_bool(self.derivative_test)} }};
static c_nlp_problem_t globl_problem = {{
    .x_size = X_SIZE,
    .rp_size = RP_SIZE,
    .g_size = G_SIZE,
    .rp = globl_rp,
    .x0 = globl_x0,
    .x_bounds = globl_x_bounds,
    .g_bounds = globl_g_bounds,
    .obj_nominal = {_num(self.obj_nominal)},
    .x_nominal = globl_x_nominal,
    .g_nominal = globl_g_nominal,
    .obj_jac = &globl_obj_jac,
    .g_jac = &globl_g_jac,
    .hes = &globl_hes,
    .callbacks = &globl_callbacks,
    .solver_ctx = &globl_solver_ctx,
    .user_data = (void*)0
}};

int main_{self.name}(int argc, char** argv) {{
    return main_nlp(argc, argv, &globl_problem);
}}
{main}
"""

    def _render_structured_c(self, emitted: dict[str, object], standalone_main: bool) -> str:
        scalar = emitted["SCALAR"]
        scalar_jac_buf = emitted["SCALAR_JAC_BUF"]
        obj_jac = emitted["OBJ_JAC"]
        obj_buf = emitted["OBJ_BUF"]
        g_jac = emitted["G_JAC"]
        g_buf = emitted["G_BUF"]
        hes = emitted["HES_SPARSITY"]
        local_objectives = emitted["LOCAL_OBJECTIVES"]
        local_constraints = emitted["LOCAL_CONSTRAINTS"]
        local_seed_size = max(
            [1 + len(self.constraints) if scalar is not None else 1]
            + [1 for _, _ in local_objectives]
            + [block.output_size for block, _ in local_constraints]
        )
        obj_buf_by_col = {col: buf for buf, (_, col) in zip(obj_buf, obj_jac)}
        hes_buf_by_pair = {pair: idx for idx, pair in enumerate(hes)}

        def global_at(block: MappedObjectiveBlock | MappedConstraintBlock, local_col: int, rep: int) -> int:
            return block.local_globals[local_col][rep]

        def hes_buf_table(block: MappedObjectiveBlock | MappedConstraintBlock, local: LocalBlockEmission) -> list[int]:
            table = []
            for rep in range(block.count):
                for row, col in local.hes_sparsity:
                    gr = global_at(block, row, rep)
                    gc = global_at(block, col, rep)
                    if gr < gc:
                        gr, gc = gc, gr
                    table.append(hes_buf_by_pair[(gr, gc)])
            return table

        def jac_buf_table(block: MappedObjectiveBlock | MappedConstraintBlock, local: LocalBlockEmission, objective: bool, base_buf: int = 0) -> list[int]:
            table = []
            for rep in range(block.count):
                for local_buf, (row, col) in enumerate(local.jac_sparsity):
                    if objective:
                        table.append(obj_buf_by_col[global_at(block, col, rep)] if row == 0 else -1)
                    else:
                        table.append(base_buf + rep * len(local.jac_sparsity) + local_buf)
            return table

        def affine(values: list[int]) -> tuple[int, int] | None:
            if not values:
                return (0, 0)
            if len(values) == 1:
                return (values[0], 0)
            step = values[1] - values[0]
            if all(values[i] == values[0] + i * step for i in range(len(values))):
                return (values[0], step)
            return None

        def gather_lines(block: MappedObjectiveBlock | MappedConstraintBlock, indent: str, prefix: str) -> list[str]:
            lines = []
            for local_idx, globals_for_reps in enumerate(block.local_globals):
                regular = affine(globals_for_reps)
                if regular is not None:
                    base, step = regular
                    if step == 0:
                        lines.append(f"{indent}xl[{local_idx}] = x[{base}];")
                    else:
                        lines.append(f"{indent}xl[{local_idx}] = x[{base} + rep * {step}];")
                else:
                    name = f"{prefix}_gidx_{local_idx}"
                    lines.append(static_int_array(name, globals_for_reps, indent))
                    lines.append(f"{indent}xl[{local_idx}] = x[{name}[rep]];")
            return lines

        def coo(name: str, pairs: list[tuple[int, int]], bufs: list[int] | None = None) -> str:
            return f"""static coo_t {name} = {{
    .row = {_arr([r for r, _ in pairs])},
    .col = {_arr([c for _, c in pairs])},
    .buf_index = {_arr(list(range(len(pairs))) if bufs is None else bufs)},
    .nnz = {len(pairs)}
}};
"""

        def bounds(name: str, values) -> str:
            if not values:
                return f"static bounds_t {name}[1];\n"
            return f"static bounds_t {name}[{len(values)}] = {{ {', '.join(f'{{ {_num(v.lb)}, {_num(v.ub)} }}' for v in values)} }};\n"

        g_bounds_values: list[NLPConstraint] = []
        g_bounds_values.extend(self.constraints)
        for block, _ in local_constraints:
            g_bounds_values.extend(
                NLPConstraint(f"{block.name}_{i}_{j}", Expr.const(0.0), block.lb, block.ub, block.nominal)
                for i in range(block.count)
                for j in range(block.output_size)
            )

        local_code = []
        if scalar is not None:
            local_code.append(scalar.value)
            local_code.append(scalar.jac)
            local_code.append(scalar.hes)
        for _, local in local_objectives:
            local_code.append(local.value)
            local_code.append(local.jac if local.jac_mode == "direct" else local.jvp)
            local_code.append(local.hes if local.hes_mode == "direct" else local.hvp)
        for _, local in local_constraints:
            local_code.append(local.value)
            local_code.append(local.jac if local.jac_mode == "direct" else local.jvp)
            local_code.append(local.hes if local.hes_mode == "direct" else local.hvp)

        eval_lines = ["    out[0] = 0.0;"]
        if scalar is not None:
            eval_lines.append(f"    f64 scalar_tmp[{max(1 + len(self.constraints), 1)}];")
            eval_lines.append("    moo_nlp_scalar_value(x, rp, scalar_tmp);")
            eval_lines.append("    out[0] += scalar_tmp[0];")
            for i in range(len(self.constraints)):
                eval_lines.append(f"    out[1 + {i}] = scalar_tmp[1 + {i}];")
        for idx, (block, local) in enumerate(local_objectives):
            eval_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            eval_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            eval_lines.extend(gather_lines(block, "        ", f"obj_{idx}"))
            eval_lines.append("        f64 tmp[1];")
            eval_lines.append(f"        moo_nlp_obj_map_{idx}_value(xl, rp, tmp);")
            eval_lines.append("        out[0] += tmp[0];")
            eval_lines.append("    }")

        g_base = len(self.constraints)
        for idx, (block, local) in enumerate(local_constraints):
            eval_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            eval_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            eval_lines.extend(gather_lines(block, "        ", f"g_{idx}"))
            eval_lines.append(f"        f64 tmp[{max(block.output_size, 1)}];")
            eval_lines.append(f"        moo_nlp_g_map_{idx}_value(xl, rp, tmp);")
            for local_row in range(block.output_size):
                eval_lines.append(f"        out[1 + {g_base} + rep * {block.output_size} + {local_row}] = tmp[{local_row}];")
            eval_lines.append("    }")
            g_base += block.count * block.output_size

        jac_lines = [f"    for (int i = 0; i < JAC_BUF_SIZE; ++i) {{ out[i] = 0.0; }}"]
        if scalar is not None and scalar.jac_sparsity:
            jac_lines.append(f"    f64 scalar_jac[{max(len(scalar.jac_sparsity), 1)}];")
            jac_lines.append("    moo_nlp_scalar_jvp_sparse(x, rp, scalar_jac);")
            jac_lines.append(static_int_array("scalar_jac_buf", scalar_jac_buf, "    "))
            for local_buf, _ in enumerate(scalar.jac_sparsity):
                jac_lines.append(f"    out[scalar_jac_buf[{local_buf}]] += scalar_jac[{local_buf}];")
        for idx, (block, local) in enumerate(local_objectives):
            if not local.jac_sparsity:
                continue
            jac_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            jac_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            jac_lines.extend(gather_lines(block, "        ", f"obj_jac_{idx}"))
            obj_jbuf = jac_buf_table(block, local, objective=True)
            jac_lines.append(static_int_array(f"obj_{idx}_jbuf", obj_jbuf, "        "))
            if local.jac_mode == "direct":
                jac_lines.extend(render_local_direct_objective_jac_lines(f"moo_nlp_obj_map_{idx}_jvp", local.jac_sparsity, f"obj_{idx}_jbuf", "        "))
            else:
                jac_lines.extend(
                    render_local_colored_jac_lines(
                        f"moo_nlp_obj_map_{idx}_jvp",
                        local.jac_sparsity,
                        local.jac_colors,
                        block.local_size,
                        1,
                        "out",
                        "rep",
                        accumulate=True,
                        indent="        ",
                        jbuf_name=f"obj_{idx}_jbuf",
                    )
                )
            jac_lines.append("    }")

        g_buf_cursor = len(obj_jac)
        if scalar is not None:
            g_buf_cursor += sum(1 for row, _ in scalar.jac_sparsity if row >= 1)
        for idx, (block, local) in enumerate(local_constraints):
            if local.jac_sparsity:
                jac_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
                jac_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
                jac_lines.extend(gather_lines(block, "        ", f"g_jac_{idx}"))
                if local.jac_mode == "direct":
                    jac_lines.extend(render_local_direct_constraint_jac_lines(f"moo_nlp_g_map_{idx}_jvp", local.jac_sparsity, g_buf_cursor, "        "))
                else:
                    jac_lines.extend(
                        render_local_colored_jac_lines(
                            f"moo_nlp_g_map_{idx}_jvp",
                            local.jac_sparsity,
                            local.jac_colors,
                            block.local_size,
                            block.output_size,
                            "out",
                            "rep",
                            accumulate=False,
                            indent="        ",
                            base_buf=g_buf_cursor,
                        )
                    )
                jac_lines.append("    }")
                g_buf_cursor += block.count * len(local.jac_sparsity)

        hes_lines = [f"    for (int i = 0; i < HES_SIZE; ++i) {{ out[i] = 0.0; }}", f"    f64 seed[{local_seed_size}] = {{0}};"]
        if scalar is not None and scalar.hes_sparsity:
            scalar_hbuf = []
            for row, col in scalar.hes_sparsity:
                r, c = (row, col) if row >= col else (col, row)
                scalar_hbuf.append(hes_buf_by_pair[(r, c)])
            hes_lines.append(f"    f64 scalar_seed[{max(1 + len(self.constraints), 1)}] = {{0}};")
            hes_lines.append("    scalar_seed[0] = obj_factor;")
            for i in range(len(self.constraints)):
                hes_lines.append(f"    scalar_seed[1 + {i}] = lambda[{i}];")
            hes_lines.append(f"    f64 scalar_hes[{max(len(scalar.hes_sparsity), 1)}];")
            hes_lines.append("    moo_nlp_scalar_hvp_sparse(x, scalar_seed, rp, scalar_hes);")
            hes_lines.append(static_int_array("scalar_hbuf", scalar_hbuf, "    "))
            for local_buf, _ in enumerate(scalar.hes_sparsity):
                hes_lines.append(f"    out[scalar_hbuf[{local_buf}]] += scalar_hes[{local_buf}];")
        for idx, (block, local) in enumerate(local_objectives):
            if not local.hes_sparsity:
                continue
            hes_lines.append("    seed[0] = obj_factor;")
            hes_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            hes_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            hes_lines.extend(gather_lines(block, "        ", f"obj_hes_{idx}"))
            if local.hes_mode == "direct":
                hbuf = hes_buf_table(block, local)
                hes_lines.append(static_int_array(f"obj_{idx}_hbuf", hbuf, "        "))
                hes_lines.extend(render_local_direct_hes_lines(f"moo_nlp_obj_map_{idx}_hvp", local.hes_sparsity, f"obj_{idx}_hbuf", "tmp", "        "))
            else:
                hes_lines.extend(render_local_colored_hes_lines(f"moo_nlp_obj_map_{idx}_hvp", local.hes_sparsity, local.hes_colors, block.local_size, hes_buf_table(block, local), "        "))
            hes_lines.append("    }")

        g_base = len(self.constraints)
        for idx, (block, local) in enumerate(local_constraints):
            if local.hes_sparsity:
                hes_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
                hes_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
                hes_lines.extend(gather_lines(block, "        ", f"g_hes_{idx}"))
                hbuf = hes_buf_table(block, local)
                hes_lines.append(static_int_array(f"g_{idx}_hbuf", hbuf, "        "))
                for local_row in range(block.output_size):
                    hes_lines.append(f"        for (int seed_i = 0; seed_i < {block.output_size}; ++seed_i) {{ seed[seed_i] = 0.0; }}")
                    hes_lines.append(f"        seed[{local_row}] = lambda[{g_base} + rep * {block.output_size} + {local_row}];")
                    if local.hes_mode == "direct":
                        hes_lines.extend(render_local_direct_hes_lines(f"moo_nlp_g_map_{idx}_hvp", local.hes_sparsity, f"g_{idx}_hbuf", f"tmp_{local_row}", "        "))
                    else:
                        hes_lines.extend(render_local_colored_hes_lines(f"moo_nlp_g_map_{idx}_hvp", local.hes_sparsity, local.hes_colors, block.local_size, hbuf, "        ", hbuf_name=f"g_{idx}_hbuf"))
                hes_lines.append("    }")
            g_base += block.count * block.output_size

        main = f"""
int main(int argc, char** argv) {{
    return main_{self.name}(argc, argv);
}}
""" if standalone_main else ""

        return f"""#include <float.h>
#include <stdbool.h>
#include <math.h>

#include <interfaces/nlp/main_nlp.h>
#include "{self.name}.h"

#define X_SIZE {self.x_size}
#define RP_SIZE {len(self.runtime_parameters)}
#define G_SIZE {self.g_size}
#define OUT_SIZE (1 + G_SIZE)
#define JAC_BUF_SIZE {len(obj_jac) + len(g_jac)}
#define HES_SIZE {len(hes)}

{''.join(local_code)}
static f64 globl_rp[RP_SIZE] = {{ {', '.join(_num(v) for _, v in self.runtime_parameters)} }};
static f64 globl_x0[X_SIZE] = {{ {', '.join(_num(v.guess) for v in self.variables)} }};
{bounds("globl_x_bounds", self.variables)}
{bounds("globl_g_bounds", g_bounds_values)}
static f64 globl_x_nominal[X_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.variables)} }};
static f64 globl_g_nominal[G_SIZE] = {{ {', '.join(_num(v.nominal) for v in g_bounds_values)} }};
{coo("globl_obj_jac", obj_jac, obj_buf)}
{coo("globl_g_jac", g_jac, g_buf)}
{coo("globl_hes", hes)}

static void eval_all(const f64* x, const f64* rp, f64* out, void* user_data) {{
{chr(10).join(eval_lines)}
}}

static void jac_all(const f64* x, const f64* rp, f64* out, void* user_data) {{
{chr(10).join(jac_lines)}
}}

static void hes_all(const f64* x, const f64* rp, const f64* lambda, f64 obj_factor, f64* out, void* user_data) {{
{chr(10).join(hes_lines)}
}}

static c_nlp_callbacks_t globl_callbacks = {{ eval_all, jac_all, hes_all }};
static solver_ctx_t globl_solver_ctx = {{ .derivative_test = {_c_bool(self.derivative_test)} }};
static c_nlp_problem_t globl_problem = {{
    .x_size = X_SIZE,
    .rp_size = RP_SIZE,
    .g_size = G_SIZE,
    .rp = globl_rp,
    .x0 = globl_x0,
    .x_bounds = globl_x_bounds,
    .g_bounds = globl_g_bounds,
    .obj_nominal = {_num(self.obj_nominal)},
    .x_nominal = globl_x_nominal,
    .g_nominal = globl_g_nominal,
    .obj_jac = &globl_obj_jac,
    .g_jac = &globl_g_jac,
    .hes = &globl_hes,
    .callbacks = &globl_callbacks,
    .solver_ctx = &globl_solver_ctx,
    .user_data = (void*)0
}};

int main_{self.name}(int argc, char** argv) {{
    return main_nlp(argc, argv, &globl_problem);
}}
{main}
"""

    def _render_h(self) -> str:
        guard = f"MOO_NLP_CODEGEN_{self.name.upper()}_H"
        return f"""#ifndef {guard}
#define {guard}

#ifdef __cplusplus
extern "C" {{
#endif

int main_{self.name}(int argc, char** argv);

#ifdef __cplusplus
}}
#endif

#endif
"""

    def _write_codegen_report(self, out_path: Path, emitted: dict[str, object], c_path: Path) -> None:
        report = emitted.get("REPORT", {})
        lines = [
            f"model={self.name}",
            "problem=NLP",
            f"derivative_strategy={self.derivative_strategy}",
            f"generated_c_bytes={c_path.stat().st_size}",
        ]
        if emitted.get("STRUCTURED"):
            jac_mode = "loop"
            hes_mode = "loop"
            lines.extend([
                "structure_mode=loop",
                f"local_jacobian_mode={emitted.get('LOCAL_JAC_MODE', 'direct')}",
                f"local_hessian_mode={emitted.get('LOCAL_HES_MODE', 'direct')}",
            ])
        else:
            jac_mode = str(emitted.get("JAC_MODE", "direct"))
            hes_mode = str(emitted.get("HES_MODE", "direct"))
            lines.append("structure_mode=scalar")
        lines.extend([f"jacobian_mode={jac_mode}", f"hessian_mode={hes_mode}"])
        if isinstance(report, dict):
            lines.extend(f"derivative_{key}={value}" for key, value in sorted(report.items()))
        (out_path / "codegen_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
