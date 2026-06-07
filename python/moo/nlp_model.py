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
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .ad_codegen import emit_function
from .common import SolverMixin, SolverSettings, derivative_mode
from .model import Expr, as_expr, _arr, _c_bool, _clean_name, _num
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
    count: int
    local_size: int
    offsets: list[int]
    expr: Expr


@dataclass
class MappedConstraintBlock:
    name: str
    count: int
    local_size: int
    offsets: list[int]
    exprs: list[Expr]
    lb: float
    ub: float
    nominal: float = 1.0

    @property
    def expr(self) -> Expr:
        return self.exprs[0]

    @property
    def output_size(self) -> int:
        return len(self.exprs)


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
    def __init__(self, offset: int = 0):
        self.offset = offset

    def __add__(self, other: int) -> "LoopIndex":
        return LoopIndex(self.offset + int(other))

    def __radd__(self, other: int) -> "LoopIndex":
        return self.__add__(other)

    def __sub__(self, other: int) -> "LoopIndex":
        return LoopIndex(self.offset - int(other))


class NLPVector:
    def __init__(self, model: "NLPModel", name: str, start: int, size: int):
        self.model = model
        self.name = name
        self.start = start
        self.size = size

    def __len__(self) -> int:
        return self.size

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            start, stop, step = idx.indices(self.size)
            if step != 1:
                raise ValueError("NLPVector slices only support step=1")
            return NLPVector(self.model, self.name, self.start + start, stop - start)
        if isinstance(idx, LoopIndex):
            return self.model._loop_var(self, idx.offset)
        idx = int(idx)
        if idx < 0:
            idx += self.size
        if idx < 0 or idx >= self.size:
            raise IndexError(idx)
        return Expr(f"x[{self.start + idx}]", f"x[{self.start + idx}]")

    def __sub__(self, other: object) -> "VectorExpr":
        return VectorExpr([self[i] - other for i in range(self.size)])

    def __add__(self, other: object) -> "VectorExpr":
        if isinstance(other, NLPVector):
            if len(self) != len(other):
                raise ValueError("vector sizes must match")
            return VectorExpr([self[i] + other[i] for i in range(self.size)])
        return VectorExpr([self[i] + other for i in range(self.size)])


class VectorExpr:
    def __init__(self, values: list[Expr]):
        self.values = values

    def __len__(self) -> int:
        return len(self.values)

    def __getitem__(self, idx):
        if isinstance(idx, slice):
            return VectorExpr(self.values[idx])
        return self.values[idx]

    def __iter__(self):
        return iter(self.values)


class NLPModel(SolverMixin):
    def __init__(self, name: str):
        self.name = _clean_name(name)
        self.solver_settings = SolverSettings()
        self.codegen_strategy = "auto"
        self.codegen_structure_strategy = "auto"
        self.codegen_local_strategy = "auto"
        self.variables: list[NLPVar] = []
        self.constraints: list[NLPConstraint] = []
        self.mapped_objectives: list[MappedObjectiveBlock] = []
        self.mapped_constraints: list[MappedConstraintBlock] = []
        self.runtime_parameters: list[tuple[str, float]] = []
        self.objective: Expr | None = None
        self.obj_nominal = 1.0
        self.derivative_test = False
        self._loop_offsets: dict[int, int] | None = None

    @property
    def x_size(self) -> int:
        return len(self.variables)

    def add_variable(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> Expr:
        idx = len(self.variables)
        self.variables.append(NLPVar(name or f"x{idx}", lb, ub, guess, nominal))
        return Expr(f"x[{idx}]", f"x[{idx}]")

    def add_variables(self, name: str, size: int, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> NLPVector:
        start = len(self.variables)
        for i in range(size):
            self.variables.append(NLPVar(f"{name}{i}", lb, ub, guess, nominal))
        return NLPVector(self, name, start, size)

    def add_runtime_parameter(self, name: str, value: float) -> Expr:
        idx = len(self.runtime_parameters)
        self.runtime_parameters.append((name, value))
        return Expr(f"rp[{idx}]", f"rp[{idx}]")

    def minimize(self, expr: object, nominal: float = 1.0) -> None:
        self.objective = as_expr(expr)
        self.obj_nominal = nominal

    def minimize_sum(self, count: int, body: Callable[[LoopIndex], object], name: str | None = None) -> None:
        idx = LoopIndex()
        self._loop_offsets = {}
        expr = as_expr(body(idx))
        offsets = self._finish_loop_offsets(count)
        self.mapped_objectives.append(MappedObjectiveBlock(name or f"obj_map{len(self.mapped_objectives)}", count, len(offsets), offsets, expr))

    def add_constraint(self, expr: object, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        idx = len(self.constraints)
        self.constraints.append(NLPConstraint(name or f"g{idx}", as_expr(expr), low, high, nominal))

    def add_constraints(self, count_or_expr: int | VectorExpr, body: Callable[[LoopIndex], object] | None = None, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        if isinstance(count_or_expr, VectorExpr):
            for i, expr in enumerate(count_or_expr):
                self.add_constraint(expr, lb=low, ub=high, name=f"{name or 'g'}_{i}", nominal=nominal)
            return
        if body is None:
            raise ValueError("add_constraints(count, body, ...) requires a body callable")
        self.add_constraints_map(int(count_or_expr), body, lb=low, ub=high, name=name, nominal=nominal)

    def add_constraints_map(self, count: int, body: Callable[[LoopIndex], object], lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        idx = LoopIndex()
        self._loop_offsets = {}
        raw = body(idx)
        if isinstance(raw, VectorExpr):
            exprs = [as_expr(expr) for expr in raw]
        elif isinstance(raw, (list, tuple)):
            exprs = [as_expr(expr) for expr in raw]
        else:
            exprs = [as_expr(raw)]
        offsets = self._finish_loop_offsets(count)
        block_name = name or f"g_map{len(self.mapped_constraints)}"
        self.mapped_constraints.append(MappedConstraintBlock(block_name, count, len(offsets), offsets, exprs, low, high, nominal))

    def sumsqr(self, values: NLPVector | VectorExpr) -> Expr:
        from .model import sum_expr
        return sum_expr(value * value for value in values)

    def _loop_var(self, vector: NLPVector, offset: int) -> Expr:
        if self._loop_offsets is None:
            if offset < 0 or offset >= vector.size:
                raise IndexError(offset)
            return vector[offset]
        local = self._loop_offsets.setdefault(vector.start + offset, len(self._loop_offsets))
        return Expr(f"xl[{local}]", f"xl[{local}]")

    def _finish_loop_offsets(self, count: int) -> list[int]:
        if self._loop_offsets is None:
            raise RuntimeError("internal loop context is not active")
        offsets = [global_offset for global_offset, _ in sorted(self._loop_offsets.items(), key=lambda item: item[1])]
        self._loop_offsets = None
        for offset in offsets:
            if offset < 0 or offset + count > self.x_size:
                raise ValueError(f"mapped block index offset {offset} is outside the variable range for count={count}")
        return offsets

    def generate(self, out_dir: str | os.PathLike[str], repo_root: str | os.PathLike[str] | None = None, standalone_main: bool = True) -> tuple[Path, Path]:
        self._validate()
        root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        emitted = self._emit_ad()
        c_path = out_path / f"{self.name}.c"
        h_path = out_path / f"{self.name}.h"
        c_path.write_text(self._render_c(emitted, standalone_main), encoding="utf-8")
        h_path.write_text(self._render_h(), encoding="utf-8")
        self._write_codegen_report(out_path, emitted, c_path)
        return c_path, h_path

    def compile(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, generate: bool = False) -> Path:
        root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
        out_path = Path(out_dir)
        abs_out_path = out_path if out_path.is_absolute() else root / out_path
        c_path = out_path / f"{self.name}.c"
        if generate or not c_path.exists():
            c_path, _ = self.generate(out_path, root, standalone_main=True)
        exe = abs_out_path / self.name
        cc = shutil.which("cc") or shutil.which("gcc")
        if cc is None:
            raise RuntimeError("No C compiler found in PATH")
        build = Path(build_dir)
        subprocess.run([
            cc, "-std=c99", "-O3",
            f"-I{root / 'src'}",
            f"-I{abs_out_path}",
            str(root / c_path if not c_path.is_absolute() else c_path),
            f"-L{root / build}",
            "-lmoo",
            f"-Wl,-rpath,{root / build}",
            "-lm",
            "-o", str(exe),
        ], cwd=root, check=True)
        return exe

    def optimize(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, solver: str | None = None, solver_args: list[str] | None = None, capture: bool = False, generate: bool = False, run_cwd: str | os.PathLike[str] | None = None) -> OptimizationRun:
        root = Path(repo_root) if repo_root is not None else Path(__file__).resolve().parents[2]
        exe = self.compile(out_dir, build_dir=build_dir, repo_root=root, generate=generate)
        cwd = Path(run_cwd) if run_cwd is not None else root / Path(build_dir)
        if not cwd.is_absolute():
            cwd = root / cwd
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
        outputs = [self.objective.lfg if self.objective is not None else "0.0"]
        outputs.extend(c.expr.lfg for c in self.constraints)
        emitted = emit_function(
            "x",
            self.x_size,
            outputs,
            [("rp", len(self.runtime_parameters))],
            "moo_nlp_value",
            "moo_nlp_jvp",
            "moo_nlp_hvp",
        )
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
        outputs = [block.expr.lfg] if isinstance(block, MappedObjectiveBlock) else [expr.lfg for expr in block.exprs]
        emitted = emit_function(
            "xl",
            block.local_size,
            outputs,
            [("rp", len(self.runtime_parameters))],
            f"{prefix}_value",
            f"{prefix}_jvp",
            f"{prefix}_hvp",
        )
        jac_mode = derivative_mode(self.codegen_local_strategy, emitted.jac_sparsity, emitted.jac_colors)
        hes_mode = derivative_mode(self.codegen_local_strategy, emitted.hes_sparsity, emitted.hes_colors)
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
        if self.objective is not None or self.constraints:
            raise NotImplementedError("structured NLP codegen currently supports mapped objectives/constraints without additional scalar terms")
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

        for block, local in local_objectives:
            for rep in range(block.count):
                for row, col in local.jac_sparsity:
                    if row == 0:
                        obj_buf_for(rep + block.offsets[col])
                for row, col in local.hes_sparsity:
                    hes_buf_for(rep + block.offsets[row], rep + block.offsets[col])

        g_base = len(self.constraints)
        g_buf_next = len(obj_cols)
        for block, local in local_constraints:
            for rep in range(block.count):
                for row, col in local.jac_sparsity:
                    global_row = g_base + rep * block.output_size + row
                    g_jac.append((global_row, rep + block.offsets[col]))
                    g_buf.append(g_buf_next)
                    g_buf_next += 1
                for row, col in local.hes_sparsity:
                    hes_buf_for(rep + block.offsets[row], rep + block.offsets[col])
            g_base += block.count * block.output_size

        obj_jac = [(0, col) for col, _ in sorted(obj_cols.items(), key=lambda item: item[1])]
        hes = [pair for pair, _ in sorted(hes_pairs.items(), key=lambda item: item[1])]
        all_local = [local for _, local in local_objectives] + [local for _, local in local_constraints]
        jac_modes = {local.jac_mode for local in all_local}
        hes_modes = {local.hes_mode for local in all_local}
        report = {
            "kernel_source": "structured_loop_blocks",
            "mapped_objective_blocks": len(local_objectives),
            "mapped_constraint_blocks": len(local_constraints),
            "mapped_repetitions": sum(block.count for block, _ in local_objectives) + sum(block.count for block, _ in local_constraints),
            "mapped_constraint_outputs": sum(block.count * block.output_size for block, _ in local_constraints),
            "jacobian_nnz": len(obj_jac) + len(g_jac),
            "hessian_nnz": len(hes),
        }
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

    def _pairs(self, text: str) -> list[tuple[int, int]]:
        if isinstance(text, list):
            return text
        pairs = []
        for line in text.splitlines():
            if line.strip():
                a, b = line.split(",", 1)
                pairs.append((int(a), int(b)))
        return pairs

    def _render_c(self, emitted: dict[str, str], standalone_main: bool) -> str:
        if emitted.get("STRUCTURED"):
            return self._render_structured_c(emitted, standalone_main)
        jac = self._pairs(emitted.get("JAC_SPARSITY", ""))
        hes = self._pairs(emitted.get("HES_SPARSITY", ""))
        jac_colors = emitted.get("JAC_COLORS", [])
        hes_colors = emitted.get("HES_COLORS", [])
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
        jac_mode = derivative_mode(self.codegen_local_strategy, jac, jac_colors if isinstance(jac_colors, list) else [])
        hes_mode = derivative_mode(self.codegen_local_strategy, hes, hes_colors if isinstance(hes_colors, list) else [])
        sections = []
        sections.append("JAC" if jac_mode == "direct" else "JVP")
        sections.append("HES" if hes_mode == "direct" else "HVP")
        derivative_sections = "\n".join(
            emitted.get(key, "")
            for key in sections
        )
        if jac_mode == "direct":
            jac_body = "    moo_nlp_jvp_sparse(x, rp, out);"
        elif jac_mode == "colored":
            jac_body = f"""    f64 v[X_SIZE] = {{0}};
    f64 tmp[OUT_SIZE] = {{0}};
{self._render_colored_jac_fill(jac, jac_colors if isinstance(jac_colors, list) else [])}"""
        else:
            jac_body = f"""    f64 v[X_SIZE] = {{0}};
    f64 tmp[OUT_SIZE] = {{0}};
{self._render_jac_fill(jac)}"""
        if hes_mode == "direct":
            hes_body = """    f64 seed[OUT_SIZE] = {0};
    seed[0] = obj_factor;
    for (int i = 0; i < G_SIZE; ++i) { seed[1 + i] = lambda[i]; }
    moo_nlp_hvp_sparse(x, seed, rp, out);"""
        elif hes_mode == "colored":
            hes_body = f"""    f64 seed[OUT_SIZE] = {{0}};
    f64 v[X_SIZE] = {{0}};
    f64 tmp[X_SIZE] = {{0}};
    moo_nlp_hvp_cache_t cache;
    seed[0] = obj_factor;
    for (int i = 0; i < G_SIZE; ++i) {{ seed[1 + i] = lambda[i]; }}
    moo_nlp_hvp_prepare(x, seed, rp, &cache);
{self._render_colored_hes_fill(hes, hes_colors if isinstance(hes_colors, list) else [])}"""
        else:
            hes_body = f"""    f64 seed[OUT_SIZE] = {{0}};
    f64 v[X_SIZE] = {{0}};
    f64 tmp[X_SIZE] = {{0}};
    moo_nlp_hvp_cache_t cache;
    seed[0] = obj_factor;
    for (int i = 0; i < G_SIZE; ++i) {{ seed[1 + i] = lambda[i]; }}
    moo_nlp_hvp_prepare(x, seed, rp, &cache);
{self._render_hes_fill(hes)}"""

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
        obj_jac = emitted["OBJ_JAC"]
        obj_buf = emitted["OBJ_BUF"]
        g_jac = emitted["G_JAC"]
        g_buf = emitted["G_BUF"]
        hes = emitted["HES_SPARSITY"]
        local_objectives = emitted["LOCAL_OBJECTIVES"]
        local_constraints = emitted["LOCAL_CONSTRAINTS"]
        local_seed_size = max(
            [1]
            + [1 for _, _ in local_objectives]
            + [block.output_size for block, _ in local_constraints]
        )
        obj_buf_by_col = {col: buf for buf, (_, col) in zip(obj_buf, obj_jac)}
        hes_buf_by_pair = {pair: idx for idx, pair in enumerate(hes)}

        def hes_buf_table(block: MappedObjectiveBlock | MappedConstraintBlock, local: LocalBlockEmission) -> list[int]:
            table = []
            for rep in range(block.count):
                for row, col in local.hes_sparsity:
                    gr = rep + block.offsets[row]
                    gc = rep + block.offsets[col]
                    if gr < gc:
                        gr, gc = gc, gr
                    table.append(hes_buf_by_pair[(gr, gc)])
            return table

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
        for block, _ in local_constraints:
            g_bounds_values.extend(
                NLPConstraint(f"{block.name}_{i}_{j}", block.exprs[j], block.lb, block.ub, block.nominal)
                for i in range(block.count)
                for j in range(block.output_size)
            )

        local_code = []
        for _, local in local_objectives:
            local_code.append(local.value)
            local_code.append(local.jac if local.jac_mode == "direct" else local.jvp)
            local_code.append(local.hes if local.hes_mode == "direct" else local.hvp)
        for _, local in local_constraints:
            local_code.append(local.value)
            local_code.append(local.jac if local.jac_mode == "direct" else local.jvp)
            local_code.append(local.hes if local.hes_mode == "direct" else local.hvp)

        eval_lines = ["    out[0] = 0.0;"]
        for idx, (block, local) in enumerate(local_objectives):
            eval_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            eval_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            for local_idx, offset in enumerate(block.offsets):
                eval_lines.append(f"        xl[{local_idx}] = x[rep + {offset}];")
            eval_lines.append("        f64 tmp[1];")
            eval_lines.append(f"        moo_nlp_obj_map_{idx}_value(xl, rp, tmp);")
            eval_lines.append("        out[0] += tmp[0];")
            eval_lines.append("    }")

        g_base = 0
        for idx, (block, local) in enumerate(local_constraints):
            eval_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            eval_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            for local_idx, offset in enumerate(block.offsets):
                eval_lines.append(f"        xl[{local_idx}] = x[rep + {offset}];")
            eval_lines.append(f"        f64 tmp[{max(block.output_size, 1)}];")
            eval_lines.append(f"        moo_nlp_g_map_{idx}_value(xl, rp, tmp);")
            for local_row in range(block.output_size):
                eval_lines.append(f"        out[1 + {g_base} + rep * {block.output_size} + {local_row}] = tmp[{local_row}];")
            eval_lines.append("    }")
            g_base += block.count * block.output_size

        jac_lines = [f"    for (int i = 0; i < JAC_BUF_SIZE; ++i) {{ out[i] = 0.0; }}"]
        for idx, (block, local) in enumerate(local_objectives):
            if not local.jac_sparsity:
                continue
            jac_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            jac_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            for local_idx, offset in enumerate(block.offsets):
                jac_lines.append(f"        xl[{local_idx}] = x[rep + {offset}];")
            if local.jac_mode == "direct":
                jac_lines.append(f"        f64 tmp[{max(len(local.jac_sparsity), 1)}];")
                jac_lines.append(f"        moo_nlp_obj_map_{idx}_jvp_sparse(xl, rp, tmp);")
                for local_buf, (row, col) in enumerate(local.jac_sparsity):
                    if row == 0:
                        jac_lines.append(f"        out[{obj_buf_by_col[block.offsets[col]]} + rep] += tmp[{local_buf}];")
            else:
                jac_lines.extend(self._render_local_colored_jac_lines(f"moo_nlp_obj_map_{idx}_jvp", local, block, obj_buf_by_col, "out", "rep", accumulate=True, indent="        "))
            jac_lines.append("    }")

        g_base = 0
        g_buf_cursor = len(obj_jac)
        for idx, (block, local) in enumerate(local_constraints):
            if local.jac_sparsity:
                jac_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
                jac_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
                for local_idx, offset in enumerate(block.offsets):
                    jac_lines.append(f"        xl[{local_idx}] = x[rep + {offset}];")
                if local.jac_mode == "direct":
                    jac_lines.append(f"        f64 tmp[{max(len(local.jac_sparsity), 1)}];")
                    jac_lines.append(f"        moo_nlp_g_map_{idx}_jvp_sparse(xl, rp, tmp);")
                    for local_buf, _ in enumerate(local.jac_sparsity):
                        jac_lines.append(f"        out[{g_buf_cursor} + rep * {len(local.jac_sparsity)} + {local_buf}] = tmp[{local_buf}];")
                else:
                    jac_lines.extend(self._render_local_colored_jac_lines(f"moo_nlp_g_map_{idx}_jvp", local, block, None, "out", "rep", accumulate=False, indent="        ", base_buf=g_buf_cursor))
                jac_lines.append("    }")
                g_buf_cursor += block.count * len(local.jac_sparsity)
            g_base += block.count * block.output_size

        hes_lines = [f"    for (int i = 0; i < HES_SIZE; ++i) {{ out[i] = 0.0; }}", f"    f64 seed[{local_seed_size}] = {{0}};"]
        for idx, (block, local) in enumerate(local_objectives):
            if not local.hes_sparsity:
                continue
            hes_lines.append("    seed[0] = obj_factor;")
            hes_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
            hes_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
            for local_idx, offset in enumerate(block.offsets):
                hes_lines.append(f"        xl[{local_idx}] = x[rep + {offset}];")
            if local.hes_mode == "direct":
                hes_lines.append(f"        f64 tmp[{max(len(local.hes_sparsity), 1)}];")
                hes_lines.append(f"        moo_nlp_obj_map_{idx}_hvp_sparse(xl, seed, rp, tmp);")
                hbuf = hes_buf_table(block, local)
                hes_lines.append(self._static_int_array(f"obj_{idx}_hbuf", hbuf, "        "))
                for local_buf, _ in enumerate(local.hes_sparsity):
                    hes_lines.append(f"        out[obj_{idx}_hbuf[rep * {len(local.hes_sparsity)} + {local_buf}]] += tmp[{local_buf}];")
            else:
                hes_lines.extend(self._render_local_colored_hes_lines(f"moo_nlp_obj_map_{idx}_hvp", local, block, hes_buf_table(block, local), indent="        "))
            hes_lines.append("    }")

        g_base = 0
        for idx, (block, local) in enumerate(local_constraints):
            if local.hes_sparsity:
                hes_lines.append(f"    for (int rep = 0; rep < {block.count}; ++rep) {{")
                hes_lines.append(f"        f64 xl[{max(block.local_size, 1)}];")
                for local_idx, offset in enumerate(block.offsets):
                    hes_lines.append(f"        xl[{local_idx}] = x[rep + {offset}];")
                hbuf = hes_buf_table(block, local)
                hes_lines.append(self._static_int_array(f"g_{idx}_hbuf", hbuf, "        "))
                for local_row in range(block.output_size):
                    hes_lines.append(f"        for (int seed_i = 0; seed_i < {block.output_size}; ++seed_i) {{ seed[seed_i] = 0.0; }}")
                    hes_lines.append(f"        seed[{local_row}] = lambda[{g_base} + rep * {block.output_size} + {local_row}];")
                    if local.hes_mode == "direct":
                        hes_lines.append(f"        f64 tmp_{local_row}[{max(len(local.hes_sparsity), 1)}];")
                        hes_lines.append(f"        moo_nlp_g_map_{idx}_hvp_sparse(xl, seed, rp, tmp_{local_row});")
                        for local_buf, _ in enumerate(local.hes_sparsity):
                            hes_lines.append(f"        out[g_{idx}_hbuf[rep * {len(local.hes_sparsity)} + {local_buf}]] += tmp_{local_row}[{local_buf}];")
                    else:
                        hes_lines.extend(self._render_local_colored_hes_lines(f"moo_nlp_g_map_{idx}_hvp", local, block, hbuf, indent="        ", seed_row=local_row, hbuf_name=f"g_{idx}_hbuf"))
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

    def _render_jac_fill(self, pairs: list[tuple[int, int]]) -> str:
        lines = []
        for idx, (row, col) in enumerate(pairs):
            lines += [f"    v[{col}] = 1.0;", "    moo_nlp_jvp(x, rp, v, tmp);", f"    out[{idx}] = tmp[{row}];", f"    v[{col}] = 0.0;"]
        return "\n".join(lines) or "    (void)out;"

    def _local_color_metadata(self, pairs: list[tuple[int, int]], colors: list[int]) -> tuple[list[int], list[int], list[int], list[int], list[int]]:
        by_color: dict[int, list[tuple[int, int, int]]] = {}
        for idx, (row, col) in enumerate(pairs):
            color = colors[col] if 0 <= col < len(colors) else col
            by_color.setdefault(color, []).append((idx, row, col))
        color_cols: list[int] = []
        color_offsets = [0]
        scatter_idx: list[int] = []
        scatter_row: list[int] = []
        scatter_col: list[int] = []
        scatter_offsets = [0]
        for color in sorted(by_color):
            entries = by_color[color]
            cols = sorted({col for _, _, col in entries})
            color_cols.extend(cols)
            color_offsets.append(len(color_cols))
            for idx, row, col in entries:
                scatter_idx.append(idx)
                scatter_row.append(row)
                scatter_col.append(col)
            scatter_offsets.append(len(scatter_idx))
        return color_offsets, color_cols, scatter_offsets, scatter_idx, scatter_row, scatter_col

    def _static_int_array(self, name: str, values: list[int], indent: str) -> str:
        return f"{indent}static const int {name}[{max(len(values), 1)}] = {{ {', '.join(map(str, values)) or '0'} }};"

    def _render_local_colored_jac_lines(self, fn: str, local: LocalBlockEmission, block: MappedObjectiveBlock | MappedConstraintBlock, obj_buf_by_col: dict[int, int] | None, out_name: str, rep_name: str, accumulate: bool, indent: str, base_buf: int = 0) -> list[str]:
        color_offsets, color_cols, scatter_offsets, scatter_idx, scatter_row, scatter_col = self._local_color_metadata(local.jac_sparsity, local.jac_colors)
        lines = [
            f"{indent}f64 v[{max(block.local_size, 1)}] = {{0}};",
            f"{indent}f64 tmp_color[{max(block.output_size if isinstance(block, MappedConstraintBlock) else 1, 1)}] = {{0}};",
            self._static_int_array("color_offsets", color_offsets, indent),
            self._static_int_array("color_cols", color_cols, indent),
            self._static_int_array("scatter_offsets", scatter_offsets, indent),
            self._static_int_array("scatter_idx", scatter_idx, indent),
            self._static_int_array("scatter_row", scatter_row, indent),
            self._static_int_array("scatter_col", scatter_col, indent),
            f"{indent}for (int color = 0; color < {max(len(color_offsets) - 1, 0)}; ++color) {{",
            f"{indent}    for (int k = color_offsets[color]; k < color_offsets[color + 1]; ++k) {{ v[color_cols[k]] = 1.0; }}",
            f"{indent}    {fn}(xl, rp, v, tmp_color);",
            f"{indent}    for (int k = scatter_offsets[color]; k < scatter_offsets[color + 1]; ++k) {{",
        ]
        if accumulate and obj_buf_by_col is not None:
            offsets = [obj_buf_by_col[offset] for offset in block.offsets]
            lines.append(self._static_int_array("obj_buf_for_local_col", offsets, indent + "        "))
            lines.append(f"{indent}        {out_name}[obj_buf_for_local_col[scatter_col[k]] + {rep_name}] += tmp_color[scatter_row[k]];")
        else:
            lines.append(f"{indent}        {out_name}[{base_buf} + {rep_name} * {len(local.jac_sparsity)} + scatter_idx[k]] = tmp_color[scatter_row[k]];")
        lines.extend([
            f"{indent}    }}",
            f"{indent}    for (int k = color_offsets[color]; k < color_offsets[color + 1]; ++k) {{ v[color_cols[k]] = 0.0; }}",
            f"{indent}}}",
        ])
        return lines

    def _render_local_colored_hes_lines(self, fn: str, local: LocalBlockEmission, block: MappedObjectiveBlock | MappedConstraintBlock, hbuf: list[int], indent: str, seed_row: int = 0, hbuf_name: str | None = None) -> list[str]:
        pairs = list(local.hes_sparsity)
        if not pairs:
            return []
        color_offsets, color_cols, scatter_offsets, scatter_idx, scatter_row, scatter_col = self._local_color_metadata(pairs, local.hes_colors)
        hbuf_name = hbuf_name or "h_buf_for_local"
        lines = [
            f"{indent}{{",
            f"{indent}f64 v_h[{max(block.local_size, 1)}] = {{0}};",
            f"{indent}f64 tmp_h[{max(block.local_size, 1)}] = {{0}};",
            f"{indent}{fn}_cache_t cache;",
            f"{indent}{fn}_prepare(xl, seed, rp, &cache);",
            self._static_int_array("h_color_offsets", color_offsets, indent),
            self._static_int_array("h_color_cols", color_cols, indent),
            self._static_int_array("h_scatter_offsets", scatter_offsets, indent),
            self._static_int_array("h_scatter_idx", scatter_idx, indent),
            self._static_int_array("h_scatter_row", scatter_row, indent),
            f"{indent}for (int color = 0; color < {max(len(color_offsets) - 1, 0)}; ++color) {{",
            f"{indent}    for (int k = h_color_offsets[color]; k < h_color_offsets[color + 1]; ++k) {{ v_h[h_color_cols[k]] = 1.0; }}",
            f"{indent}    {fn}_apply(&cache, v_h, tmp_h);",
        ]
        if hbuf_name == "h_buf_for_local":
            lines.append(self._static_int_array("h_buf_for_local", hbuf, indent))
        lines.extend([
            f"{indent}    for (int k = h_scatter_offsets[color]; k < h_scatter_offsets[color + 1]; ++k) {{ out[{hbuf_name}[rep * {len(local.hes_sparsity)} + h_scatter_idx[k]]] += tmp_h[h_scatter_row[k]]; }}",
            f"{indent}    for (int k = h_color_offsets[color]; k < h_color_offsets[color + 1]; ++k) {{ v_h[h_color_cols[k]] = 0.0; }}",
            f"{indent}}}",
            f"{indent}}}",
        ])
        return lines

    def _render_hes_fill(self, pairs: list[tuple[int, int]]) -> str:
        lines = []
        for idx, (row, col) in enumerate(pairs):
            lines += [f"    v[{col}] = 1.0;", "    moo_nlp_hvp_apply(&cache, v, tmp);", f"    out[{idx}] = tmp[{row}];", f"    v[{col}] = 0.0;"]
        return "\n".join(lines) or "    (void)out;"

    def _render_colored_jac_fill(self, pairs: list[tuple[int, int]], colors: list[int]) -> str:
        return self._render_colored_fill(
            pairs,
            colors,
            "moo_nlp_jvp(x, rp, v, tmp);",
            lambda idx, row, col: f"    out[{idx}] = tmp[{row}];",
        )

    def _render_colored_hes_fill(self, pairs: list[tuple[int, int]], colors: list[int]) -> str:
        return self._render_colored_fill(
            pairs,
            colors,
            "moo_nlp_hvp_apply(&cache, v, tmp);",
            lambda idx, row, col: f"    out[{idx}] = tmp[{row}];",
        )

    def _render_colored_fill(self, pairs: list[tuple[int, int]], colors: list[int], call: str, scatter) -> str:
        if not pairs:
            return "    (void)out;"
        by_color: dict[int, list[tuple[int, int, int]]] = {}
        for idx, (row, col) in enumerate(pairs):
            color = colors[col] if 0 <= col < len(colors) else col
            by_color.setdefault(color, []).append((idx, row, col))
        ordered_colors = sorted(by_color)
        color_cols: list[int] = []
        color_offsets = [0]
        scatter_buf: list[int] = []
        scatter_row: list[int] = []
        scatter_offsets = [0]
        for color in ordered_colors:
            entries = by_color[color]
            cols = sorted({col for _, _, col in entries})
            color_cols.extend(cols)
            color_offsets.append(len(color_cols))
            for idx, row, _ in entries:
                scatter_buf.append(idx)
                scatter_row.append(row)
            scatter_offsets.append(len(scatter_buf))
        return f"""    static const int color_offsets[{len(color_offsets)}] = {{ {', '.join(map(str, color_offsets))} }};
    static const int color_cols[{max(len(color_cols), 1)}] = {{ {', '.join(map(str, color_cols)) or '0'} }};
    static const int scatter_offsets[{len(scatter_offsets)}] = {{ {', '.join(map(str, scatter_offsets))} }};
    static const int scatter_buf[{max(len(scatter_buf), 1)}] = {{ {', '.join(map(str, scatter_buf)) or '0'} }};
    static const int scatter_row[{max(len(scatter_row), 1)}] = {{ {', '.join(map(str, scatter_row)) or '0'} }};
    for (int color = 0; color < {len(ordered_colors)}; ++color) {{
        for (int i = color_offsets[color]; i < color_offsets[color + 1]; ++i) {{ v[color_cols[i]] = 1.0; }}
        {call}
        for (int i = scatter_offsets[color]; i < scatter_offsets[color + 1]; ++i) {{ out[scatter_buf[i]] = tmp[scatter_row[i]]; }}
        for (int i = color_offsets[color]; i < color_offsets[color + 1]; ++i) {{ v[color_cols[i]] = 0.0; }}
    }}"""

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
            f"strategy={self.codegen_strategy}",
            f"structure_strategy={self.codegen_structure_strategy}",
            f"local_strategy={self.codegen_local_strategy}",
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
            jac_mode = derivative_mode(str(self.codegen_local_strategy), emitted.get("JAC_SPARSITY", []), emitted.get("JAC_COLORS", []))
            hes_mode = derivative_mode(str(self.codegen_local_strategy), emitted.get("HES_SPARSITY", []), emitted.get("HES_COLORS", []))
            lines.append("structure_mode=scalar")
        lines.extend([f"jacobian_mode={jac_mode}", f"hessian_mode={hes_mode}"])
        if isinstance(report, dict):
            lines.extend(f"derivative_{key}={value}" for key, value in sorted(report.items()))
        (out_path / "codegen_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
