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

from .ad_codegen import emit_function
from .common import SolverMixin, SolverSettings
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


class NLPModel(SolverMixin):
    def __init__(self, name: str):
        self.name = _clean_name(name)
        self.solver_settings = SolverSettings()
        self.variables: list[NLPVar] = []
        self.constraints: list[NLPConstraint] = []
        self.runtime_parameters: list[tuple[str, float]] = []
        self.objective: Expr | None = None
        self.obj_nominal = 1.0
        self.derivative_test = False

    @property
    def x_size(self) -> int:
        return len(self.variables)

    def add_variable(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> Expr:
        idx = len(self.variables)
        self.variables.append(NLPVar(name or f"x{idx}", lb, ub, guess, nominal))
        return Expr(f"x[{idx}]", f"x[{idx}]")

    def add_runtime_parameter(self, name: str, value: float) -> Expr:
        idx = len(self.runtime_parameters)
        self.runtime_parameters.append((name, value))
        return Expr(f"rp[{idx}]", f"rp[{idx}]")

    def minimize(self, expr: object, nominal: float = 1.0) -> None:
        self.objective = as_expr(expr)
        self.obj_nominal = nominal

    def add_constraint(self, expr: object, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None, name: str | None = None, nominal: float = 1.0) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        idx = len(self.constraints)
        self.constraints.append(NLPConstraint(name or f"g{idx}", as_expr(expr), low, high, nominal))

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
            [constraint.name for constraint in self.constraints],
        )
        return OptimizationRun(process, raw_results, result_view)

    def run(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, solver: str | None = None, solver_args: list[str] | None = None, capture: bool = False, generate: bool = True, run_cwd: str | os.PathLike[str] | None = None) -> OptimizationRun:
        return self.optimize(out_dir, build_dir, repo_root, solver, solver_args, capture, generate, run_cwd)

    def _solver_args(self, solver: str | None, solver_args: list[str] | None) -> list[str]:
        return self.solver_settings.cli_args(solver, solver_args)

    def _validate(self) -> None:
        if not self.variables:
            raise ValueError("NLPModel requires at least one variable")
        if self.objective is None:
            raise ValueError("NLPModel requires an objective; call minimize(...)")

    def _emit_ad(self) -> dict[str, str | list[tuple[int, int]]]:
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
            "JAC_SPARSITY": emitted.jac_sparsity,
            "HES_SPARSITY": emitted.hes_sparsity,
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
        jac = self._pairs(emitted.get("JAC_SPARSITY", ""))
        hes = self._pairs(emitted.get("HES_SPARSITY", ""))
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

        return f"""#include <float.h>
#include <stdbool.h>

#include <interfaces/nlp/main_nlp.h>
#include "{self.name}.h"

#define X_SIZE {self.x_size}
#define RP_SIZE {len(self.runtime_parameters)}
#define G_SIZE {len(self.constraints)}
#define OUT_SIZE (1 + G_SIZE)

{emitted.get("VALUE", "")}
{emitted.get("JVP", "")}
{emitted.get("HVP", "")}
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
    f64 v[X_SIZE] = {{0}};
    f64 tmp[OUT_SIZE] = {{0}};
{self._render_jac_fill(jac)}
}}

static void hes_all(const f64* x, const f64* rp, const f64* lambda, f64 obj_factor, f64* out, void* user_data) {{
    f64 seed[OUT_SIZE] = {{0}};
    f64 v[X_SIZE] = {{0}};
    f64 tmp[X_SIZE] = {{0}};
    moo_nlp_hvp_cache_t cache;
    seed[0] = obj_factor;
    for (int i = 0; i < G_SIZE; ++i) {{ seed[1 + i] = lambda[i]; }}
    moo_nlp_hvp_prepare(x, seed, rp, &cache);
{self._render_hes_fill(hes)}
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

    def _render_hes_fill(self, pairs: list[tuple[int, int]]) -> str:
        lines = []
        for idx, (row, col) in enumerate(pairs):
            lines += [f"    v[{col}] = 1.0;", "    moo_nlp_hvp_apply(&cache, v, tmp);", f"    out[{idx}] = tmp[{row}];", f"    v[{col}] = 0.0;"]
        return "\n".join(lines) or "    (void)out;"

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
