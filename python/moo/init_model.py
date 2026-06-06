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
from typing import Iterable

from .ad_codegen import emit_function
from .common import SolverMixin, SolverSettings
from .model import Expr, as_expr, _arr, _c_bool, _clean_name, _num
from .results import InitResult, OptimizationRun, read_results


@dataclass
class InitVar:
    name: str
    lb: float = -math.inf
    ub: float = math.inf
    guess: float = 0.0
    nominal: float = 1.0
    base: float = 0.0


@dataclass
class InitConstraint:
    expr: Expr
    lb: float
    ub: float


class InitParameter(Expr):
    def __init__(self, idx: int, base: float):
        super().__init__(f"z[{idx}]", f"z[{idx}]")
        self._base = base

    @property
    def delta(self) -> Expr:
        return self - self._base

    def __iter__(self):
        yield self
        yield self.delta


class InitModel(SolverMixin):
    def __init__(self, name: str):
        self.name = _clean_name(name)
        self.solver_settings = SolverSettings()
        self.variables: list[InitVar] = []
        self.parameters: list[InitVar] = []
        self.runtime_parameters: list[tuple[str, float]] = []
        self.f_constraints: list[Expr] = []
        self.g_constraints: list[InitConstraint] = []
        self.objective: Expr | None = None
        self.derivative_test = False

    @property
    def y_size(self) -> int:
        return len(self.variables)

    @property
    def p_size(self) -> int:
        return len(self.parameters)

    @property
    def z_size(self) -> int:
        return self.y_size + self.p_size

    def add_variable(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> Expr:
        idx = len(self.variables)
        self.variables.append(InitVar(name or f"y{idx}", lb, ub, guess, nominal))
        return Expr(f"z[{idx}]", f"z[{idx}]")

    def add_parameter(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float | None = None, base: float | None = None, nominal: float = 1.0) -> InitParameter:
        if guess is None and base is None:
            guess = 0.0
            base = 0.0
        elif guess is None:
            guess = base
        elif base is None:
            base = guess

        idx = len(self.parameters)
        self.parameters.append(InitVar(name or f"p{idx}", lb, ub, guess, nominal, base))
        return InitParameter(self.y_size + idx, base)

    def delta(self, parameter: InitParameter | Expr) -> Expr:
        if isinstance(parameter, InitParameter):
            return parameter.delta
        raise TypeError("delta(...) expects a parameter returned by add_parameter")

    def add_runtime_parameter(self, name: str, value: float) -> Expr:
        idx = len(self.runtime_parameters)
        self.runtime_parameters.append((name, value))
        return Expr(f"rp[{idx}]", f"rp[{idx}]")

    def set_objective(self, expr: object) -> None:
        self.objective = as_expr(expr)

    def add_f(self, expr: object) -> None:
        self.f_constraints.append(as_expr(expr))

    def add_g(self, expr: object, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        self.g_constraints.append(InitConstraint(as_expr(expr), low, high))

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
        run_cmd = [str(exe)] + self._solver_args(solver, solver_args)
        if capture:
            process = subprocess.run(run_cmd, cwd=cwd, text=True, capture_output=True, check=False)
        else:
            process = subprocess.run(run_cmd, cwd=cwd, text=True, check=False)
        raw_results = read_results(cwd)
        result_view = InitResult(
            raw_results,
            [var.name for var in self.variables],
            [var.name for var in self.parameters],
        )
        return OptimizationRun(process, raw_results, result_view)

    def run(self, out_dir: str | os.PathLike[str], build_dir: str | os.PathLike[str] = "build", repo_root: str | os.PathLike[str] | None = None, solver: str | None = None, solver_args: list[str] | None = None, capture: bool = False, generate: bool = True, run_cwd: str | os.PathLike[str] | None = None) -> OptimizationRun:
        return self.optimize(out_dir, build_dir, repo_root, solver, solver_args, capture, generate, run_cwd)

    def _solver_args(self, solver: str | None, solver_args: list[str] | None) -> list[str]:
        return self.solver_settings.cli_args(solver, solver_args)

    def _validate(self) -> None:
        if self.objective is None:
            raise ValueError("InitModel requires an objective")
        if not self.variables and not self.parameters:
            raise ValueError("InitModel requires at least one variable or parameter")

    def _pairs(self, text: str) -> list[tuple[int, int]]:
        if isinstance(text, list):
            return text
        out = []
        for line in text.splitlines():
            if line.strip():
                a, b = line.split(",", 1)
                out.append((int(a), int(b)))
        return out

    def _emit_ad(self) -> dict[str, str | list[tuple[int, int]]]:
        outputs = [self.objective.lfg if self.objective is not None else "0.0"]
        outputs.extend(expr.lfg for expr in self.f_constraints)
        outputs.extend(c.expr.lfg for c in self.g_constraints)
        emitted = emit_function(
            "z",
            self.z_size,
            outputs,
            [("rp", len(self.runtime_parameters))],
            "moo_init_value",
            "moo_init_jvp",
            "moo_init_hvp",
        )
        return {
            "VALUE": emitted.value,
            "JVP": emitted.jvp,
            "HVP": emitted.hvp,
            "JAC_SPARSITY": emitted.jac_sparsity,
            "HES_SPARSITY": emitted.hes_sparsity,
        }

    def _render_c(self, emitted: dict[str, str], standalone_main: bool) -> str:
        jac = self._pairs(emitted.get("JAC_SPARSITY", ""))
        hes = self._pairs(emitted.get("HES_SPARSITY", ""))
        obj_jac = [(0, col) for row, col in jac if row == 0]
        f_jac = [(row - 1, col) for row, col in jac if 1 <= row < 1 + len(self.f_constraints)]
        g_jac = [(row - 1 - len(self.f_constraints), col) for row, col in jac if row >= 1 + len(self.f_constraints)]

        obj_buf = [idx for idx, (row, _) in enumerate(jac) if row == 0]
        f_buf = [idx for idx, (row, _) in enumerate(jac) if 1 <= row < 1 + len(self.f_constraints)]
        g_buf = [idx for idx, (row, _) in enumerate(jac) if row >= 1 + len(self.f_constraints)]

        def coo(name: str, pairs: list[tuple[int, int]], bufs: list[int] | None = None) -> str:
            return f"""static coo_t {name} = {{
    .row = {_arr([r for r, _ in pairs])},
    .col = {_arr([c for _, c in pairs])},
    .buf_index = {_arr(list(range(len(pairs))) if bufs is None else bufs)},
    .nnz = {len(pairs)}
}};
"""

        def bounds(name: str, values: list[InitVar] | list[InitConstraint]) -> str:
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
#include <string.h>

#include <interfaces/init/main_init.h>
#include "{self.name}.h"

#define Y_SIZE {self.y_size}
#define P_SIZE {self.p_size}
#define Z_SIZE (Y_SIZE + P_SIZE)
#define RP_SIZE {len(self.runtime_parameters)}
#define F_SIZE {len(self.f_constraints)}
#define G_SIZE {len(self.g_constraints)}
#define OUT_SIZE (1 + F_SIZE + G_SIZE)

{emitted.get("VALUE", "")}
{emitted.get("JVP", "")}
{emitted.get("HVP", "")}
static f64 globl_rp[RP_SIZE] = {{ {', '.join(_num(v) for _, v in self.runtime_parameters)} }};
static f64 globl_y0[Y_SIZE] = {{ {', '.join(_num(v.guess) for v in self.variables)} }};
static f64 globl_p0[P_SIZE] = {{ {', '.join(_num(v.base) for v in self.parameters)} }};
static f64 globl_dp0[P_SIZE] = {{ {', '.join(_num(v.guess - v.base) for v in self.parameters)} }};
{bounds("globl_y_bounds", self.variables)}
{bounds("globl_p_bounds", self.parameters)}
static bounds_t globl_dp_bounds[P_SIZE] = {{ {', '.join('{ -DBL_MAX, DBL_MAX }' for _ in self.parameters)} }};
{bounds("globl_g_bounds", self.g_constraints)}
static f64 globl_y_nominal[Y_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.variables)} }};
static f64 globl_p_nominal[P_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.parameters)} }};
static f64 globl_dp_nominal[P_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.parameters)} }};
static f64 globl_f_nominal[F_SIZE] = {{ {', '.join('1.0' for _ in self.f_constraints)} }};
static f64 globl_g_nominal[G_SIZE] = {{ {', '.join('1.0' for _ in self.g_constraints)} }};
{coo("globl_obj_jac", obj_jac, obj_buf)}
{coo("globl_f_jac", f_jac, f_buf)}
{coo("globl_g_jac", g_jac, g_buf)}
{coo("globl_hes", hes)}

static void eval_all(const f64* z, const f64* rp, f64* out, void* user_data) {{
    moo_init_value(z, rp, out);
}}

static void jac_all(const f64* z, const f64* rp, f64* out, void* user_data) {{
    f64 v[Z_SIZE] = {{0}};
    f64 tmp[OUT_SIZE] = {{0}};
{self._render_jac_fill(jac)}
}}

static void hes_all(const f64* z, const f64* rp, const f64* lambda, f64 obj_factor, f64* out, void* user_data) {{
    f64 seed[OUT_SIZE] = {{0}};
    f64 v[Z_SIZE] = {{0}};
    f64 tmp[Z_SIZE] = {{0}};
    moo_init_hvp_cache_t cache;
    seed[0] = obj_factor;
    for (int i = 0; i < F_SIZE + G_SIZE; ++i) {{ seed[1 + i] = lambda[i]; }}
    moo_init_hvp_prepare(z, seed, rp, &cache);
{self._render_hes_fill(hes)}
}}

static c_init_callbacks_t globl_callbacks = {{ eval_all, jac_all, hes_all }};
static solver_ctx_t globl_solver_ctx = {{ .derivative_test = {_c_bool(self.derivative_test)} }};
static c_init_problem_t globl_problem = {{
    .y_size = Y_SIZE,
    .p_size = P_SIZE,
    .rp_size = RP_SIZE,
    .f_size = F_SIZE,
    .g_size = G_SIZE,
    .rp = globl_rp,
    .y0 = globl_y0,
    .p0 = globl_p0,
    .dp0 = globl_dp0,
    .y_bounds = globl_y_bounds,
    .p_bounds = globl_p_bounds,
    .dp_bounds = globl_dp_bounds,
    .g_bounds = globl_g_bounds,
    .y_nominal = globl_y_nominal,
    .p_nominal = globl_p_nominal,
    .dp_nominal = globl_dp_nominal,
    .f_nominal = globl_f_nominal,
    .g_nominal = globl_g_nominal,
    .obj_nominal = 1.0,
    .obj_jac = &globl_obj_jac,
    .f_jac = &globl_f_jac,
    .g_jac = &globl_g_jac,
    .hes = &globl_hes,
    .callbacks = &globl_callbacks,
    .solver_ctx = &globl_solver_ctx,
    .user_data = (void*)0
}};

int main_{self.name}(int argc, char** argv) {{
    return main_init(argc, argv, &globl_problem);
}}
{main}
"""

    def _render_jac_fill(self, pairs: list[tuple[int, int]]) -> str:
        lines = []
        for idx, (row, col) in enumerate(pairs):
            lines += [f"    v[{col}] = 1.0;", "    moo_init_jvp(z, rp, v, tmp);", f"    out[{idx}] = tmp[{row}];", f"    v[{col}] = 0.0;"]
        return "\n".join(lines) or "    (void)out;"

    def _render_hes_fill(self, pairs: list[tuple[int, int]]) -> str:
        lines = []
        for idx, (row, col) in enumerate(pairs):
            lines += [f"    v[{col}] = 1.0;", "    moo_init_hvp_apply(&cache, v, tmp);", f"    out[{idx}] = tmp[{row}];", f"    v[{col}] = 0.0;"]
        return "\n".join(lines) or "    (void)out;"

    def _render_h(self) -> str:
        guard = f"MOO_INIT_CODEGEN_{self.name.upper()}_H"
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
