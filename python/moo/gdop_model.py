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
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .callback_codegen import render_hessian_callback_body, render_jacobian_callback_body
from .common import select_derivative_callback_mode, parse_sparsity_pairs
from .expressions import Expr, as_expr, cos, exp, log, pow_const, sin, sum_expr, tan
from .local_function import InputGroup, LocalGraphFunction, LocalValueFunction
from .model import BaseModel, _arr, _c_bool, _double_arr, _num
from . import paths
from .results import GDOPResult, OptimizationRun, read_results


@dataclass
class Variable:
    name: str
    lb: float = -math.inf
    ub: float = math.inf
    guess: float = 0.0
    nominal: float = 1.0
    start: float | None = None
    final: float | None = None


@dataclass
class Constraint:
    expr: Expr
    lb: float
    ub: float


class GDOPModel(BaseModel):
    def __init__(self, name: str):
        super().__init__(name)
        self.states: list[Variable] = []
        self.controls: list[Variable] = []
        self.parameters: list[Variable] = []
        self.dynamics: list[Expr | None] = []
        self.path_constraints: list[Constraint] = []
        self.boundary_constraints: list[Constraint] = []
        self.mayer: Expr | None = None
        self.lagrange: Expr | None = None
        self.t0_bounds = (0.0, 0.0)
        self.tf_bounds = (1.0, 1.0)
        self.mesh_intervals = 25
        self.mesh_nodes = 3
        self.mesh_l2bn_p1_it = 0
        self.mesh_l2bn_p2_it = 0
        self.mesh_l2bn_p2_lvl = 0.0

    @property
    def x_size(self) -> int:
        return len(self.states)

    @property
    def u_size(self) -> int:
        return len(self.controls)

    @property
    def p_size(self) -> int:
        return len(self.parameters)

    @property
    def z_size(self) -> int:
        return self.x_size + self.u_size + self.p_size

    @property
    def mr_size(self) -> int:
        return 2 * (self.x_size + self.u_size) + self.p_size + 2

    def set_time_fixed(self, t0: float = 0.0, tf: float = 1.0) -> None:
        self.t0_bounds = (t0, t0)
        self.tf_bounds = (tf, tf)

    def set_time_free(
        self,
        t0_guess: float = 0.0,
        tf_guess: float = 1.0,
        t0_bounds: tuple[float, float] = (-math.inf, math.inf),
        tf_bounds: tuple[float, float] = (-math.inf, math.inf),
    ) -> None:
        self.t0_bounds = t0_bounds
        self.tf_bounds = tf_bounds
        self.t0_guess = t0_guess
        self.tf_guess = tf_guess

    def mesh(self, intervals: int = 25, nodes: int = 3):
        self.mesh_intervals = intervals
        self.mesh_nodes = nodes
        return self

    def mesh_refinement(self, l2bn_p1_it: int = 0, l2bn_p2_it: int = 0, l2bn_p2_lvl: float = 0.0):
        self.mesh_l2bn_p1_it = l2bn_p1_it
        self.mesh_l2bn_p2_it = l2bn_p2_it
        self.mesh_l2bn_p2_lvl = l2bn_p2_lvl
        return self

    def add_state(self, name: str | None = None, start: float = 0.0, lb: float = -math.inf, ub: float = math.inf, final: float | None = None, nominal: float = 1.0) -> Expr:
        idx = len(self.states)
        self.states.append(Variable(name or f"x{idx}", lb, ub, start, nominal, start, final))
        self.dynamics.append(None)
        return Expr.variable("z", idx, mr_group="xf", mr_idx=idx)

    def initial(self, var: Expr) -> Expr:
        symbol = var.mr_symbol
        if symbol is None or symbol.group not in {"xf", "uf"}:
            raise ValueError("initial() currently expects a single state/control variable")
        group = "x0" if symbol.group == "xf" else "u0"
        return Expr.variable(group, int(symbol.index))

    def add_control(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> Expr:
        idx = len(self.controls)
        z_idx = self.x_size + idx
        self.controls.append(Variable(name or f"u{idx}", lb, ub, guess, nominal))
        return Expr.variable("z", z_idx, mr_group="uf", mr_idx=idx)

    def add_parameter(self, name: str | None = None, lb: float = -math.inf, ub: float = math.inf, guess: float = 0.0, nominal: float = 1.0) -> Expr:
        idx = len(self.parameters)
        z_idx = self.x_size + self.u_size + idx
        self.parameters.append(Variable(name or f"p{idx}", lb, ub, guess, nominal))
        return Expr.variable("z", z_idx, mr_group="p", mr_idx=idx)

    def add_runtime_parameter(self, name: str, value: float) -> Expr:
        idx = len(self.runtime_parameters)
        self.runtime_parameters.append((name, value))
        return Expr.variable("rp", idx)

    @property
    def time(self) -> Expr:
        return Expr.variable("tau", 0)

    @property
    def t0(self) -> Expr:
        return Expr.variable("t0", 0)

    @property
    def tf(self) -> Expr:
        return Expr.variable("tf", 0)

    def set_dynamics(self, state: Expr, rhs: object) -> None:
        symbol = state.symbol
        if symbol is None or symbol.group != "z":
            raise ValueError("set_dynamics expects a state variable returned by add_state")
        idx = int(symbol.index)
        if idx >= self.x_size:
            raise ValueError("set_dynamics expects a state variable")
        self.dynamics[idx] = as_expr(rhs)

    def add_lagrange(self, expr: object) -> None:
        self.lagrange = as_expr(expr)

    def add_mayer(self, expr: object) -> None:
        self.mayer = as_expr(expr)

    def add_path(self, expr: object, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        self.path_constraints.append(Constraint(as_expr(expr), low, high))

    def add_boundary(self, expr: object, lb: float = -math.inf, ub: float = math.inf, eq: float | None = None) -> None:
        low, high = (eq, eq) if eq is not None else (lb, ub)
        self.boundary_constraints.append(Constraint(as_expr(expr), low, high))

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

    def compile(
        self,
        out_dir: str | os.PathLike[str],
        build_dir: str | os.PathLike[str] = "build",
        repo_root: str | os.PathLike[str] | None = None,
        generate: bool = False,
    ) -> Path:
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
        cmd = [
            cc,
            "-std=c99", "-O3",
            f"-I{include}",
            f"-I{abs_out_path}",
            str(c_path),
            f"-L{lib_dir}",
            "-lmoo",
            f"-Wl,-rpath,{lib_dir}",
            "-lm",
            "-o",
            str(exe),
        ]
        subprocess.run(cmd, check=True)
        return exe

    def optimize(
        self,
        out_dir: str | os.PathLike[str],
        build_dir: str | os.PathLike[str] = "build",
        repo_root: str | os.PathLike[str] | None = None,
        solver: str | None = None,
        solver_args: list[str] | None = None,
        capture: bool = False,
        generate: bool = False,
        run_cwd: str | os.PathLike[str] | None = None
    ) -> OptimizationRun:
        exe = self.compile(out_dir, build_dir=build_dir, repo_root=repo_root, generate=generate)
        cwd = Path(run_cwd) if run_cwd is not None else paths.library_dir(build_dir)
        if not cwd.is_absolute():
            cwd = cwd.resolve()
        cwd.mkdir(parents=True, exist_ok=True)
        run_cmd = [str(exe)] + self._solver_args(solver, solver_args)
        if capture:
            result = subprocess.run(run_cmd, cwd=cwd, text=True, capture_output=True, check=False)
        else:
            result = subprocess.run(run_cmd, cwd=cwd, text=True, check=False)
        raw_results = read_results(cwd)
        result_view = GDOPResult(
            raw_results,
            [var.name for var in self.states],
            [var.name for var in self.controls],
        )
        return OptimizationRun(result, raw_results, result_view)

    def run(
        self,
        out_dir: str | os.PathLike[str],
        build_dir: str | os.PathLike[str] = "build",
        repo_root: str | os.PathLike[str] | None = None,
        solver: str | None = None,
        solver_args: list[str] | None = None,
        capture: bool = False,
        generate: bool = True,
        run_cwd: str | os.PathLike[str] | None = None
    ) -> OptimizationRun:
        return self.optimize(
            out_dir,
            build_dir=build_dir,
            repo_root=repo_root,
            solver=solver,
            solver_args=solver_args,
            capture=capture,
            generate=generate,
            run_cwd=run_cwd
        )

    def _solver_args(self, solver: str | None, solver_args: list[str] | None) -> list[str]:
        return self.solver_settings.cli_args(solver, solver_args)

    def _validate(self) -> None:
        if not self.states:
            raise ValueError("A GDOP needs at least one state")
        missing = [idx for idx, rhs in enumerate(self.dynamics) if rhs is None]
        if missing:
            raise ValueError(f"Missing dynamics for state indices: {missing}")

    def _resolve_lfg(self, expr: Expr | None) -> Expr | None:
        return expr

    def _mr_offset(self, group: str, idx: int) -> int:
        offsets = {
            "x0": 0,
            "u0": self.x_size,
            "xf": self.x_size + self.u_size,
            "uf": 2 * self.x_size + self.u_size,
            "p": 2 * (self.x_size + self.u_size),
            "t0": 2 * (self.x_size + self.u_size) + self.p_size,
            "tf": 2 * (self.x_size + self.u_size) + self.p_size + 1,
        }
        return offsets[group] + (0 if group in {"t0", "tf"} else idx)

    def _map_mr_variable(self, group: str, idx: int) -> tuple[str, int] | None:
        if group == "rp":
            return None
        if group in {"x0", "u0", "xf", "uf", "p", "t0", "tf"}:
            return ("b", self._mr_offset(group, idx))
        raise ValueError(f"Unsupported GDOP boundary variable group {group!r}")

    def _resolve_mr_expr(self, expr: Expr | None) -> Expr | None:
        if expr is None:
            return None
        return expr.remap_mr(self._map_mr_variable)

    def _emit_ad(self) -> dict[str, str | list[tuple[int, int]]]:
        lfg_outputs = []
        if self.lagrange is not None:
            lfg_outputs.append(self._resolve_lfg(self.lagrange))
        lfg_outputs.extend(self._resolve_lfg(rhs) for rhs in self.dynamics if rhs is not None)
        lfg_outputs.extend(self._resolve_lfg(c.expr) for c in self.path_constraints)

        mr_outputs = []
        if self.mayer is not None:
            mr_outputs.append(self._resolve_mr_expr(self.mayer))
        mr_outputs.extend(self._resolve_mr_expr(c.expr) for c in self.boundary_constraints)

        ode_outputs = [self._resolve_lfg(rhs) for rhs in self.dynamics if rhs is not None]

        lfg = LocalGraphFunction(
            name="moo_lfg",
            input=InputGroup("z", self.z_size),
            outputs=lfg_outputs,
            params=[InputGroup("rp", len(self.runtime_parameters)), InputGroup("tau", 1)],
        ).emit()
        mr = LocalGraphFunction(
            name="moo_mr",
            input=InputGroup("b", self.mr_size),
            outputs=mr_outputs,
            params=[InputGroup("rp", len(self.runtime_parameters))],
        ).emit()
        ode_value, ode_jac = LocalValueFunction(
            name="moo_ode",
            input=InputGroup("z", self.z_size),
            outputs=ode_outputs,
            params=[InputGroup("rp", len(self.runtime_parameters)), InputGroup("tau", 1)],
        ).emit()
        lfg_jac_mode = select_derivative_callback_mode(self.derivative_strategy, lfg.jac_sparsity, lfg.jac_colors)
        lfg_hes_mode = select_derivative_callback_mode(self.derivative_strategy, lfg.hes_sparsity, lfg.hes_colors)
        mr_jac_mode = select_derivative_callback_mode(self.derivative_strategy, mr.jac_sparsity, mr.jac_colors)
        mr_hes_mode = select_derivative_callback_mode(self.derivative_strategy, mr.hes_sparsity, mr.hes_colors)
        return {
            "LFG_VALUE": lfg.value,
            "LFG_JVP": lfg.jvp,
            "LFG_HVP": lfg.hvp,
            "LFG_JAC": lfg.jac,
            "LFG_HES": lfg.hes,
            "LFG_JAC_SPARSITY": lfg.jac_sparsity,
            "LFG_HES_SPARSITY": lfg.hes_sparsity,
            "LFG_JAC_COLORS": lfg.jac_colors,
            "LFG_HES_COLORS": lfg.hes_colors,
            "LFG_JAC_MODE": lfg_jac_mode,
            "LFG_HES_MODE": lfg_hes_mode,
            "LFG_REPORT": lfg.report,
            "MR_VALUE": mr.value,
            "MR_JVP": mr.jvp,
            "MR_HVP": mr.hvp,
            "MR_JAC": mr.jac,
            "MR_HES": mr.hes,
            "MR_JAC_SPARSITY": mr.jac_sparsity,
            "MR_HES_SPARSITY": mr.hes_sparsity,
            "MR_JAC_COLORS": mr.jac_colors,
            "MR_HES_COLORS": mr.hes_colors,
            "MR_JAC_MODE": mr_jac_mode,
            "MR_HES_MODE": mr_hes_mode,
            "MR_REPORT": mr.report,
            "ODE_VALUE": ode_value,
            "ODE_JAC_SPARSITY": ode_jac,
        }

    def _render_c(self, emitted: dict[str, str], standalone_main: bool) -> str:
        lfg_eval_count = (1 if self.lagrange is not None else 0) + self.x_size + len(self.path_constraints)
        mr_eval_count = (1 if self.mayer is not None else 0) + len(self.boundary_constraints)
        lfg_jac = parse_sparsity_pairs(emitted.get("LFG_JAC_SPARSITY", ""))
        lfg_hes = parse_sparsity_pairs(emitted.get("LFG_HES_SPARSITY", ""))
        mr_jac = parse_sparsity_pairs(emitted.get("MR_JAC_SPARSITY", ""))
        mr_hes = parse_sparsity_pairs(emitted.get("MR_HES_SPARSITY", ""))
        lfg_jac_colors = emitted.get("LFG_JAC_COLORS", [])
        lfg_hes_colors = emitted.get("LFG_HES_COLORS", [])
        mr_jac_colors = emitted.get("MR_JAC_COLORS", [])
        mr_hes_colors = emitted.get("MR_HES_COLORS", [])
        ode_jac = [(row, col) for row, col in parse_sparsity_pairs(emitted.get("ODE_JAC_SPARSITY", "")) if col < self.x_size]

        def coo(name: str, pairs: list[tuple[int, int]], buf_indices: list[int] | None = None) -> str:
            rows = [r for r, _ in pairs]
            cols = [c for _, c in pairs]
            bufs = list(range(len(pairs))) if buf_indices is None else buf_indices
            return f"""static coo_t {name} = {{
    .row = {_arr(rows)},
    .col = {_arr(cols)},
    .buf_index = {_arr(bufs)},
    .nnz = {len(pairs)}
}};
"""

        lfg_hes_buf_indices = []
        lfg_hes_non_pp = 0
        lfg_hes_pp = 0
        for row, col in lfg_hes:
            if row >= self.x_size + self.u_size and col >= self.x_size + self.u_size:
                lfg_hes_buf_indices.append(lfg_hes_pp)
                lfg_hes_pp += 1
            else:
                lfg_hes_buf_indices.append(lfg_hes_non_pp)
                lfg_hes_non_pp += 1

        def bounds(name: str, values: list[Variable] | list[Constraint] | list[tuple[float, float]]) -> str:
            if not values:
                return f"static bounds_t {name}[1];\n"
            items = []
            for v in values:
                if isinstance(v, Variable):
                    items.append(f"{{ {_num(v.lb)}, {_num(v.ub)} }}")
                elif isinstance(v, Constraint):
                    items.append(f"{{ {_num(v.lb)}, {_num(v.ub)} }}")
                else:
                    items.append(f"{{ {_num(v[0])}, {_num(v[1])} }}")
            return f"static bounds_t {name}[{len(items)}] = {{ {', '.join(items)} }};\n"

        xu0 = [f"{{ {_num(v.start)}, true }}" for v in self.states] + [f"{{ 0.0, false }}" for _ in self.controls]
        xuf = [f"{{ {_num(v.final)}, {_c_bool(v.final is not None)} }}" for v in self.states] + [f"{{ 0.0, false }}" for _ in self.controls]
        rp_values = [value for _, value in self.runtime_parameters]

        lfg_lambda_size = self.x_size + len(self.path_constraints)
        mr_lambda_size = len(self.boundary_constraints)
        lfg_jac_mode = str(emitted.get("LFG_JAC_MODE", "direct"))
        lfg_hes_mode = str(emitted.get("LFG_HES_MODE", "direct"))
        mr_jac_mode = str(emitted.get("MR_JAC_MODE", "direct"))
        mr_hes_mode = str(emitted.get("MR_HES_MODE", "direct"))

        section_keys = ["LFG_VALUE"]
        section_keys.append("LFG_JAC" if lfg_jac_mode == "direct" else "LFG_JVP")
        section_keys.append("LFG_HES" if lfg_hes_mode == "direct" else "LFG_HVP")
        section_keys.append("MR_VALUE")
        section_keys.append("MR_JAC" if mr_jac_mode == "direct" else "MR_JVP")
        section_keys.append("MR_HES" if mr_hes_mode == "direct" else "MR_HVP")
        section_keys.append("ODE_VALUE")
        deduped_section_keys = list(dict.fromkeys(section_keys))
        code_sections = "\n".join(emitted.get(key, "") for key in deduped_section_keys)
        main = f"""
int main(int argc, char** argv) {{
    return main_{self.name}(argc, argv);
}}
""" if standalone_main else ""
        lfg_jac_body = render_jacobian_callback_body(
            lfg_jac_mode,
            "moo_lfg_jvp_sparse(z, globl_rp, tau, out);",
            "Z_SIZE",
            str(max(lfg_eval_count, 1)),
            lfg_jac,
            lfg_jac_colors if isinstance(lfg_jac_colors, list) else [],
            "moo_lfg_jvp(z, globl_rp, tau, v, tmp);",
        )
        lfg_hes_body = render_hessian_callback_body(
            lfg_hes_mode,
            f"""    f64 tmp[{max(len(lfg_hes), 1)}] = {{0}};
    moo_lfg_hvp_sparse(z, seed, globl_rp, tau, tmp);
{self._render_sparse_hes_scatter(lfg_hes, lfg_hes_buf_indices, True)}""",
            "Z_SIZE",
            "Z_SIZE",
            lfg_hes,
            lfg_hes_colors if isinstance(lfg_hes_colors, list) else [],
            """    moo_lfg_hvp_cache_t cache;
    moo_lfg_hvp_prepare(z, seed, globl_rp, tau, &cache);""",
            "moo_lfg_hvp_apply(&cache, v, tmp);",
            buf_indices=lfg_hes_buf_indices,
            split_pp=True,
            pp_start=self.x_size + self.u_size,
        )
        mr_jac_body = render_jacobian_callback_body(
            mr_jac_mode,
            "moo_mr_jvp_sparse(b, globl_rp, out);",
            str(max(self.mr_size, 1)),
            str(max(mr_eval_count, 1)),
            mr_jac,
            mr_jac_colors if isinstance(mr_jac_colors, list) else [],
            "moo_mr_jvp(b, globl_rp, v, tmp);",
        )
        mr_hes_body = render_hessian_callback_body(
            mr_hes_mode,
            "    moo_mr_hvp_sparse(b, seed, globl_rp, out);",
            str(max(self.mr_size, 1)),
            str(max(self.mr_size, 1)),
            mr_hes,
            mr_hes_colors if isinstance(mr_hes_colors, list) else [],
            """    moo_mr_hvp_cache_t cache;
    moo_mr_hvp_prepare(b, seed, globl_rp, &cache);""",
            "moo_mr_hvp_apply(&cache, v, tmp);",
            buf_indices=list(range(len(mr_hes))),
        )

        return f"""#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include <interfaces/gdop/structures.h>
#include <interfaces/gdop/main_gdop.h>
#include "{self.name}.h"

#define X_SIZE {self.x_size}
#define U_SIZE {self.u_size}
#define P_SIZE {self.p_size}
#define XU_SIZE (X_SIZE + U_SIZE)
#define Z_SIZE (X_SIZE + U_SIZE + P_SIZE)
#define RP_SIZE {len(self.runtime_parameters)}
#define R_SIZE {len(self.boundary_constraints)}
#define G_SIZE {len(self.path_constraints)}
#define HAS_MAYER {_c_bool(self.mayer is not None)}
#define HAS_LAGRANGE {_c_bool(self.lagrange is not None)}

{code_sections}
{bounds("globl_x_bounds", self.states)}
{bounds("globl_u_bounds", self.controls)}
{bounds("globl_p_bounds", self.parameters)}
{bounds("globl_T_bounds", [self.t0_bounds, self.tf_bounds])}
{bounds("globl_g_bounds", self.path_constraints)}
{bounds("globl_r_bounds", self.boundary_constraints)}
static optional_value_t globl_xu0_fixed[XU_SIZE] = {{ {', '.join(xu0)} }};
static optional_value_t globl_xuf_fixed[XU_SIZE] = {{ {', '.join(xuf)} }};
static f64 globl_x_nominal[X_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.states)} }};
static f64 globl_u_nominal[U_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.controls)} }};
static f64 globl_p_nominal[P_SIZE] = {{ {', '.join(_num(v.nominal) for v in self.parameters)} }};
static f64 globl_obj_nominal = 1.0;
static f64 globl_f_nominal[X_SIZE] = {{ {', '.join('1.0' for _ in self.states)} }};
static f64 globl_g_nominal[G_SIZE] = {{ {', '.join('1.0' for _ in self.path_constraints)} }};
static f64 globl_r_nominal[R_SIZE] = {{ {', '.join('1.0' for _ in self.boundary_constraints)} }};
static f64 globl_rp[RP_SIZE] = {{ {', '.join(_num(v) for v in rp_values)} }};
static const char* data_files[1];

static eval_structure_t globl_lfg_eval = {{ .buf_index = {_arr(range(lfg_eval_count))} }};
{coo("globl_lfg_jac", lfg_jac)}
{coo("globl_lfg_lt_hes", lfg_hes, lfg_hes_buf_indices)}
static eval_structure_t globl_mr_eval = {{ .buf_index = {_arr(range(mr_eval_count))} }};
{coo("globl_mr_jac", mr_jac)}
{coo("globl_mr_lt_hes", mr_hes)}
{coo("globl_ode_jac", ode_jac)}

static void pack_z(const f64* xu, const f64* p, f64* z) {{
    memcpy(z, xu, XU_SIZE * sizeof(f64));
    memcpy(z + XU_SIZE, p, P_SIZE * sizeof(f64));
}}

static void eval_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {{
    f64 z[Z_SIZE];
    f64 tau[1] = {{ t }};
    pack_z(xu, p, z);
    moo_lfg_value(z, globl_rp, tau, out);
}}

static void jac_lfg(const f64* xu, const f64* p, f64 t, const f64* data, f64* out, void* user_data) {{
    f64 z[Z_SIZE];
    f64 tau[1] = {{ t }};
    pack_z(xu, p, z);
{lfg_jac_body}
}}

static void hes_lfg(const f64* xu, const f64* p, const f64* lambda, const f64 obj_factor, f64 t, const f64* data, f64* out, f64* out_pp, void* user_data) {{
    f64 z[Z_SIZE];
    f64 tau[1] = {{ t }};
    f64 seed[{max(lfg_eval_count, 1)}] = {{0}};
    pack_z(xu, p, z);
    int seed_offset = 0;
    if (HAS_LAGRANGE) {{ seed[0] = obj_factor; seed_offset = 1; }}
    for (int i = 0; i < {lfg_lambda_size}; ++i) {{ seed[seed_offset + i] = lambda[i]; }}
{lfg_hes_body}
}}

static void eval_mr(const f64* xu0, const f64* xuf, const f64* p, f64 t0, f64 tf, const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {{
    if ({mr_eval_count} == 0) {{ (void)out; return; }}
    f64 b[{max(self.mr_size, 1)}];
    memcpy(b, xu0, XU_SIZE * sizeof(f64));
    memcpy(b + XU_SIZE, xuf, XU_SIZE * sizeof(f64));
    memcpy(b + 2 * XU_SIZE, p, P_SIZE * sizeof(f64));
    b[2 * XU_SIZE + P_SIZE] = t0;
    b[2 * XU_SIZE + P_SIZE + 1] = tf;
    moo_mr_value(b, globl_rp, out);
}}

static void jac_mr(const f64* xu0, const f64* xuf, const f64* p, f64 t0, f64 tf, const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {{
    f64 b[{max(self.mr_size, 1)}];
    memcpy(b, xu0, XU_SIZE * sizeof(f64));
    memcpy(b + XU_SIZE, xuf, XU_SIZE * sizeof(f64));
    memcpy(b + 2 * XU_SIZE, p, P_SIZE * sizeof(f64));
    b[2 * XU_SIZE + P_SIZE] = t0;
    b[2 * XU_SIZE + P_SIZE + 1] = tf;
{mr_jac_body}
}}

static void hes_mr(const f64* xu0, const f64* xuf, const f64* p, const f64* lambda, const f64 obj_factor, f64 t0, f64 tf, const f64* data_t0, const f64* data_tf, f64* out, void* user_data) {{
    f64 b[{max(self.mr_size, 1)}];
    f64 seed[{max(mr_eval_count, 1)}] = {{0}};
    memcpy(b, xu0, XU_SIZE * sizeof(f64));
    memcpy(b + XU_SIZE, xuf, XU_SIZE * sizeof(f64));
    memcpy(b + 2 * XU_SIZE, p, P_SIZE * sizeof(f64));
    b[2 * XU_SIZE + P_SIZE] = t0;
    b[2 * XU_SIZE + P_SIZE + 1] = tf;
    int seed_offset = 0;
    if (HAS_MAYER) {{ seed[0] = obj_factor; seed_offset = 1; }}
    for (int i = 0; i < {mr_lambda_size}; ++i) {{ seed[seed_offset + i] = lambda[i]; }}
{mr_hes_body}
}}

static void ode_eval_f(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* f, void* user_data) {{
    f64 z[Z_SIZE];
    f64 tau[1] = {{ t }};
    memcpy(z, x, X_SIZE * sizeof(f64));
    memcpy(z + X_SIZE, u, U_SIZE * sizeof(f64));
    memcpy(z + XU_SIZE, p, P_SIZE * sizeof(f64));
    moo_ode_value(z, globl_rp, tau, f);
}}

static void ode_jac_f(const f64* x, const f64* u, const f64* p, f64 t, const f64* data, f64* dfdx, void* user_data) {{
    f64 z[Z_SIZE];
    f64 tau[1] = {{ t }};
    f64 all[{max(len(lfg_jac), 1)}];
    f64 xu[XU_SIZE];
    memcpy(xu, x, X_SIZE * sizeof(f64));
    memcpy(xu + X_SIZE, u, U_SIZE * sizeof(f64));
    pack_z(xu, p, z);
    jac_lfg(xu, p, t, data, all, user_data);
{self._render_ode_jac_copy(ode_jac, lfg_jac)}
}}

static c_callbacks_t globl_callbacks = {{
    eval_lfg,
    jac_lfg,
    hes_lfg,
    eval_mr,
    jac_mr,
    hes_mr,
    ode_eval_f,
    ode_jac_f
}};

static mesh_ref_ctx_t globl_mesh_ctx = {{
    .initial_intervals = {self.mesh_intervals},
    .nodes_per_interval = {self.mesh_nodes},
    .l2bn_p1_it = {self.mesh_l2bn_p1_it},
    .l2bn_p2_it = {self.mesh_l2bn_p2_it},
    .l2bn_p2_lvl = {_num(self.mesh_l2bn_p2_lvl)}
}};

static solver_ctx_t globl_solver_ctx = {{ .derivative_test = {_c_bool(self.derivative_test)} }};

static c_problem_t globl_c_problem = {{
    .x_size = X_SIZE,
    .u_size = U_SIZE,
    .xu_size = XU_SIZE,
    .p_size = P_SIZE,
    .rp_size = RP_SIZE,
    .r_size = R_SIZE,
    .g_size = G_SIZE,
    .has_mayer = HAS_MAYER,
    .has_lagrange = HAS_LAGRANGE,
    .data_filepath = data_files,
    .data_file_count = 0,
    .rp = globl_rp,
    .x_bounds = globl_x_bounds,
    .u_bounds = globl_u_bounds,
    .p_bounds = globl_p_bounds,
    .T_bounds = globl_T_bounds,
    .r_bounds = globl_r_bounds,
    .g_bounds = globl_g_bounds,
    .xu0_fixed = globl_xu0_fixed,
    .xuf_fixed = globl_xuf_fixed,
    .x_nominal = globl_x_nominal,
    .u_nominal = globl_u_nominal,
    .p_nominal = globl_p_nominal,
    .obj_nominal = &globl_obj_nominal,
    .f_nominal = globl_f_nominal,
    .g_nominal = globl_g_nominal,
    .r_nominal = globl_r_nominal,
    .lfg_eval = &globl_lfg_eval,
    .lfg_jac = &globl_lfg_jac,
    .lfg_lt_hes = &globl_lfg_lt_hes,
    .mr_eval = &globl_mr_eval,
    .mr_jac = &globl_mr_jac,
    .mr_lt_hes = &globl_mr_lt_hes,
    .ode_jac = &globl_ode_jac,
    .callbacks = &globl_callbacks,
    .mesh_ctx = &globl_mesh_ctx,
    .solver_ctx = &globl_solver_ctx,
    .user_data = (void*)0
}};

int main_{self.name}(int argc, char** argv) {{
    return main_gdop(argc, argv, &globl_c_problem);
}}
{main}
"""

    def _render_sparse_hes_scatter(self, pairs: list[tuple[int, int]], buf_indices: list[int], split_pp: bool) -> str:
        lines = []
        for idx, ((row, col), buf) in enumerate(zip(pairs, buf_indices)):
            target = f"out_pp[{buf}] +=" if split_pp and row >= self.x_size + self.u_size and col >= self.x_size + self.u_size else f"out[{buf}] ="
            lines.append(f"    {target} tmp[{idx}];")
        return "\n".join(lines) or "    (void)out;"

    def _render_ode_jac_copy(self, ode_jac: list[tuple[int, int]], lfg_jac: list[tuple[int, int]]) -> str:
        lines = []
        f_offset = 1 if self.lagrange is not None else 0
        lfg_map = {pair: idx for idx, pair in enumerate(lfg_jac)}
        for out_idx, (row, col) in enumerate(ode_jac):
            lfg_idx = lfg_map.get((f_offset + row, col))
            if lfg_idx is None:
                lines.append(f"    dfdx[{out_idx}] = 0.0;")
            else:
                lines.append(f"    dfdx[{out_idx}] = all[{lfg_idx}];")
        return "\n".join(lines) or "    (void)dfdx;"

    def _render_h(self) -> str:
        guard = f"moo_{self.name.upper()}_H"
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
        lines = [
            f"model={self.name}",
            "problem=GDOP",
            f"derivative_strategy={self.derivative_strategy}",
            f"generated_c_bytes={c_path.stat().st_size}",
        ]
        for prefix in ("LFG", "MR"):
            report = emitted.get(f"{prefix}_REPORT", {})
            jac_mode = str(emitted.get(f"{prefix}_JAC_MODE", "direct"))
            hes_mode = str(emitted.get(f"{prefix}_HES_MODE", "direct"))
            lines.extend([f"{prefix.lower()}_jacobian_mode={jac_mode}", f"{prefix.lower()}_hessian_mode={hes_mode}"])
            if isinstance(report, dict):
                lines.extend(f"{prefix.lower()}_{key}={value}" for key, value in sorted(report.items()))
        (out_path / "codegen_report.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
