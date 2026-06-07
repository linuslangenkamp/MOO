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

from dataclasses import dataclass
from typing import Iterable


@dataclass
class SolverSettings:
    backend: str = "Ipopt"
    tolerance: float = 1e-10
    iterations: int = 5000
    derivative_test: bool = False
    linear_solver: str = "MUMPS"
    hessian: str = "Exact"
    jacobian: str = "Exact"
    gradient: str = "Exact"
    uno_preset: str = "ipopt"
    cpu_time: float = 3600.0
    warm_start: bool = False
    qp: bool = False

    def cli_args(self, solver: str | None = None, solver_args: list[str] | None = None) -> list[str]:
        args = list(solver_args or [])
        backend = solver or self.backend
        normalized = backend.strip().lower()
        if normalized == "ipopt":
            backend_value = "Ipopt"
        elif normalized == "uno":
            backend_value = "Uno"
        else:
            raise ValueError(f"Unknown solver {backend!r}; expected 'Ipopt' or 'Uno'")

        defaults = {
            "NLPSolver": backend_value,
            "Tolerance": self.tolerance,
            "Iterations": self.iterations,
            "CPUTime": self.cpu_time,
            "LinearSolver": self.linear_solver,
            "Hessian": self.hessian,
            "Jacobian": self.jacobian,
            "Gradient": self.gradient,
            "IpoptDerivativeTest": self.derivative_test,
            "WarmStart": self.warm_start,
            "QP": self.qp,
            "UnoPreset": self.uno_preset,
        }

        existing = {
            arg[2:].split("=", 1)[0]
            for arg in args
            if arg.startswith("--") and "=" in arg
        }
        for key, value in defaults.items():
            if key in existing:
                continue
            if isinstance(value, bool):
                text = "true" if value else "false"
            else:
                text = str(value)
            args.insert(0, f"--{key}={text}")
        return args


class SolverMixin:
    def codegen(self, strategy: str = "auto", structure: str | None = None, local: str | None = None):
        structure_mode, local_mode = parse_codegen_strategy(strategy, structure, local)
        self.codegen_strategy = strategy.strip().lower() if isinstance(strategy, str) else "auto"
        self.codegen_structure_strategy = structure_mode
        self.codegen_local_strategy = local_mode
        return self

    def solver(
        self,
        backend: str = "Ipopt",
        tolerance: float = 1e-10,
        iterations: int = 5000,
        derivative_test: bool = False,
        linear_solver: str = "MUMPS",
        hessian: str = "Exact",
        jacobian: str = "Exact",
        gradient: str = "Exact",
        uno_preset: str = "ipopt",
        cpu_time: float = 3600.0,
        warm_start: bool = False,
        qp: bool = False,
    ):
        self.solver_settings = SolverSettings(
            backend=backend,
            tolerance=tolerance,
            iterations=iterations,
            derivative_test=derivative_test,
            linear_solver=linear_solver,
            hessian=hessian,
            jacobian=jacobian,
            gradient=gradient,
            uno_preset=uno_preset,
            cpu_time=cpu_time,
            warm_start=warm_start,
            qp=qp,
        )
        if hasattr(self, "derivative_test"):
            self.derivative_test = derivative_test
        return self


def derivative_mode(strategy: str, pairs: Iterable[tuple[int, int]], colors: list[int]) -> str:
    sparse_pairs = list(pairs)
    if not sparse_pairs:
        return "direct"
    if strategy == "basis":
        return "basis"
    if strategy == "direct":
        return "direct"
    if strategy == "colored":
        return "colored"

    used_colors = {
        colors[col]
        for _, col in sparse_pairs
        if 0 <= col < len(colors) and colors[col] >= 0
    }
    if len(sparse_pairs) >= 256 and 0 < 4 * len(used_colors) < len(sparse_pairs):
        return "colored"
    return "direct"


def parse_codegen_strategy(strategy: str = "auto", structure: str | None = None, local: str | None = None) -> tuple[str, str]:
    normalized = strategy.strip().lower()
    aliases = {
        "auto": ("auto", "auto"),
        "loop": ("loop", "auto"),
        "direct": ("auto", "direct"),
        "colored": ("auto", "colored"),
        "basis": ("auto", "basis"),
        "loop-auto": ("loop", "auto"),
        "loop-direct": ("loop", "direct"),
        "loop-colored": ("loop", "colored"),
        "loop-basis": ("loop", "basis"),
        "scalar-auto": ("scalar", "auto"),
        "scalar-direct": ("scalar", "direct"),
        "scalar-colored": ("scalar", "colored"),
        "scalar-basis": ("scalar", "basis"),
    }
    if normalized not in aliases:
        raise ValueError(f"Unknown codegen strategy {strategy!r}; expected one of {sorted(aliases)}")
    structure_mode, local_mode = aliases[normalized]
    if structure is not None:
        structure_mode = structure.strip().lower()
    if local is not None:
        local_mode = local.strip().lower()
    if structure_mode not in {"auto", "loop", "scalar"}:
        raise ValueError("codegen structure must be one of 'auto', 'loop', or 'scalar'")
    if local_mode not in {"auto", "direct", "colored", "basis"}:
        raise ValueError("codegen local must be one of 'auto', 'direct', 'colored', or 'basis'")
    return structure_mode, local_mode
