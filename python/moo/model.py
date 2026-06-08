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

import math
import re
from typing import Iterable

from .common import SolverMixin, SolverSettings


def _clean_name(name: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z_]", "_", name)
    if not cleaned or cleaned[0].isdigit():
        cleaned = f"model_{cleaned}"
    return cleaned


def _num(value: float | int | None) -> str:
    if value is None:
        return "0.0"
    if value == math.inf:
        return "DBL_MAX"
    if value == -math.inf:
        return "-DBL_MAX"
    return repr(float(value))


def _c_bool(value: bool) -> str:
    return "true" if value else "false"


def _arr(values: Iterable[int], ctype: str = "int") -> str:
    vals = list(values)
    if not vals:
        return f"({ctype}[]){{}}"
    return f"({ctype}[]){{{', '.join(map(str, vals))}}}"


def _double_arr(values: Iterable[float]) -> str:
    vals = list(values)
    if not vals:
        return "(f64[]){}"
    return f"(f64[]){{{', '.join(_num(v) for v in vals)}}}"


class BaseModel(SolverMixin):
    def __init__(self, name: str):
        self.name = _clean_name(name)
        self.solver_settings = SolverSettings()
        self.derivative_strategy = "auto"
        self.runtime_parameters: list[tuple[str, float]] = []
        self.derivative_test = False
