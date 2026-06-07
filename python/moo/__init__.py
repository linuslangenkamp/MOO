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

from pkgutil import extend_path
from pathlib import Path

__path__ = extend_path(__path__, __name__)
_build_ext_dir = Path(__file__).resolve().parents[2] / "build" / "python" / "moo"
if _build_ext_dir.exists():
    __path__.append(str(_build_ext_dir))

from .model import (
    Expr,
    GDOPModel,
    Model,
    cos,
    exp,
    log,
    pow_const,
    sin,
    sum_expr,
    tan,
)

from .init_model import InitModel
from .nlp_model import NLPModel
from .nlp_model import dot, matrix, vec, vector
from .radauIIA import RADAU, RadauConstants, radauIIA
from .results import OptimizationRun, ResultSet, ResultTable, read_results


def gdop_model(name: str) -> GDOPModel:
    return GDOPModel(name)


def init_model(name: str) -> InitModel:
    return InitModel(name)


def nlp_model(name: str) -> NLPModel:
    return NLPModel(name)

__version__ = "0.1.0"

__all__ = [
    "Expr",
    "Model",
    "GDOPModel",
    "InitModel",
    "NLPModel",
    "gdop_model",
    "init_model",
    "nlp_model",
    "OptimizationRun",
    "ResultSet",
    "ResultTable",
    "read_results",
    "vec",
    "vector",
    "matrix",
    "dot",
    "RADAU",
    "RadauConstants",
    "radauIIA",
    "sin",
    "sum_expr",
    "cos",
    "tan",
    "exp",
    "log",
    "pow_const",
]
