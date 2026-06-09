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

from dataclasses import dataclass

from . import ad
from .ad_codegen import EmittedFunction, emit_function, emit_value_function


@dataclass(frozen=True)
class InputGroup:
    name: str
    size: int


@dataclass
class EmittedLocalGraphFunction:
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


def with_derivative_modes(emitted: EmittedFunction, strategy: str) -> EmittedLocalGraphFunction:
    return EmittedLocalGraphFunction(
        value=emitted.value,
        jvp=emitted.jvp,
        hvp=emitted.hvp,
        jac=emitted.jac,
        hes=emitted.hes,
        jac_sparsity=emitted.jac_sparsity,
        hes_sparsity=emitted.hes_sparsity,
        jac_colors=emitted.jac_colors,
        hes_colors=emitted.hes_colors,
        jac_mode=ad.derivative_callback_mode(strategy, emitted.jac_sparsity, emitted.jac_colors),
        hes_mode=ad.derivative_callback_mode(strategy, emitted.hes_sparsity, emitted.hes_colors),
        report=emitted.report,
    )


@dataclass
class LocalGraphFunction:
    name: str
    input: InputGroup
    outputs: list[object]
    params: list[InputGroup]

    def emit(self) -> EmittedFunction:
        return emit_function(
            self.input.name,
            self.input.size,
            self.outputs,
            [(param.name, param.size) for param in self.params],
            f"{self.name}_value",
            f"{self.name}_jvp",
            f"{self.name}_hvp",
        )

    def emit_with_modes(self, strategy: str) -> EmittedLocalGraphFunction:
        return with_derivative_modes(self.emit(), strategy)


@dataclass
class LocalValueFunction:
    name: str
    input: InputGroup
    outputs: list[object]
    params: list[InputGroup]

    def emit(self) -> tuple[str, list[tuple[int, int]]]:
        return emit_value_function(
            self.input.name,
            self.input.size,
            self.outputs,
            [(param.name, param.size) for param in self.params],
            f"{self.name}_value",
        )
