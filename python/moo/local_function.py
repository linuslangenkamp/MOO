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

from .ad_codegen import EmittedFunction, emit_function, emit_value_function


@dataclass(frozen=True)
class InputGroup:
    name: str
    size: int


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
