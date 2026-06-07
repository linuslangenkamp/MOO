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
from pathlib import Path


@dataclass(frozen=True)
class RadauConstants:
    """Radau IIA constants for one collocation stage count."""

    stages: int
    B: tuple[float, ...]
    C: tuple[float, ...]
    C0: tuple[float, ...]
    D: tuple[tuple[float, ...], ...]
    W: tuple[float, ...]
    W0: tuple[float, ...]

    @property
    def b(self) -> tuple[float, ...]:
        return self.B

    @property
    def c(self) -> tuple[float, ...]:
        return self.C

    @property
    def c0(self) -> tuple[float, ...]:
        return self.C0

    @property
    def d(self) -> tuple[tuple[float, ...], ...]:
        return self.D

    @property
    def w(self) -> tuple[float, ...]:
        return self.W

    @property
    def w0(self) -> tuple[float, ...]:
        return self.W0

    def as_dict(self) -> dict[str, object]:
        return {
            "stages": self.stages,
            "B": self.B,
            "C": self.C,
            "C0": self.C0,
            "D": self.D,
            "W": self.W,
            "W0": self.W0,
        }


def _data_dir() -> Path:
    return Path(__file__).resolve().parents[2] / "data/fLGR"


def _numbers_from_line(line: str) -> tuple[float, ...]:
    cleaned = line.strip().replace("{", "").replace("}", "")
    cleaned = cleaned.replace(";", "")
    cleaned = cleaned.rstrip(",")
    if not cleaned:
        return ()
    return tuple(float(part.strip()) for part in cleaned.split(",") if part.strip())


def _read_vector_rows(filename: str) -> list[tuple[float, ...]]:
    path = _data_dir() / filename
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        numbers = _numbers_from_line(line)
        if numbers:
            rows.append(numbers)
    return rows


def _read_d_matrices(max_stages: int) -> dict[int, tuple[tuple[float, ...], ...]]:
    rows = _read_vector_rows("radauConstantsD.data")
    matrices = {}
    cursor = 0
    for stages in range(max_stages + 1):
        width = stages + 1
        block = rows[cursor:cursor + width]
        if len(block) != width or any(len(row) != width for row in block):
            raise ValueError(f"Invalid Radau D block for {stages} stages")
        matrices[stages] = tuple(tuple(row) for row in block)
        cursor += width
    return matrices


def _load_radau(max_stages: int = 10) -> dict[int, RadauConstants]:
    b_rows = _read_vector_rows("radauConstantsB.data")
    c_rows = _read_vector_rows("radauConstantsC.data")
    c0_rows = _read_vector_rows("radauConstantsC0.data")
    w_rows = _read_vector_rows("radauConstantsW.data")
    w0_rows = _read_vector_rows("radauConstantsW0.data")
    d_matrices = _read_d_matrices(max_stages)

    constants = {}
    for stages in range(1, max_stages + 1):
        constants[stages] = RadauConstants(
            stages=stages,
            B=b_rows[stages - 1],
            C=c_rows[stages - 1],
            C0=c0_rows[stages],
            D=d_matrices[stages],
            W=w_rows[stages - 1],
            W0=w0_rows[stages],
        )
    return constants


RADAU: dict[int, RadauConstants] = _load_radau(25)


def radauIIA(stages: int) -> RadauConstants:
    """Return Radau IIA constants for ``stages`` collocation points."""

    try:
        return RADAU[stages]
    except KeyError as exc:
        available = f"{min(RADAU)}..{max(RADAU)}"
        raise ValueError(f"Radau constants are available for stages {available}") from exc
