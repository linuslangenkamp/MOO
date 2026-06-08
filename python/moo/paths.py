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

import os
from pathlib import Path


def package_root() -> Path:
    return Path(__file__).resolve().parent


def source_root() -> Path:
    env = os.environ.get("MOO_REPO_ROOT")
    if env:
        return Path(env).expanduser().resolve()
    root = package_root().parents[1]
    if (root / "src").is_dir():
        return root
    return root


def include_dir() -> Path:
    env = os.environ.get("MOO_INCLUDE_DIR")
    if env:
        return Path(env).expanduser().resolve()
    root = source_root()
    src = root / "src"
    if src.is_dir():
        return src
    installed = package_root() / "include"
    return installed.resolve()


def library_dir(build_dir: str | os.PathLike[str] = "build") -> Path:
    env = os.environ.get("MOO_LIBRARY_DIR")
    if env:
        return Path(env).expanduser().resolve()
    build = Path(build_dir).expanduser()
    if build.is_absolute():
        return build.resolve()
    return (source_root() / build).resolve()


def library_path(build_dir: str | os.PathLike[str] = "build") -> Path:
    lib_dir = library_dir(build_dir)
    for name in ("libmoo.so", "libmoo.dylib", "moo.dll"):
        candidate = lib_dir / name
        if candidate.exists():
            return candidate
    return lib_dir / "libmoo.so"


def data_dir() -> Path:
    env = os.environ.get("MOO_DATA_DIR")
    if env:
        return Path(env).expanduser().resolve()
    root = source_root()
    data = root / "data"
    if data.is_dir():
        return data
    return (package_root() / "data").resolve()
