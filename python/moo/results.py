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

import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


class ResultTable:
    def __init__(self, name: str, path: Path, header: list[str], rows: list[dict[str, str]]):
        self.name = name
        self.path = path
        self.header = header
        self.rows = rows

    def column(self, name: str, numeric: bool = True) -> list[float] | list[str]:
        values = [row[name] for row in self.rows]
        if numeric:
            return [float(value) for value in values]
        return values

    def dataframe(self):
        try:
            import pandas as pd
        except ImportError as exc:
            raise RuntimeError("pandas is not installed; use table(...).rows or install pandas") from exc
        return pd.read_csv(self.path)


class ResultSet:
    def __init__(self, directory: Path, tables: dict[str, ResultTable]):
        self.directory = directory
        self.tables = tables

    def table(self, name: str) -> ResultTable:
        key = name.removesuffix(".csv")
        if key not in self.tables:
            available = ", ".join(sorted(self.tables))
            raise KeyError(f"Unknown result table {name!r}; available: {available}")
        return self.tables[key]

    def dataframe(self, name: str):
        return self.table(name).dataframe()


def _line_marker_style(
    plottype: str,
    linestyle: str,
    nodestyle: str,
) -> tuple[str, str]:
    if plottype not in {"line", "nodes", "both"}:
        raise ValueError("plottype must be one of 'line', 'nodes', or 'both'")

    if plottype == "line":
        return linestyle, ""

    if plottype == "nodes":
        return "", nodestyle

    return linestyle, nodestyle


class PlotAccessor:
    def __init__(self, result: "BaseResult"):
        self._result = result

    def _finish(self, fig, show: bool = False, save: str | Path | None = None):
        if save is not None:
            path = Path(save)
            path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(path, bbox_inches="tight")
        if show:
            import matplotlib.pyplot as plt
            plt.show()
        return fig

    def states(
        self,
        show: bool = False,
        save: str | Path | None = None,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        return self._finish(
            self._result._plot_kind(
                "states",
                plottype=plottype,
                linestyle=linestyle,
                nodestyle=nodestyle,
            ),
            show,
            save,
        )

    def controls(
        self,
        show: bool = False,
        save: str | Path | None = None,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        return self._finish(
            self._result._plot_kind(
                "controls",
                plottype=plottype,
                linestyle=linestyle,
                nodestyle=nodestyle,
            ),
            show,
            save,
        )

    def costates(
        self,
        show: bool = False,
        save: str | Path | None = None,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        return self._finish(
            self._result._plot_kind(
                "costates",
                plottype=plottype,
                linestyle=linestyle,
                nodestyle=nodestyle,
            ),
            show,
            save,
        )

    def variables(self, show: bool = False, save: str | Path | None = None):
        return self._finish(self._result._plot_kind("variables"), show, save)

    def constraints(self, show: bool = False, save: str | Path | None = None):
        return self._finish(self._result._plot_kind("constraints"), show, save)

    def all(
        self,
        show: bool = False,
        save: str | Path | None = None,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        return self._finish(
            self._result._plot_all(
                plottype=plottype,
                linestyle=linestyle,
                nodestyle=nodestyle,
            ),
            show,
            save,
        )


class BaseResult:
    def __init__(self, raw: ResultSet):
        self.raw = raw
        self.plot = PlotAccessor(self)

    def table(self, name: str) -> ResultTable:
        return self.raw.table(name)

    def dataframe(self, name: str):
        return self.raw.dataframe(name)

    def _plot_kind(
        self,
        kind: str,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        raise NotImplementedError(f"{type(self).__name__} does not support plot.{kind}()")

    def _plot_all(
        self,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        raise NotImplementedError(f"{type(self).__name__} does not support plot.all()")


class GDOPResult(BaseResult):
    def __init__(self, raw: ResultSet, state_names: list[str], control_names: list[str]):
        super().__init__(raw)
        self.state_names = state_names
        self.control_names = control_names
        self.time = self._column("primals_optimal_solution", "time")
        self.states = {
            name: self._column("primals_optimal_solution", f"x_{idx}")
            for idx, name in enumerate(state_names)
            if f"x_{idx}" in self.raw.table("primals_optimal_solution").header
        }
        self.controls = {
            name: self._column("primals_optimal_solution", f"u_{idx}")
            for idx, name in enumerate(control_names)
            if f"u_{idx}" in self.raw.table("primals_optimal_solution").header
        }
        self.costates = self._prefixed_columns("costates_optimal_solution")
        self.lower_costates = self._prefixed_columns("lower_costates_optimal_solution")
        self.upper_costates = self._prefixed_columns("upper_costates_optimal_solution")

    def _column(self, table: str, col: str) -> list[float]:
        return self.raw.table(table).column(col)

    def _prefixed_columns(self, table: str) -> dict[str, list[float]]:
        try:
            tbl = self.raw.table(table)
        except KeyError:
            return {}
        return {col: tbl.column(col) for col in tbl.header if col != "time"}

    def _plot_kind(
        self,
        kind: str,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        import matplotlib.pyplot as plt

        line, marker = _line_marker_style(plottype, linestyle, nodestyle)

        fig, ax = plt.subplots()
        data = getattr(self, kind)

        for name, values in data.items():
            ax.plot(
                self.time,
                values,
                linestyle=line,
                marker=marker,
                label=name,
            )

        ax.set_xlabel("time")
        ax.set_ylabel(kind)
        ax.grid(True)
        ax.legend()
        return fig

    def _plot_all(
        self,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        import matplotlib.pyplot as plt

        line, marker = _line_marker_style(plottype, linestyle, nodestyle)

        fig, axes = plt.subplots(2, 1, sharex=True)

        for name, values in self.states.items():
            axes[0].plot(
                self.time,
                values,
                linestyle=line,
                marker=marker,
                label=name,
            )

        for name, values in self.controls.items():
            axes[1].plot(
                self.time,
                values,
                linestyle=line,
                marker=marker,
                label=name,
            )

        axes[0].set_ylabel("states")
        axes[1].set_ylabel("controls")
        axes[1].set_xlabel("time")

        for ax in axes:
            ax.grid(True)
            ax.legend()

        return fig


class InitResult(BaseResult):
    def __init__(self, raw: ResultSet, variable_names: list[str], parameter_names: list[str]):
        super().__init__(raw)
        row = raw.table("init_optimal_solution").rows[0] if raw.table("init_optimal_solution").rows else {}
        self.objective = float(row.get("objective", 0.0))
        self.variables = {
            name: float(row[f"y_{idx}"])
            for idx, name in enumerate(variable_names)
            if f"y_{idx}" in row
        }
        self.parameters = {
            name: float(row[f"p_{idx}"])
            for idx, name in enumerate(parameter_names)
            if f"p_{idx}" in row
        }
        self.deltas = {
            name: float(row[f"dp_{idx}"])
            for idx, name in enumerate(parameter_names)
            if f"dp_{idx}" in row
        }
        self.f = {key: float(value) for key, value in row.items() if key.startswith("f_")}
        self.g = {key: float(value) for key, value in row.items() if key.startswith("g_")}
        self.lambdas = {key: float(value) for key, value in row.items() if key.startswith("lambda_")}
        self.bound_duals = {
            key: float(value)
            for key, value in row.items()
            if key.startswith("z_lb_") or key.startswith("z_ub_")
        }

    def _plot_kind(
        self,
        kind: str,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        import matplotlib.pyplot as plt

        data = getattr(self, kind)
        fig, ax = plt.subplots()
        ax.bar(list(data), list(data.values()))
        ax.set_ylabel(kind)
        ax.grid(True, axis="y")
        return fig

    def _plot_all(
        self,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        return self._plot_kind("variables")


class NLPResult(BaseResult):
    def __init__(self, raw: ResultSet, variable_names: list[str], constraint_names: list[str]):
        super().__init__(raw)
        row = raw.table("nlp_optimal_solution").rows[0] if raw.table("nlp_optimal_solution").rows else {}
        self.objective = float(row.get("objective", 0.0))
        self.variables = {
            name: float(row[f"x_{idx}"])
            for idx, name in enumerate(variable_names)
            if f"x_{idx}" in row
        }
        self.constraints = {
            name: float(row[f"g_{idx}"])
            for idx, name in enumerate(constraint_names)
            if f"g_{idx}" in row
        }
        self.lambdas = {key: float(value) for key, value in row.items() if key.startswith("lambda_")}
        self.bound_duals = {
            key: float(value)
            for key, value in row.items()
            if key.startswith("z_lb_") or key.startswith("z_ub_")
        }

    def _plot_kind(
        self,
        kind: str,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        import matplotlib.pyplot as plt

        data = getattr(self, kind)
        fig, ax = plt.subplots()
        ax.bar(list(data), list(data.values()))
        ax.set_ylabel(kind)
        ax.grid(True, axis="y")
        return fig

    def _plot_all(
        self,
        plottype: str = "both",
        linestyle: str = "-",
        nodestyle: str = "",
    ):
        return self._plot_kind("variables")


@dataclass
class OptimizationRun:
    process: subprocess.CompletedProcess[str]
    results: ResultSet
    result: BaseResult | None = None

    def read_results(self) -> ResultSet:
        self.results = read_results(self.results.directory)
        return self.results

    def table(self, name: str) -> ResultTable:
        return self.results.table(name)

    def dataframe(self, name: str):
        return self.results.dataframe(name)

    @property
    def returncode(self) -> int:
        return self.process.returncode


def read_results(directory: str | Path) -> ResultSet:
    base = Path(directory)
    tables: dict[str, ResultTable] = {}
    for path in sorted(base.glob("*.csv")):
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            header = list(reader.fieldnames or [])
        tables[path.stem] = ResultTable(path.stem, path, header, rows)
    return ResultSet(base, tables)
