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

from pathlib import Path
from time import perf_counter

from moo import nlp_model


def build_model(name: str, n: int = 500):
    model = nlp_model(name)
    x = model.add_variables("x", n, lb=-10.0, ub=10.0, guess=0.1, nominal=1.0)

    model.minimize_sum(n, lambda i: (x[i] - 1.0) ** 2, name="tracking")
    model.minimize_sum(n - 1, lambda i: 0.01 * x[i] * x[i + 1], name="coupling")
    model.add_constraints(n - 1, lambda i: x[i] + x[i + 1], lb=-5.0, ub=5.0, name="band1")
    model.add_constraints(n - 1, lambda i: x[i] * (x[i] + x[i + 1]), lb=-5.0, ub=5.0, name="band2")
    model.add_constraints(n - 1, lambda i: x[i] - x[i + 1] * x[i + 1], lb=-5.0, ub=5.0, name="band3")
    model.solver(tolerance=1e-10, derivative_test=False)

    return model


def print_table(title: str, rows: list[dict]):
    if not rows:
        return

    headers = list(rows[0].keys())

    widths = {
        header: max(len(str(header)), *(len(str(row.get(header, ""))) for row in rows))
        for header in headers
    }

    def separator(char: str = "-") -> str:
        return "+" + "+".join(char * (widths[header] + 2) for header in headers) + "+"

    def format_row(row: dict) -> str:
        return "| " + " | ".join(str(row.get(header, "")).ljust(widths[header]) for header in headers) + " |"

    print(f"\n{title}")
    print(separator("-"))
    print(format_row({header: header for header in headers}))
    print(separator("="))

    for row in rows:
        print(format_row(row))

    print(separator("-"))


def read_codegen_report(path: Path) -> dict:
    report = {}

    for line in path.read_text().splitlines():
        line = line.strip()

        if not line or "=" not in line:
            continue

        key, value = line.split("=", 1)
        report[key.strip()] = value.strip()

    return report


def report_get(report: dict, key: str, default: str = "-") -> str:
    return report.get(key, default)


def apply_codegen(model, strategy):
    if isinstance(strategy, str):
        return model.codegen(strategy)

    return model.codegen(**strategy)


def run_status(result):
    if hasattr(result, "status"):
        return result.status

    if hasattr(result, "returncode"):
        return result.returncode

    if hasattr(result, "process") and hasattr(result.process, "returncode"):
        return result.process.returncode

    return None


def codegen_summary_row(label: str, report: dict, n: int) -> dict:
    return {
        "formulation": label,
        "n": n,
        "strategy": report_get(report, "strategy"),
        "structure": report_get(report, "structure_mode"),
        "local jac": report_get(report, "local_jacobian_mode"),
        "local hess": report_get(report, "local_hessian_mode"),
        "jacobian": report_get(report, "jacobian_mode"),
        "hessian": report_get(report, "hessian_mode"),
        "C bytes": report_get(report, "generated_c_bytes"),
        "Jac nnz": report_get(report, "derivative_jacobian_nnz"),
        "Hess nnz": report_get(report, "derivative_hessian_nnz"),
        "graph": report_get(report, "derivative_graph_source"),
    }


def local_block_row(label: str, report: dict, n: int, kind: str, index: int) -> dict:
    prefix = f"derivative_local_{kind}_{index}"

    return {
        "formulation": label,
        "n": n,
        "block": f"{kind}_{index}",
        "repetitions": report_get(report, f"{prefix}_repetitions"),
        "input": report_get(report, f"{prefix}_input_size"),
        "output": report_get(report, f"{prefix}_output_size"),
        "Jac mode": report_get(report, f"{prefix}_jacobian_mode"),
        "Jac nnz": report_get(report, f"{prefix}_jacobian_nnz"),
        "Jac colors": report_get(report, f"{prefix}_jacobian_colors"),
        "Hess mode": report_get(report, f"{prefix}_hessian_mode"),
        "Hess nnz": report_get(report, f"{prefix}_hessian_nnz"),
        "Hess colors": report_get(report, f"{prefix}_hessian_colors"),
    }


def mapped_summary_row(label: str, report: dict, n: int) -> dict:
    return {
        "formulation": label,
        "n": n,
        "constraint blocks": report_get(report, "derivative_mapped_constraint_blocks"),
        "constraint outputs": report_get(report, "derivative_mapped_constraint_outputs"),
        "objective blocks": report_get(report, "derivative_mapped_objective_blocks"),
        "mapped repetitions": report_get(report, "derivative_mapped_repetitions"),
    }


if __name__ == "__main__":
    n = 500

    cases = [
        {
            "label": "auto",
            "model_name": "nlp_sparse_auto",
            "out_dir": Path("build/moo/nlp_sparse_auto"),
            "strategy": "auto",
        },
        {
            "label": "loop",
            "model_name": "nlp_sparse_loop",
            "out_dir": Path("build/moo/nlp_sparse_loop"),
            "strategy": "loop",
        },
        {
            "label": "direct",
            "model_name": "nlp_sparse_direct",
            "out_dir": Path("build/moo/nlp_sparse_direct"),
            "strategy": "direct",
        },
        {
            "label": "colored",
            "model_name": "nlp_sparse_colored",
            "out_dir": Path("build/moo/nlp_sparse_colored"),
            "strategy": "colored",
        },
        {
            "label": "basis",
            "model_name": "nlp_sparse_basis",
            "out_dir": Path("build/moo/nlp_sparse_basis"),
            "strategy": "basis",
        },
        {
            "label": "loop-direct",
            "model_name": "nlp_sparse_loop_direct",
            "out_dir": Path("build/moo/nlp_sparse_loop_direct"),
            "strategy": "loop-direct",
        },
        {
            "label": "loop-colored",
            "model_name": "nlp_sparse_loop_colored",
            "out_dir": Path("build/moo/nlp_sparse_loop_colored"),
            "strategy": "loop-colored",
        },
        {
            "label": "loop-colored-kwargs",
            "model_name": "nlp_sparse_loop_colored_kwargs",
            "out_dir": Path("build/moo/nlp_sparse_loop_colored_kwargs"),
            "strategy": {
                "structure": "loop",
                "local": "colored",
            },
        },
    ]

    generated_cases = []
    timing_rows = []

    for case in cases:
        t0 = perf_counter()
        model = build_model(case["model_name"], n=n)
        build_elapsed = perf_counter() - t0

        t0 = perf_counter()
        generated = apply_codegen(model, case["strategy"])
        c_path, _ = generated.generate(case["out_dir"])
        codegen_elapsed = perf_counter() - t0

        generated_cases.append(
            {
                **case,
                "generated": generated,
                "c_path": c_path,
                "build_time": build_elapsed,
                "codegen_time": codegen_elapsed,
            }
        )

        print(
            f"{case['label']:20s} build_model: {build_elapsed:.6f} s, "
            f"codegen: {codegen_elapsed:.6f} s",
            flush=True,
        )

    for case in generated_cases:
        t0 = perf_counter()
        result = case["generated"].optimize(case["out_dir"], solver="Ipopt")
        optimize_elapsed = perf_counter() - t0

        timing_rows.append(
            {
                "formulation": case["label"],
                "n": n,
                "build_model [s]": f"{case['build_time']:.6f}",
                "codegen [s]": f"{case['codegen_time']:.6f}",
                "optimization [s]": f"{optimize_elapsed:.6f}",
                "C bytes": case["c_path"].stat().st_size,
                "status": run_status(result),
            }
        )

        print(
            f"{case['label']:20s} optimization time (n={n}): {optimize_elapsed:.6f} s",
            flush=True,
        )

    reports = {
        case["label"]: read_codegen_report(case["out_dir"] / "codegen_report.txt")
        for case in generated_cases
    }

    codegen_summary_rows = [
        codegen_summary_row(case["label"], reports[case["label"]], n)
        for case in generated_cases
    ]

    local_derivative_rows = []

    for case in generated_cases:
        label = case["label"]
        report = reports[label]

        local_derivative_rows.extend(
            [
                local_block_row(label, report, n, "constraint", 0),
                local_block_row(label, report, n, "objective", 0),
                local_block_row(label, report, n, "objective", 1),
            ]
        )

    mapped_summary_rows = [
        mapped_summary_row(case["label"], reports[case["label"]], n)
        for case in generated_cases
    ]

    print_table("Timings", timing_rows)

    print()
    print_table("Codegen summary", codegen_summary_rows)

    print()
    print_table("Local derivative blocks", local_derivative_rows)

    print()
    print_table("Mapped derivative summary", mapped_summary_rows)
