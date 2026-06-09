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


def static_int_array(name: str, values: list[int], indent: str) -> str:
    return f"{indent}static const int {name}[{max(len(values), 1)}] = {{ {', '.join(map(str, values)) or '0'} }};"


def color_metadata(pairs: list[tuple[int, int]], colors: list[int]) -> tuple[list[int], list[int], list[int], list[int], list[int], list[int]]:
    by_color: dict[int, list[tuple[int, int, int]]] = {}
    for idx, (row, col) in enumerate(pairs):
        color = colors[col] if 0 <= col < len(colors) else col
        by_color.setdefault(color, []).append((idx, row, col))

    color_cols: list[int] = []
    color_offsets = [0]
    scatter_idx: list[int] = []
    scatter_row: list[int] = []
    scatter_col: list[int] = []
    scatter_offsets = [0]
    for color in sorted(by_color):
        entries = by_color[color]
        cols = sorted({col for _, _, col in entries})
        color_cols.extend(cols)
        color_offsets.append(len(color_cols))
        for idx, row, col in entries:
            scatter_idx.append(idx)
            scatter_row.append(row)
            scatter_col.append(col)
        scatter_offsets.append(len(scatter_idx))
    return color_offsets, color_cols, scatter_offsets, scatter_idx, scatter_row, scatter_col


def render_local_colored_jac_lines(
    fn: str,
    pairs: list[tuple[int, int]],
    colors: list[int],
    input_size: int,
    output_size: int,
    out_name: str,
    rep_name: str,
    accumulate: bool,
    indent: str,
    base_buf: int = 0,
    jbuf_name: str | None = None,
) -> list[str]:
    color_offsets, color_cols, scatter_offsets, scatter_idx, scatter_row, _ = color_metadata(pairs, colors)
    lines = [
        f"{indent}f64 v[{max(input_size, 1)}] = {{0}};",
        f"{indent}f64 tmp_color[{max(output_size, 1)}] = {{0}};",
        static_int_array("color_offsets", color_offsets, indent),
        static_int_array("color_cols", color_cols, indent),
        static_int_array("scatter_offsets", scatter_offsets, indent),
        static_int_array("scatter_idx", scatter_idx, indent),
        static_int_array("scatter_row", scatter_row, indent),
        f"{indent}for (int color = 0; color < {max(len(color_offsets) - 1, 0)}; ++color) {{",
        f"{indent}    for (int k = color_offsets[color]; k < color_offsets[color + 1]; ++k) {{ v[color_cols[k]] = 1.0; }}",
        f"{indent}    {fn}(xl, rp, v, tmp_color);",
        f"{indent}    for (int k = scatter_offsets[color]; k < scatter_offsets[color + 1]; ++k) {{",
    ]
    if accumulate and jbuf_name is not None:
        lines.append(f"{indent}        {out_name}[{jbuf_name}[{rep_name} * {len(pairs)} + scatter_idx[k]]] += tmp_color[scatter_row[k]];")
    elif accumulate:
        lines.append(f"{indent}        {out_name}[scatter_idx[k]] += tmp_color[scatter_row[k]];")
    else:
        lines.append(f"{indent}        {out_name}[{base_buf} + {rep_name} * {len(pairs)} + scatter_idx[k]] = tmp_color[scatter_row[k]];")
    lines.extend([
        f"{indent}    }}",
        f"{indent}    for (int k = color_offsets[color]; k < color_offsets[color + 1]; ++k) {{ v[color_cols[k]] = 0.0; }}",
        f"{indent}}}",
    ])
    return lines


def render_local_direct_objective_jac_lines(fn: str, pairs: list[tuple[int, int]], jbuf_name: str, indent: str, out_name: str = "out") -> list[str]:
    if not pairs:
        return []
    lines = [
        f"{indent}f64 tmp[{max(len(pairs), 1)}];",
        f"{indent}{fn}_sparse(xl, rp, tmp);",
    ]
    if all(row == 0 for row, _ in pairs):
        lines.extend([
            f"{indent}for (int local_buf = 0; local_buf < {len(pairs)}; ++local_buf) {{",
            f"{indent}    {out_name}[{jbuf_name}[rep * {len(pairs)} + local_buf]] += tmp[local_buf];",
            f"{indent}}}",
        ])
    else:
        for local_buf, (row, _) in enumerate(pairs):
            if row == 0:
                lines.append(f"{indent}{out_name}[{jbuf_name}[rep * {len(pairs)} + {local_buf}]] += tmp[{local_buf}];")
    return lines


def render_local_direct_constraint_jac_lines(fn: str, pairs: list[tuple[int, int]], base_buf: int, indent: str, out_name: str = "out") -> list[str]:
    if not pairs:
        return []
    return [
        f"{indent}f64 tmp[{max(len(pairs), 1)}];",
        f"{indent}{fn}_sparse(xl, rp, tmp);",
        f"{indent}for (int local_buf = 0; local_buf < {len(pairs)}; ++local_buf) {{",
        f"{indent}    {out_name}[{base_buf} + rep * {len(pairs)} + local_buf] = tmp[local_buf];",
        f"{indent}}}",
    ]


def render_local_direct_hes_lines(fn: str, pairs: list[tuple[int, int]], hbuf_name: str, tmp_name: str, indent: str) -> list[str]:
    if not pairs:
        return []
    return [
        f"{indent}f64 {tmp_name}[{max(len(pairs), 1)}];",
        f"{indent}{fn}_sparse(xl, seed, rp, {tmp_name});",
        f"{indent}for (int local_buf = 0; local_buf < {len(pairs)}; ++local_buf) {{",
        f"{indent}    out[{hbuf_name}[rep * {len(pairs)} + local_buf]] += {tmp_name}[local_buf];",
        f"{indent}}}",
    ]


def render_local_colored_hes_lines(
    fn: str,
    pairs: list[tuple[int, int]],
    colors: list[int],
    input_size: int,
    hbuf: list[int],
    indent: str,
    hbuf_name: str | None = None,
) -> list[str]:
    if not pairs:
        return []
    color_offsets, color_cols, scatter_offsets, scatter_idx, scatter_row, _ = color_metadata(pairs, colors)
    hbuf_name = hbuf_name or "h_buf_for_local"
    lines = [
        f"{indent}{{",
        f"{indent}f64 v_h[{max(input_size, 1)}] = {{0}};",
        f"{indent}f64 tmp_h[{max(input_size, 1)}] = {{0}};",
        static_int_array("h_color_offsets", color_offsets, indent),
        static_int_array("h_color_cols", color_cols, indent),
        static_int_array("h_scatter_offsets", scatter_offsets, indent),
        static_int_array("h_scatter_idx", scatter_idx, indent),
        static_int_array("h_scatter_row", scatter_row, indent),
        f"{indent}for (int color = 0; color < {max(len(color_offsets) - 1, 0)}; ++color) {{",
        f"{indent}    for (int k = h_color_offsets[color]; k < h_color_offsets[color + 1]; ++k) {{ v_h[h_color_cols[k]] = 1.0; }}",
        f"{indent}    {fn}(xl, seed, rp, v_h, tmp_h);",
    ]
    if hbuf_name == "h_buf_for_local":
        lines.append(static_int_array("h_buf_for_local", hbuf, indent))
    lines.extend([
        f"{indent}    for (int k = h_scatter_offsets[color]; k < h_scatter_offsets[color + 1]; ++k) {{ out[{hbuf_name}[rep * {len(pairs)} + h_scatter_idx[k]]] += tmp_h[h_scatter_row[k]]; }}",
        f"{indent}    for (int k = h_color_offsets[color]; k < h_color_offsets[color + 1]; ++k) {{ v_h[h_color_cols[k]] = 0.0; }}",
        f"{indent}}}",
        f"{indent}}}",
    ])
    return lines
