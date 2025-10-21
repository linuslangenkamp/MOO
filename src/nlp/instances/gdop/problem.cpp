// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2025 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "problem.h"

namespace GDOP {

ProblemConstants::ProblemConstants(bool has_mayer,
                                   bool has_lagrange,
                                   FixedVector<Bounds>&& x_bounds,
                                   FixedVector<Bounds>&& u_bounds,
                                   FixedVector<Bounds>&& p_bounds,
                                   std::array<Bounds, 2> T_bounds,
                                   FixedVector<std::optional<f64>> x0_fixed,
                                   FixedVector<std::optional<f64>> xf_fixed,
                                   std::array<std::optional<f64>, 2> T_fixed,
                                   FixedVector<Bounds>&& r_bounds,
                                   FixedVector<Bounds>&& g_bounds,
                                   Mesh& mesh)
: x_size(static_cast<int>(x_bounds.size())),
  u_size(static_cast<int>(u_bounds.size())),
  p_size(static_cast<int>(p_bounds.size())),
  xu_size(x_size + u_size),
  has_mayer(has_mayer),
  has_lagrange(has_lagrange),
  f_size(x_size),
  g_size(static_cast<int>(g_bounds.size())),
  r_size(static_cast<int>(r_bounds.size())),
  fg_size(f_size + g_size),
  x_bounds(std::move(x_bounds)),
  u_bounds(std::move(u_bounds)),
  p_bounds(std::move(p_bounds)),
  T_bounds(std::move(T_bounds)),
  x0_fixed(std::move(x0_fixed)),
  xf_fixed(std::move(xf_fixed)),
  T_fixed(std::move(T_fixed)),
  free_time(!(T_fixed[0] && T_fixed[1])),
  r_bounds(std::move(r_bounds)),
  g_bounds(std::move(g_bounds)),
  mesh(mesh.shared_from_this()) {}

int FullSweepLayout::compute_jac_nnz() {
    int sum = 0;

    if (L) sum += L->jac.nnz();

    for (auto const& fn : f) {
        sum += fn.jac.nnz();
    }

    for (auto const& fn : g) {
        sum += fn.jac.nnz();
    }

    return sum;
}

int BoundarySweepLayout::compute_jac_nnz() {
    int sum = 0;

    if (M) sum += M->jac.nnz();

    for (auto const& fn : r) {
        sum += fn.jac.nnz();
    }

    return sum;
}

void FullSweepBuffers::resize(const Mesh& mesh) {
    eval = FixedVector<f64>(mesh.node_count * eval_size);
    jac = FixedVector<f64>(mesh.node_count * jac_size);
    hes = FixedVector<f64>(mesh.node_count * hes_size);
}

FunctionLFG& access_fLg_from_row(FullSweepLayout& layout_lfg, int row) {
    int f_size = layout_lfg.f.int_size();
    int L_size = layout_lfg.L ? 1 : 0;

    if (row < f_size) {
        return layout_lfg.f[row];
    }
    else if (layout_lfg.L && row == f_size) {
        return *layout_lfg.L;
    }
    else {
        return layout_lfg.g[row - f_size - L_size];
    }
}

FunctionLFG& access_Lfg_from_row(FullSweepLayout& layout_lfg, int row) {
    int f_size = layout_lfg.f.int_size();
    int L_size = layout_lfg.L ? 1 : 0;

    if (layout_lfg.L && row == 0) {
        return *layout_lfg.L;
    }
    else if (row < f_size + L_size) {
        return layout_lfg.f[row - L_size];
    }
    else {
        return layout_lfg.g[row - f_size - L_size];
    }
}

void FullSweep::print_jacobian_sparsity_pattern() {
    FixedTableFormat<4> table_format = {{18, 15, 16, 14}, {Align::Center, Align::Center, Align::Center, Align::Center}};

    Log::start_module(table_format, "FullSweep (LFG) - Jacobian Sparsity");
    Log::row(table_format, "Function", "Variable Type", "Variable Index", "Buffer Index");
    Log::dashes(table_format);

    if (layout.L) {
        for (const auto& entry : layout.L->jac.dx) {
            Log::row(table_format, "Lagrange - L", "dx", entry.col, entry.buf_index);
        }
        for (const auto& entry : layout.L->jac.du) {
            Log::row(table_format, "Lagrange - L", "du", entry.col, entry.buf_index);
        }
        for (const auto& entry : layout.L->jac.dp) {
            Log::row(table_format, "Lagrange - L", "dp", entry.col, entry.buf_index);
        }
        Log::dashes(table_format);
    }

    for (size_t f = 0; f < layout.f.size(); f++) {
        const auto& current_f = layout.f[f];

        for (const auto& entry : current_f.jac.dx) {
            Log::row(table_format, fmt::format("Dynamic - f[{}]", f), "dx", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_f.jac.du) {
            Log::row(table_format, fmt::format("Dynamic - f[{}]", f), "du", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_f.jac.dp) {
            Log::row(table_format, fmt::format("Dynamic - f[{}]", f), "dp", entry.col, entry.buf_index);
        }
        Log::dashes(table_format);
    }

    for (size_t g = 0; g < layout.g.size(); g++) {
        const auto& current_g = layout.g[g];

        for (const auto& entry : current_g.jac.dx) {
            Log::row(table_format, fmt::format("Path - g[{}]", g), "dx", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_g.jac.du) {
            Log::row(table_format, fmt::format("Path - g[{}]", g), "du", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_g.jac.dp) {
            Log::row(table_format, fmt::format("Path - g[{}]", g), "dp", entry.col, entry.buf_index);
        }
        Log::dashes(table_format);
    }
    Log::info("");
}

void BoundarySweep::print_jacobian_sparsity_pattern() {
    FixedTableFormat<4> table_format = {{18, 15, 16, 14}, {Align::Center, Align::Center, Align::Center, Align::Center}};

    Log::start_module(table_format, "BoundarySweep (MR) - Jacobian Sparsity");
    Log::row(table_format, "Function", "Variable Type", "Variable Index", "Buffer Index");
    Log::dashes(table_format);

    if (layout.M) {
        for (const auto& entry : layout.M->jac.dx0) {
            Log::row(table_format, "Mayer - M", "dx0", entry.col, entry.buf_index);
        }
        for (const auto& entry : layout.M->jac.dxf) {
            Log::row(table_format, "Mayer - M", "dxf", entry.col, entry.buf_index);
        }
        for (const auto& entry : layout.M->jac.dp) {
            Log::row(table_format, "Mayer - M", "dp", entry.col, entry.buf_index);
        }
        for (const auto& entry : layout.M->jac.dT) {
            Log::row(table_format, "Mayer - M", "dT", entry.col, entry.buf_index);
        }
        Log::dashes(table_format);
    }

    for (size_t r = 0; r < layout.r.size(); r++) {
        const auto& current_r = layout.r[r];

        for (const auto& entry : current_r.jac.dx0) {
            Log::row(table_format, fmt::format("Boundary - r[{}]", r), "dx0", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_r.jac.dxf) {
            Log::row(table_format, fmt::format("Boundary - r[{}]", r), "dxf", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_r.jac.dp) {
            Log::row(table_format, fmt::format("Boundary - r[{}]", r), "dp", entry.col, entry.buf_index);
        }
        for (const auto& entry : current_r.jac.dT) {
            Log::row(table_format, fmt::format("Boundary - r[{}]", r), "dT", entry.col, entry.buf_index);
        }
        Log::dashes(table_format);
    }
    Log::info("");
}

} // namespace GDOP
