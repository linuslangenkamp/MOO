// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2026 University of Applied Sciences and Arts
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

#include <algorithm>

namespace CInit {

namespace {

FixedVector<f64> copy_values(f64 *values, int size) {
    FixedVector<f64> out(size);
    for (int i = 0; i < size; i++) {
        out[i] = values ? values[i] : 0.0;
    }
    return out;
}

FixedVector<Bounds> copy_bounds(bounds_t *bounds, int size) {
    FixedVector<Bounds> out(size);
    for (int i = 0; i < size; i++) {
        out[i].lb = bounds ? bounds[i].lb : MINUS_INFINITY;
        out[i].ub = bounds ? bounds[i].ub : PLUS_INFINITY;
    }
    return out;
}

FixedVector<int> copy_ints(int *values, int size) {
    FixedVector<int> out(size);
    for (int i = 0; i < size; i++) {
        out[i] = values[i];
    }
    return out;
}

} // namespace

Problem::Problem(c_init_problem_t *c_problem)
    : c_problem(c_problem),
      z_buffer(c_problem->y_size + c_problem->p_size),
      eval_buffer(1 + c_problem->f_size + c_problem->g_size),
      jac_buffer((c_problem->obj_jac ? c_problem->obj_jac->nnz : 0) + (c_problem->f_jac ? c_problem->f_jac->nnz : 0) + (c_problem->g_jac ? c_problem->g_jac->nnz : 0)),
      lambda_buffer(c_problem->f_size + c_problem->g_size),
      hes_buffer(c_problem->hes ? c_problem->hes->nnz : 0) {
    formulation.y_size = c_problem->y_size;
    formulation.p_size = c_problem->p_size;
    formulation.f_size = c_problem->f_size;
    formulation.g_size = c_problem->g_size;
    formulation.objective = Init::Objective::USER;

    formulation.y0 = copy_values(c_problem->y0, c_problem->y_size);
    formulation.p0 = copy_values(c_problem->p0, c_problem->p_size);
    formulation.dp0 = copy_values(c_problem->dp0, c_problem->p_size);

    formulation.y_bounds = copy_bounds(c_problem->y_bounds, c_problem->y_size);
    formulation.p_bounds = copy_bounds(c_problem->p_bounds, c_problem->p_size);
    formulation.dp_bounds = copy_bounds(c_problem->dp_bounds, c_problem->p_size);
    formulation.g_bounds = copy_bounds(c_problem->g_bounds, c_problem->g_size);

    formulation.y_nominal = copy_values(c_problem->y_nominal, c_problem->y_size);
    formulation.p_nominal = copy_values(c_problem->p_nominal, c_problem->p_size);
    formulation.dp_nominal = copy_values(c_problem->dp_nominal, c_problem->p_size);
    formulation.f_nominal = copy_values(c_problem->f_nominal, c_problem->f_size);
    formulation.g_nominal = copy_values(c_problem->g_nominal, c_problem->g_size);
    formulation.obj_nominal = c_problem->obj_nominal;

    formulation.jac_f_rows = copy_ints(c_problem->f_jac->row, c_problem->f_jac->nnz);
    formulation.jac_f_cols = copy_ints(c_problem->f_jac->col, c_problem->f_jac->nnz);
    formulation.jac_g_rows = copy_ints(c_problem->g_jac->row, c_problem->g_jac->nnz);
    formulation.jac_g_cols = copy_ints(c_problem->g_jac->col, c_problem->g_jac->nnz);
    formulation.hes_rows = copy_ints(c_problem->hes->row, c_problem->hes->nnz);
    formulation.hes_cols = copy_ints(c_problem->hes->col, c_problem->hes->nnz);
}

const f64 *Problem::pack_z(const f64 *y, const f64 *p) {
    for (int i = 0; i < c_problem->y_size; i++) {
        z_buffer[i] = y[i];
    }
    for (int i = 0; i < c_problem->p_size; i++) {
        z_buffer[c_problem->y_size + i] = p[i];
    }
    return z_buffer.raw();
}

void Problem::eval_objective(const f64 *y, const f64 *p, f64 &obj) {
    c_problem->callbacks->eval_all(pack_z(y, p), c_problem->rp, eval_buffer.raw(), c_problem->user_data);
    obj = eval_buffer[0];
}

void Problem::eval_grad_objective(const f64 *y, const f64 *p, f64 *grad_y, f64 *grad_p) {
    std::fill(grad_y, grad_y + c_problem->y_size, 0.0);
    std::fill(grad_p, grad_p + c_problem->p_size, 0.0);
    c_problem->callbacks->jac_all(pack_z(y, p), c_problem->rp, jac_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->obj_jac->nnz; nz++) {
        const int col = c_problem->obj_jac->col[nz];
        if (col < c_problem->y_size) {
            grad_y[col] = jac_buffer[c_problem->obj_jac->buf_index[nz]];
        } else {
            grad_p[col - c_problem->y_size] = jac_buffer[c_problem->obj_jac->buf_index[nz]];
        }
    }
}

void Problem::eval_hessian_objective(const f64 *y, const f64 *p, f64 obj_factor, f64 *hes_values) {
    lambda_buffer.fill_zero();
    c_problem->callbacks->hes_all(pack_z(y, p), c_problem->rp, lambda_buffer.raw(), obj_factor, hes_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->hes->nnz; nz++) {
        hes_values[nz] += hes_buffer[nz];
    }
}

void Problem::eval_f(const f64 *y, const f64 *p, f64 *f) {
    c_problem->callbacks->eval_all(pack_z(y, p), c_problem->rp, eval_buffer.raw(), c_problem->user_data);
    for (int i = 0; i < c_problem->f_size; i++) {
        f[i] = eval_buffer[1 + i];
    }
}

void Problem::eval_g(const f64 *y, const f64 *p, f64 *g) {
    c_problem->callbacks->eval_all(pack_z(y, p), c_problem->rp, eval_buffer.raw(), c_problem->user_data);
    for (int i = 0; i < c_problem->g_size; i++) {
        g[i] = eval_buffer[1 + c_problem->f_size + i];
    }
}

void Problem::eval_jacobian_f(const f64 *y, const f64 *p, f64 *jac_f_values) {
    c_problem->callbacks->jac_all(pack_z(y, p), c_problem->rp, jac_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->f_jac->nnz; nz++) {
        jac_f_values[nz] = jac_buffer[c_problem->f_jac->buf_index[nz]];
    }
}

void Problem::eval_jacobian_g(const f64 *y, const f64 *p, f64 *jac_g_values) {
    c_problem->callbacks->jac_all(pack_z(y, p), c_problem->rp, jac_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->g_jac->nnz; nz++) {
        jac_g_values[nz] = jac_buffer[c_problem->g_jac->buf_index[nz]];
    }
}

void Problem::eval_hessian_constraints(const f64 *y, const f64 *p, const f64 *lambda_f, const f64 *lambda_g, f64 *hes_values) {
    for (int i = 0; i < c_problem->f_size; i++) {
        lambda_buffer[i] = lambda_f[i];
    }
    for (int i = 0; i < c_problem->g_size; i++) {
        lambda_buffer[c_problem->f_size + i] = lambda_g ? lambda_g[i] : 0.0;
    }
    c_problem->callbacks->hes_all(pack_z(y, p), c_problem->rp, lambda_buffer.raw(), 0.0, hes_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->hes->nnz; nz++) {
        hes_values[nz] += hes_buffer[nz];
    }
}

} // namespace CInit
