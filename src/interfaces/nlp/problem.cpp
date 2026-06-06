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
#include <cmath>
#include <nlp/nlp_scaling.h>

namespace CNLP {

namespace {

Bounds bounds_or_unbounded(bounds_t* bounds, int idx)
{
    if (!bounds) {
        return Bounds{};
    }
    return Bounds{bounds[idx].lb, bounds[idx].ub};
}

f64 value_or_zero(f64* values, int idx)
{
    return values ? values[idx] : 0.0;
}

f64 nominal_or_one(f64* values, int idx)
{
    if (!values || !std::isfinite(values[idx]) || values[idx] == 0.0) {
        return 1.0;
    }
    return std::abs(values[idx]);
}

} // namespace

Problem::Problem(c_nlp_problem_t* c_problem)
    : c_problem(c_problem),
      eval_buffer(1 + c_problem->g_size),
      jac_buffer((c_problem->obj_jac ? c_problem->obj_jac->nnz : 0) + (c_problem->g_jac ? c_problem->g_jac->nnz : 0)),
      hes_buffer(c_problem->hes ? c_problem->hes->nnz : 0)
{
}

void Problem::get_sizes(int& number_vars, int& number_constraints)
{
    number_vars = c_problem->x_size;
    number_constraints = c_problem->g_size;
}

void Problem::get_nnz(int& nnz_jac, int& nnz_hes)
{
    nnz_jac = c_problem->g_jac ? c_problem->g_jac->nnz : 0;
    nnz_hes = c_problem->hes ? c_problem->hes->nnz : 0;
}

std::shared_ptr<::NLP::Scaling> Problem::get_scaling()
{
    FixedVector<f64> x_nominal(c_problem->x_size);
    FixedVector<f64> g_nominal(c_problem->g_size);
    bool has_scaling = false;

    for (int i = 0; i < c_problem->x_size; i++) {
        x_nominal[i] = nominal_or_one(c_problem->x_nominal, i);
        has_scaling = has_scaling || x_nominal[i] != 1.0;
    }
    for (int i = 0; i < c_problem->g_size; i++) {
        g_nominal[i] = nominal_or_one(c_problem->g_nominal, i);
        has_scaling = has_scaling || g_nominal[i] != 1.0;
    }
    const f64 obj_nominal = std::isfinite(c_problem->obj_nominal) && c_problem->obj_nominal != 0.0
                          ? std::abs(c_problem->obj_nominal)
                          : 1.0;
    has_scaling = has_scaling || obj_nominal != 1.0;

    if (!has_scaling) {
        return std::make_shared<::NLP::NoScaling>();
    }
    return std::make_shared<::NLP::NominalScaling>(std::move(x_nominal), std::move(g_nominal), obj_nominal);
}

void Problem::get_bounds(FixedVector<f64>& x_lb, FixedVector<f64>& x_ub, FixedVector<f64>& g_lb, FixedVector<f64>& g_ub)
{
    for (int i = 0; i < c_problem->x_size; i++) {
        const Bounds b = bounds_or_unbounded(c_problem->x_bounds, i);
        x_lb[i] = b.lb;
        x_ub[i] = b.ub;
    }
    for (int i = 0; i < c_problem->g_size; i++) {
        const Bounds b = bounds_or_unbounded(c_problem->g_bounds, i);
        g_lb[i] = b.lb;
        g_ub[i] = b.ub;
    }
}

void Problem::get_initial_guess(bool init_x, FixedVector<f64>& x_init, bool init_lambda, FixedVector<f64>& lambda_init, bool init_z, FixedVector<f64>& z_lb_init, FixedVector<f64>& z_ub_init)
{
    if (init_x) {
        for (int i = 0; i < c_problem->x_size; i++) {
            x_init[i] = value_or_zero(c_problem->x0, i);
        }
    }
    if (init_lambda) lambda_init.fill_zero();
    if (init_z) {
        z_lb_init.fill_zero();
        z_ub_init.fill_zero();
    }
}

void Problem::get_jac_sparsity(FixedVector<int>& i_row_jac, FixedVector<int>& j_col_jac)
{
    for (int nz = 0; nz < c_problem->g_jac->nnz; nz++) {
        i_row_jac[nz] = c_problem->g_jac->row[nz];
        j_col_jac[nz] = c_problem->g_jac->col[nz];
    }
}

void Problem::get_hes_sparsity(FixedVector<int>& i_row_hes, FixedVector<int>& j_col_hes)
{
    for (int nz = 0; nz < c_problem->hes->nnz; nz++) {
        i_row_hes[nz] = c_problem->hes->row[nz];
        j_col_hes[nz] = c_problem->hes->col[nz];
    }
}

void Problem::eval_f(bool new_x, const FixedVector<f64>& curr_x, f64& curr_obj)
{
    c_problem->callbacks->eval_all(curr_x.raw(), c_problem->rp, eval_buffer.raw(), c_problem->user_data);
    curr_obj = eval_buffer[0];
}

void Problem::eval_g(bool new_x, const FixedVector<f64>& curr_x, FixedVector<f64>& curr_g)
{
    c_problem->callbacks->eval_all(curr_x.raw(), c_problem->rp, eval_buffer.raw(), c_problem->user_data);
    for (int i = 0; i < c_problem->g_size; i++) {
        curr_g[i] = eval_buffer[1 + i];
    }
}

void Problem::eval_grad_f(bool new_x, const FixedVector<f64>& curr_x, FixedVector<f64>& curr_grad_f)
{
    curr_grad_f.fill_zero();
    c_problem->callbacks->jac_all(curr_x.raw(), c_problem->rp, jac_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->obj_jac->nnz; nz++) {
        curr_grad_f[c_problem->obj_jac->col[nz]] = jac_buffer[c_problem->obj_jac->buf_index[nz]];
    }
}

void Problem::eval_jac_g(bool new_x, const FixedVector<f64>& curr_x, const FixedVector<int>& i_row_jac, const FixedVector<int>& j_col_jac, FixedVector<f64>& curr_jac)
{
    c_problem->callbacks->jac_all(curr_x.raw(), c_problem->rp, jac_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->g_jac->nnz; nz++) {
        curr_jac[nz] = jac_buffer[c_problem->g_jac->buf_index[nz]];
    }
}

void Problem::eval_hes(bool new_x, const FixedVector<f64>& curr_x, bool new_lambda, const FixedVector<f64>& curr_lambda, f64 curr_obj_factor, const FixedVector<int>& i_row_hes, const FixedVector<int>& j_col_hes, FixedVector<f64>& curr_hes)
{
    curr_hes.fill_zero();
    c_problem->callbacks->hes_all(curr_x.raw(), c_problem->rp, curr_lambda.raw(), curr_obj_factor, hes_buffer.raw(), c_problem->user_data);
    for (int nz = 0; nz < c_problem->hes->nnz; nz++) {
        curr_hes[nz] = hes_buffer[nz];
    }
}

void Problem::finalize_solution(::NLP::ReturnCode ret, f64 opt_obj, const FixedVector<f64>& opt_x, const FixedVector<f64>& opt_lambda, const FixedVector<f64>& opt_z_lb, const FixedVector<f64>& opt_z_ub)
{
    result.objective = opt_obj;
    result.x = FixedVector<f64>(opt_x);
    result.lambda = FixedVector<f64>(opt_lambda);
    result.z_lb = FixedVector<f64>(opt_z_lb);
    result.z_ub = FixedVector<f64>(opt_z_ub);
    result.g = FixedVector<f64>(c_problem->g_size);
    if (c_problem->g_size > 0) {
        eval_g(false, opt_x, result.g);
    }
}

} // namespace CNLP
