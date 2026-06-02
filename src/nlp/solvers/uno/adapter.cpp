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

#include "adapter.h"

#include <base/log.h>

#include <chrono>

namespace UnoSolver {

namespace {

template <typename Fn>
uno_int timed_callback(f64& total_nano, int& count, Fn&& fn)
{
    const auto start = std::chrono::high_resolution_clock::now();
    fn();
    total_nano += std::chrono::duration<f64, std::nano>(
        std::chrono::high_resolution_clock::now() - start).count();
    count++;
    return 0;
}

} // namespace

void UnoCallbackTiming::reset()
{
    *this = {};
}

f64 UnoCallbackTiming::total_nano() const
{
    return objective_nano + objective_gradient_nano + constraints_nano + jacobian_nano + hessian_nano;
}

UnoAdapter::UnoAdapter(NLP::NLP& nlp)
    : nlp(nlp)
{}

UnoAdapter& UnoAdapter::from_user_data(void* user_data)
{
    return *static_cast<UnoAdapter*>(user_data);
}

void* UnoAdapter::create_model()
{
    nlp.solver_get_info(number_vars, number_constraints, nnz_jac, nnz_hes);

    x_lb = FixedVector<f64>(number_vars);
    x_ub = FixedVector<f64>(number_vars);
    x_init = FixedVector<f64>(number_vars);
    z_lb_solution = FixedVector<f64>(number_vars);
    z_ub_solution = FixedVector<f64>(number_vars);
    x_solution = FixedVector<f64>(number_vars);

    g_lb = FixedVector<f64>(number_constraints);
    g_ub = FixedVector<f64>(number_constraints);
    lambda_init = FixedVector<f64>(number_constraints);
    lambda_solution = FixedVector<f64>(number_constraints);

    i_row_jac = FixedVector<uno_int>(nnz_jac);
    j_col_jac = FixedVector<uno_int>(nnz_jac);
    i_row_hes = FixedVector<uno_int>(nnz_hes);
    j_col_hes = FixedVector<uno_int>(nnz_hes);

    nlp.solver_get_bounds(x_lb.raw(), x_ub.raw(), g_lb.raw(), g_ub.raw());
    nlp.solver_get_initial_guess(true, x_init.raw(), true, lambda_init.raw(), false, nullptr, nullptr);
    nlp.solver_get_jac_sparsity(i_row_jac.raw(), j_col_jac.raw());
    nlp.solver_get_hes_sparsity(i_row_hes.raw(), j_col_hes.raw());

    void* model = uno_create_model(UNO_PROBLEM_NONLINEAR, number_vars, x_lb.raw(), x_ub.raw(), UNO_ZERO_BASED_INDEXING);
    if (!model) {
        Log::error("[Uno Interface] Failed to create Uno model.");
        return nullptr;
    }

    bool ok = true;
    ok = ok && uno_set_user_data(model, this);
    ok = ok && uno_set_objective(model, UNO_MINIMIZE, &UnoAdapter::eval_f, &UnoAdapter::eval_grad_f);
    ok = ok && uno_set_lagrangian_sign_convention(model, UNO_MULTIPLIER_POSITIVE);
    ok = ok && uno_set_initial_primal_iterate(model, x_init.raw());
    ok = ok && uno_set_initial_dual_iterate(model, lambda_init.raw());

    if (number_constraints > 0) {
        ok = ok && uno_set_constraints(model, number_constraints, &UnoAdapter::eval_g, g_lb.raw(), g_ub.raw(),
            nnz_jac, i_row_jac.raw(), j_col_jac.raw(), &UnoAdapter::eval_jac);
    }

    if (nnz_hes > 0) {
        ok = ok && uno_set_lagrangian_hessian(model, nnz_hes, UNO_LOWER_TRIANGLE,
            i_row_hes.raw(), j_col_hes.raw(), &UnoAdapter::eval_hes);
    }

    if (!ok) {
        Log::error("[Uno Interface] Failed to populate Uno model.");
        uno_destroy_model(model);
        return nullptr;
    }

    return model;
}

void UnoAdapter::reset_timing()
{
    callback_timing.reset();
}

void UnoAdapter::finalize_solution(void* solver, NLP::ReturnCode return_code)
{
    uno_get_primal_solution(solver, x_solution.raw());
    uno_get_constraint_dual_solution(solver, lambda_solution.raw());
    uno_get_lower_bound_dual_solution(solver, z_lb_solution.raw());
    uno_get_upper_bound_dual_solution(solver, z_ub_solution.raw());

    const double objective = uno_get_solution_objective(solver);
    nlp.solver_finalize_solution(return_code, objective, x_solution.raw(), lambda_solution.raw(), z_lb_solution.raw(), z_ub_solution.raw());
}

uno_int UnoAdapter::eval_f(uno_int, const double* x, double* objective_value, void* user_data)
{
    auto& adapter = from_user_data(user_data);
    return timed_callback(adapter.callback_timing.objective_nano, adapter.callback_timing.objective_count, [&] {
        adapter.nlp.solver_eval_f(true, x, *objective_value);
    });
}

uno_int UnoAdapter::eval_grad_f(uno_int, const double* x, double* gradient, void* user_data)
{
    auto& adapter = from_user_data(user_data);
    return timed_callback(adapter.callback_timing.objective_gradient_nano, adapter.callback_timing.objective_gradient_count, [&] {
        adapter.nlp.solver_eval_grad_f(true, x, gradient);
    });
}

uno_int UnoAdapter::eval_g(uno_int, uno_int, const double* x, double* constraint_values, void* user_data)
{
    auto& adapter = from_user_data(user_data);
    return timed_callback(adapter.callback_timing.constraints_nano, adapter.callback_timing.constraints_count, [&] {
        adapter.nlp.solver_eval_g(true, x, constraint_values);
    });
}

uno_int UnoAdapter::eval_jac(uno_int, uno_int, const double* x, double* jacobian_values, void* user_data)
{
    auto& adapter = from_user_data(user_data);
    return timed_callback(adapter.callback_timing.jacobian_nano, adapter.callback_timing.jacobian_count, [&] {
        adapter.nlp.solver_eval_jac(true, x, jacobian_values);
    });
}

uno_int UnoAdapter::eval_hes(uno_int, uno_int, uno_int, const double* x, double objective_multiplier,
    const double* multipliers, double* hessian_values, void* user_data)
{
    auto& adapter = from_user_data(user_data);
    return timed_callback(adapter.callback_timing.hessian_nano, adapter.callback_timing.hessian_count, [&] {
        adapter.nlp.solver_eval_hes(true, x, true, multipliers, objective_multiplier, hessian_values);
    });
}

} // namespace UnoSolver
