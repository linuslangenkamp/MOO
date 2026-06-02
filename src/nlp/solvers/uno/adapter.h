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

#ifndef MOO_UNO_ADAPTER_H
#define MOO_UNO_ADAPTER_H

#include <base/fixed_vector.h>
#include <interfaces/C/Uno_C_API.h>

#include <nlp/nlp.h>

namespace UnoSolver {

struct UnoCallbackTiming {
    f64 objective_nano = 0.0;
    f64 objective_gradient_nano = 0.0;
    f64 constraints_nano = 0.0;
    f64 jacobian_nano = 0.0;
    f64 hessian_nano = 0.0;

    int objective_count = 0;
    int objective_gradient_count = 0;
    int constraints_count = 0;
    int jacobian_count = 0;
    int hessian_count = 0;

    void reset();
    f64 total_nano() const;
};

class UnoAdapter {
public:
    explicit UnoAdapter(NLP::NLP& nlp);

    void* create_model();
    void finalize_solution(void* solver, NLP::ReturnCode return_code);
    void reset_timing();

    int get_number_vars() const { return number_vars; }
    int get_number_constraints() const { return number_constraints; }
    int get_nnz_jac() const { return nnz_jac; }
    int get_nnz_hes() const { return nnz_hes; }
    const UnoCallbackTiming& get_callback_timing() const { return callback_timing; }

private:
    NLP::NLP& nlp;
    UnoCallbackTiming callback_timing;

    int number_vars = 0;
    int number_constraints = 0;
    int nnz_jac = 0;
    int nnz_hes = 0;

    FixedVector<f64> x_lb;
    FixedVector<f64> x_ub;
    FixedVector<f64> g_lb;
    FixedVector<f64> g_ub;
    FixedVector<f64> x_init;
    FixedVector<f64> lambda_init;
    FixedVector<f64> z_lb_solution;
    FixedVector<f64> z_ub_solution;
    FixedVector<f64> x_solution;
    FixedVector<f64> lambda_solution;
    FixedVector<uno_int> i_row_jac;
    FixedVector<uno_int> j_col_jac;
    FixedVector<uno_int> i_row_hes;
    FixedVector<uno_int> j_col_hes;

    static UnoAdapter& from_user_data(void* user_data);

    static uno_int eval_f(uno_int number_variables, const double* x, double* objective_value, void* user_data);
    static uno_int eval_grad_f(uno_int number_variables, const double* x, double* gradient, void* user_data);
    static uno_int eval_g(uno_int number_variables, uno_int number_constraints, const double* x, double* constraint_values, void* user_data);
    static uno_int eval_jac(uno_int number_variables, uno_int number_jacobian_nonzeros, const double* x, double* jacobian_values, void* user_data);
    static uno_int eval_hes(uno_int number_variables, uno_int number_constraints, uno_int number_hessian_nonzeros,
        const double* x, double objective_multiplier, const double* multipliers, double* hessian_values, void* user_data);
};

} // namespace UnoSolver

#endif // MOO_UNO_ADAPTER_H
