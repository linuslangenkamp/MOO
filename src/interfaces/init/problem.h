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

#ifndef MOO_INIT_C_PROBLEM_H
#define MOO_INIT_C_PROBLEM_H

#include <interfaces/init/structures.h>
#include <nlp/instances/init/problem.h>

namespace CInit {

class Problem final : public Init::Problem {
public:
    explicit Problem(c_init_problem_t *c_problem);

    void eval_objective(const f64 *y, const f64 *p, f64 &obj) override;
    void eval_grad_objective(const f64 *y, const f64 *p, f64 *grad_y, f64 *grad_p) override;
    void eval_hessian_objective(const f64 *y, const f64 *p, f64 obj_factor, f64 *hes_values) override;
    void eval_f(const f64 *y, const f64 *p, f64 *f) override;
    void eval_g(const f64 *y, const f64 *p, f64 *g) override;
    void eval_jacobian_f(const f64 *y, const f64 *p, f64 *jac_f_values) override;
    void eval_jacobian_g(const f64 *y, const f64 *p, f64 *jac_g_values) override;
    void eval_hessian_constraints(const f64 *y, const f64 *p, const f64 *lambda_f, const f64 *lambda_g, f64 *hes_values) override;

private:
    c_init_problem_t *c_problem;
    FixedVector<f64> z_buffer;
    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> lambda_buffer;
    FixedVector<f64> hes_buffer;

    const f64 *pack_z(const f64 *y, const f64 *p);
};

} // namespace CInit

#endif
