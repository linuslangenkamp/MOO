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

#include <cmath>
#include <iostream>

#include <base/util.h>
#include <nlp/nlp.h>
#include <nlp/solvers/nlp_solver_settings.h>

#if defined(MOO_TEST_SOLVER_IPOPT)
#include <nlp/solvers/ipopt/solver.h>
#elif defined(MOO_TEST_SOLVER_UNO)
#include <nlp/solvers/uno/solver.h>
#else
#error "Define one solver for this test."
#endif

class TrivialQP : public NLP::NLP {
public:
    ::NLP::ReturnCode ret = ::NLP::ReturnCode::GENERIC_FAILURE;
    f64 objective = 0.0;
    FixedVector<f64> solution;

    TrivialQP()
        : solution(2)
    {}

    void get_sizes(int& number_vars, int& number_constraints) override
    {
        number_vars = 2;
        number_constraints = 1;
    }

    void get_nnz(int& nnz_jac, int& nnz_hes) override
    {
        nnz_jac = 2;
        nnz_hes = 2;
    }

    void get_bounds(FixedVector<f64>& x_lb, FixedVector<f64>& x_ub,
                    FixedVector<f64>& g_lb, FixedVector<f64>& g_ub) override
    {
        x_lb[0] = 0.0;
        x_lb[1] = 0.0;
        x_ub[0] = PLUS_INFINITY;
        x_ub[1] = PLUS_INFINITY;
        g_lb[0] = 1.0;
        g_ub[0] = 1.0;
    }

    void get_initial_guess(bool init_x, FixedVector<f64>& x_init,
                           bool init_lambda, FixedVector<f64>& lambda_init,
                           bool init_z, FixedVector<f64>& z_lb_init,
                           FixedVector<f64>& z_ub_init) override
    {
        if (init_x) {
            x_init[0] = 0.25;
            x_init[1] = 0.75;
        }
        if (init_lambda) {
            lambda_init[0] = 0.0;
        }
        if (init_z) {
            z_lb_init[0] = 0.0;
            z_lb_init[1] = 0.0;
            z_ub_init[0] = 0.0;
            z_ub_init[1] = 0.0;
        }
    }

    void get_jac_sparsity(FixedVector<int>& i_row_jac, FixedVector<int>& j_col_jac) override
    {
        i_row_jac[0] = 0;
        j_col_jac[0] = 0;
        i_row_jac[1] = 0;
        j_col_jac[1] = 1;
    }

    void get_hes_sparsity(FixedVector<int>& i_row_hes, FixedVector<int>& j_col_hes) override
    {
        i_row_hes[0] = 0;
        j_col_hes[0] = 0;
        i_row_hes[1] = 1;
        j_col_hes[1] = 1;
    }

    void eval_f(bool, const FixedVector<f64>& x, f64& obj) override
    {
        obj = x[0] * x[0] + x[1] * x[1];
    }

    void eval_g(bool, const FixedVector<f64>& x, FixedVector<f64>& g) override
    {
        g[0] = x[0] + x[1];
    }

    void eval_grad_f(bool, const FixedVector<f64>& x, FixedVector<f64>& grad_f) override
    {
        grad_f[0] = 2.0 * x[0];
        grad_f[1] = 2.0 * x[1];
    }

    void eval_jac_g(bool, const FixedVector<f64>&,
                    const FixedVector<int>&, const FixedVector<int>&,
                    FixedVector<f64>& jac) override
    {
        jac[0] = 1.0;
        jac[1] = 1.0;
    }

    void eval_hes(bool, const FixedVector<f64>&,
                  bool, const FixedVector<f64>&,
                  f64 obj_factor,
                  const FixedVector<int>&, const FixedVector<int>&,
                  FixedVector<f64>& hes) override
    {
        hes[0] = 2.0 * obj_factor;
        hes[1] = 2.0 * obj_factor;
    }

    void finalize_solution(::NLP::ReturnCode return_code, f64 opt_obj,
                           const FixedVector<f64>& opt_x,
                           const FixedVector<f64>&,
                           const FixedVector<f64>&,
                           const FixedVector<f64>&) override
    {
        ret = return_code;
        objective = opt_obj;
        solution = opt_x;
    }
};

bool close(f64 value, f64 expected, f64 tolerance)
{
    return std::abs(value - expected) <= tolerance;
}

int main(int, char**)
{
    TrivialQP qp;

    char argv0[] = "test_solver_trivial_qp";
    char* argv[] = {argv0};
    NLP::NLPSolverSettings settings(1, argv);
    settings.set(NLP::Option::Tolerance, 1e-12);
    settings.set(NLP::Option::Iterations, 100);
    settings.set(NLP::Option::CPUTime, 30.0);
    settings.set(NLP::Option::Hessian, NLP::HessianOption::Exact);

#if defined(MOO_TEST_SOLVER_IPOPT)
    settings.set(NLP::Option::NLPSolver, NLP::NLPSolverOption::Ipopt);
    IpoptSolver::IpoptSolver solver(qp, settings);
#elif defined(MOO_TEST_SOLVER_UNO)
    settings.set(NLP::Option::NLPSolver, NLP::NLPSolverOption::Uno);
    settings.set(NLP::Option::UnoPreset, std::string("ipopt"));
    UnoSolver::UnoSolver solver(qp, settings);
#endif

    solver.optimize();

    const bool status_ok = qp.ret == ::NLP::ReturnCode::OPTIMAL || qp.ret == ::NLP::ReturnCode::ACCEPTABLE;
    const bool solution_ok = close(qp.solution[0], 0.5, 1e-6)
                          && close(qp.solution[1], 0.5, 1e-6)
                          && close(qp.objective, 0.5, 1e-6);

    if (!status_ok || !solution_ok) {
        std::cerr << "Unexpected trivial QP result: ret=" << static_cast<int>(qp.ret)
                  << ", objective=" << qp.objective
                  << ", x=" << qp.solution[0]
                  << ", y=" << qp.solution[1] << '\n';
        return 1;
    }

    return 0;
}
