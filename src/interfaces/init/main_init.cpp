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

#include <interfaces/init/main_init.h>

#include <base/log.h>
#include <nlp/instances/init/init.h>
#include <nlp/solvers/ipopt/solver.h>
#include <nlp/solvers/nlp_solver_settings.h>

#if MOO_WITH_UNO_ENABLED
#include <nlp/solvers/uno/solver.h>
#endif

#include "problem.h"

#include <fstream>
#include <iomanip>
#include <memory>

namespace {

void write_init_result_csv(const Init::Result &result, const Init::ProblemFormulation &form) {
    std::ofstream out("init_optimal_solution.csv");
    out << std::setprecision(17);

    out << "objective,f_l2_norm,f_max_error,g_max_violation";
    for (int i = 0; i < form.y_size; i++) {
        out << ",y_" << i;
    }
    for (int i = 0; i < form.p_size; i++) {
        out << ",p_" << i;
    }
    for (int i = 0; i < form.p_size; i++) {
        out << ",dp_" << i;
    }
    for (int i = 0; i < form.f_size; i++) {
        out << ",f_" << i;
    }
    for (int i = 0; i < form.g_size; i++) {
        out << ",g_" << i;
    }
    for (int i = 0; i < result.lambda.int_size(); i++) {
        out << ",lambda_" << i;
    }
    for (int i = 0; i < result.z_lb.int_size(); i++) {
        out << ",z_lb_" << i;
    }
    for (int i = 0; i < result.z_ub.int_size(); i++) {
        out << ",z_ub_" << i;
    }
    out << "\n";

    out << result.objective << "," << result.f_l2_norm << "," << result.f_max_error << "," << result.g_max_violation;
    for (int i = 0; i < form.y_size; i++) {
        out << "," << result.y[i];
    }
    for (int i = 0; i < form.p_size; i++) {
        out << "," << result.p[i];
    }
    for (int i = 0; i < form.p_size; i++) {
        out << "," << result.dp[i];
    }
    for (int i = 0; i < form.f_size; i++) {
        out << "," << result.constraints[i];
    }
    for (int i = 0; i < form.g_size; i++) {
        out << "," << result.constraints[form.f_size + i];
    }
    for (int i = 0; i < result.lambda.int_size(); i++) {
        out << "," << result.lambda[i];
    }
    for (int i = 0; i < result.z_lb.int_size(); i++) {
        out << "," << result.z_lb[i];
    }
    for (int i = 0; i < result.z_ub.int_size(); i++) {
        out << "," << result.z_ub[i];
    }
    out << "\n";
}

} // namespace

int main_init(int argc, char **argv, c_init_problem_t *c_problem) {
    Log::prefixed('*', "Entry point [INIT] - main_init\n");

    auto nlp_solver_settings = NLP::NLPSolverSettings(argc, argv);
    if (c_problem->solver_ctx && c_problem->solver_ctx->derivative_test) {
        nlp_solver_settings.set(NLP::Option::IpoptDerivativeTest, true);
    }
    nlp_solver_settings.print();

    CInit::Problem problem(c_problem);
    Init::Init init(problem);

    std::unique_ptr<NLP::NLPSolver> nlp_solver;
    switch (nlp_solver_settings.get_or_default<NLP::NLPSolverOption>(NLP::Option::NLPSolver)) {
        case NLP::NLPSolverOption::Ipopt:
            nlp_solver = std::make_unique<IpoptSolver::IpoptSolver>(init, nlp_solver_settings);
            break;
        case NLP::NLPSolverOption::Uno:
#if MOO_WITH_UNO_ENABLED
            nlp_solver = std::make_unique<UnoSolver::UnoSolver>(init, nlp_solver_settings);
#else
            Log::error("Uno solver requested, but MOO was built with MOO_WITH_UNO=OFF.");
            std::abort();
#endif
            break;
    }

    nlp_solver->optimize();
    write_init_result_csv(init.get_result(), problem.formulation);
    return 0;
}
