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
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <interfaces/nlp/main_nlp.h>

#include <base/log.h>
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

void write_nlp_result_csv(const CNLP::Result &result, int x_size, int g_size) {
    std::ofstream out("nlp_optimal_solution.csv");
    out << std::setprecision(17);
    out << "objective";
    for (int i = 0; i < x_size; i++) {
        out << ",x_" << i;
    }
    for (int i = 0; i < g_size; i++) {
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

    out << result.objective;
    for (int i = 0; i < x_size; i++) {
        out << "," << result.x[i];
    }
    for (int i = 0; i < g_size; i++) {
        out << "," << result.g[i];
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

int main_nlp(int argc, char **argv, c_nlp_problem_t *c_problem) {
    Log::prefixed('*', "Entry point [NLP] - main_nlp\n");

    auto nlp_solver_settings = NLP::NLPSolverSettings(argc, argv);
    if (c_problem->solver_ctx && c_problem->solver_ctx->derivative_test) {
        nlp_solver_settings.set(NLP::Option::IpoptDerivativeTest, true);
    }
    nlp_solver_settings.print();

    CNLP::Problem problem(c_problem);

    std::unique_ptr<NLP::NLPSolver> nlp_solver;
    switch (nlp_solver_settings.get_or_default<NLP::NLPSolverOption>(NLP::Option::NLPSolver)) {
        case NLP::NLPSolverOption::Ipopt:
            nlp_solver = std::make_unique<IpoptSolver::IpoptSolver>(problem, nlp_solver_settings);
            break;
        case NLP::NLPSolverOption::Uno:
#if MOO_WITH_UNO_ENABLED
            nlp_solver = std::make_unique<UnoSolver::UnoSolver>(problem, nlp_solver_settings);
#else
            Log::error("Uno solver requested, but MOO was built with MOO_WITH_UNO=OFF.");
            std::abort();
#endif
            break;
    }

    nlp_solver->optimize();
    write_nlp_result_csv(problem.result, c_problem->x_size, c_problem->g_size);
    return 0;
}
