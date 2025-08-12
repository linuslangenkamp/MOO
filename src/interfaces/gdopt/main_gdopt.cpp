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

#include <nlp/instances/gdop/orchestrator.h>
#include <nlp/solvers/ipopt/solver.h>
#include <nlp/instances/gdop/strategies.h>
#include <base/log.h>
#include <interfaces/c/problem.h>

#include "main_gdopt.h"

extern f64* _rp;

// create config for the algorithm (for now basic)
class Config {

};

Config read_yaml() {
    return Config();
}

void set_global_configuration(Config& config) {

}

void set_global_runtime_parameters(Config& config) {
    /* set _rp */
}

int main_gdopt(int argc, char** argv, c_problem_t* c_problem) {
    LOG_PREFIX('*', "Entry point [OPT] - main_gdopt\n");

    auto config = read_yaml();
    set_global_configuration(config);

    set_global_runtime_parameters(config);

    c_problem->callbacks->update_c_problem(NULL);

    auto nlp_solver_settings = NLP::NLPSolverSettings(argc, argv);
    nlp_solver_settings.print();

    //nlp_solver_settings.set(NLP::Option::IpoptDerivativeTest, true);

    auto mesh = Mesh::create_equidistant_fixed_stages(1 /* tf */, 25 /* intervals */, 3 /* stages */);
    auto problem = C::create_gdop(c_problem, *mesh);

    auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    auto gdop = GDOP::GDOP(problem);

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_settings);

    auto orchestrator = GDOP::MeshRefinementOrchestrator(gdop, std::move(strategies), ipopt_solver);

    orchestrator.optimize();

    return 0;
}
