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

#ifndef OPT_NLP_SOLVER_H
#define OPT_NLP_SOLVER_H

#include "nlp/solvers/nlp_solver_settings.h"
#include "nlp.h"

namespace NLP {

class NLPSolver {
public:
    NLPSolver(NLP& nlp, NLPSolverSettings& solver_settings);

    virtual ~NLPSolver() = default;

    NLP& nlp;
    NLPSolverSettings& solver_settings;

    // optimization entry point
    virtual void optimize() = 0;

    // receive optimizer iteration count
    virtual int get_iterations() const;

    // times in nanoseconds
    virtual f64 get_total_time() const;
    virtual f64 get_solver_time() const;
    virtual f64 get_callback_time() const;
};

}

#endif // OPT_NLP_SOLVER_H
