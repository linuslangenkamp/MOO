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

#ifndef MOO_UNO_SOLVER_H
#define MOO_UNO_SOLVER_H

#include <memory>

#include <base/export.h>
#include <base/timing.h>
#include <nlp/nlp_solver.h>

namespace UnoSolver {

struct UnoSolverData;
class UnoTimingNode;

class MOO_EXPORT UnoSolver : public NLP::NLPSolver {
    friend UnoTimingNode;

public:
    UnoSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings);
    ~UnoSolver() override;

    void optimize() override;
    void set_settings();

    int get_iterations() const override;
    f64 get_total_time() const override;
    f64 get_solver_time() const override;
    f64 get_callback_time() const override;

private:
    UnoSolverData* udata;

    NLP::ReturnCode get_return_code() const;
    void log_status() const;
};

class UnoTimingNode : public TimingNode {
    UnoSolver* uno_solver;

public:
    UnoTimingNode(std::string n, TimingNode* p = nullptr, UnoSolver* uno_solver = nullptr);

    void finalize() override;
};

} // namespace UnoSolver

#endif // MOO_UNO_SOLVER_H
