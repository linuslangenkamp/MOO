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

#ifndef MOO_AD_VM_H
#define MOO_AD_VM_H

#include "optimize.h"

namespace ad {

struct EvalEnv {
    std::unordered_map<std::string, const double *> inputs;
    std::unordered_map<std::string, const double *> params;

    EvalEnv &input(const std::string &name, const double *data);
    EvalEnv &param(const std::string &name, const double *data);
};

struct VM {
    GraphFunction f;

    explicit VM(GraphFunction fn);
    void evaluate(const EvalEnv &env, double *out) const;

private:
    double eval_scalar_node(NodeId i, const EvalEnv &env, std::vector<double> &mem, std::vector<char> &seen) const;
    std::vector<double> eval_vector_node(int vector_id, const EvalEnv &env, std::vector<double> &mem, std::vector<char> &seen) const;
};

struct StagedVM {
    GraphFunction f;
    std::string direction_name;
    std::vector<char> depends_on_direction;

    explicit StagedVM(GraphFunction fn, std::string direction = "v");

    void analyze();

    struct Prepared {
        const StagedVM *vm = nullptr;
        std::vector<double> mem;

        void apply(const double *direction, double *out) const;
    };

    Prepared prepare(const EvalEnv &env) const;
    void eval_dependent(const EvalEnv &env, std::vector<double> &mem, double *out) const;
    void eval_node(NodeId i, const EvalEnv &env, std::vector<double> &mem) const;
};

} // namespace ad

#endif // MOO_AD_VM_H
