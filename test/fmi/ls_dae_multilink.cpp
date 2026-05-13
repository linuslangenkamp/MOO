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

#include <interfaces/fmi/problem.h>

int main() {
    FMI::FMISettings settings;
    settings.path = MULTILINK_FMU_PATH;
    settings.modelname = "Multilink";
    settings.t0 = 0.0;
    settings.tf = 100.0;
    settings.intervals = 200;
    settings.stage = 1;
    // 2e6 min max

    f64 factor_u = 10;
    f64 factor_phi = 1;

    settings.control_vrefs.push_back( { 620756992, -2e6, 2e6, 2e6 } );

    settings.parameter_vrefs = {};

    settings.lagrange_vref = nullptr;

    FMI::main_fmi(settings);

    return EXIT_SUCCESS;
}
