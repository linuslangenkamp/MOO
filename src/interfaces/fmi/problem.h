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

#ifndef MOO_FMI_PROBLEM_H
#define MOO_FMI_PROBLEM_H

#include <nlp/instances/gdop/problem.h>
#include <nlp/instances/gdop/gdop.h>

namespace FMI {

class FMIData {
public:
    FMIData(const char* path, const char* modelname);

    void print();
    void eval_point_lfg(const f64* xu, const f64* p, f64 time, f64* out);
    void jac_point_lfg(const f64* xu, const f64* p, f64 time, f64* out);

    std::unique_ptr<struct FMIData_priv> priv;
};

class FullSweep : public GDOP::FullSweep {
public:
    FullSweep(GDOP::FullSweepLayout&& layout_in,
              const GDOP::ProblemConstants& pc,
              FMIData& fmi_data);

    void callback_eval(const f64* xu_nlp, const f64* p) override;
    void callback_jac(const f64* xu_nlp, const f64* p) override;
    void callback_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, const f64* lambda) override;

    FMIData& fmi_data;
};

class BoundarySweep : public GDOP::BoundarySweep {
public:
    BoundarySweep(GDOP::BoundarySweepLayout&& layout_in,
                  const GDOP::ProblemConstants& pc,
                  FMIData& fmi_data);

    void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, f64 t0, f64 tf) override;
    void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, f64 t0, f64 tf) override;
    void callback_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, f64 t0, f64 tf, const f64 mayer_factor, const f64* lambda) override;

    FMIData& fmi_data;
};

class Problem : public GDOP::Problem {
public:
    Problem(FMIData& fmi_data);
};

MOO_EXPORT void main_fmi(const char* path, const char* modelname);

} // namespace FMI

#endif // MOO_FMI_PROBLEM_H
