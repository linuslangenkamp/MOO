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

#include <base/log.h>
#include <interfaces/fmi/problem.h>

#include <fmi4c.h>
#include <fmi4c_lsdae.h>

namespace FMI {

struct FMIData_priv {
    fmuHandle* fmu;
    fmi3InstanceHandle* instance;

    std::vector<fmi3ValueReference> state_vrefs;
    std::vector<fmi3ValueReference> deriv_vrefs;
    std::vector<fmi3ValueReference> input_vrefs;
    std::vector<fmi3ValueReference> alg_vrefs;
    std::vector<fmi3ValueReference> res_vrefs;
    std::vector<fmi3ValueReference> output_vrefs;

    FMIData_priv(const char* path, const char* modelname);
};

FMIData_priv::FMIData_priv(const char* path, const char* modelname)
    : fmu(fmi4c_loadFmu(path, modelname)),
      instance(fmi3_instantiateModelExchange(fmu, true, true, nullptr, nullptr))
{}

int FMIData::nx()   const { return priv->state_vrefs.size(); }
int FMIData::ndx()  const { return priv->deriv_vrefs.size(); }
int FMIData::nu()   const { return priv->input_vrefs.size(); }
int FMIData::nw()   const { return priv->alg_vrefs.size(); }
int FMIData::nres() const { return priv->res_vrefs.size(); }
int FMIData::ny()   const { return priv->output_vrefs.size(); }

FMIData::FMIData(const char* path, const char* modelname)
    : priv(std::make_unique<FMIData_priv>(path, modelname))
{
    if (priv->fmu == nullptr) {
        Log::error("fmi4c_loadFmu failed");
        throw std::runtime_error("Failed to load FMU");
    }
    if (!fmiLsDae_isPresent(priv->fmu)) {
        Log::error("fmi-ls-dae manifest not found in FMU");
        fmi4c_freeFmu(priv->fmu);
        throw std::runtime_error("Missing fmi-ls-dae manifest");
    }

    {
        size_t count;
        fmi3_getNumberOfContinuousStates(priv->instance, &count);
        priv->deriv_vrefs.resize(count);
        priv->state_vrefs.resize(count);

        for (int idx = 1; idx < static_cast<int>(count) + 1; idx++) {
            fmi3VariableHandle* var = fmi3_getVariableByIndex(priv->fmu, idx);
            int der = fmi3_getVariableDerivativeIndex(var);
            if (der != 0) {
                priv->state_vrefs[idx - 1] = (fmi3_getVariableValueReference(var));

                fmi3VariableHandle* der_var = fmi3_getVariableByIndex(priv->fmu, der);
                priv->state_vrefs[idx - 1] = (fmi3_getVariableValueReference(der_var));
            }
        }
    }

    {
        int count = fmiLsDae_getNumberOfAlgebraicVariables(priv->fmu);
        priv->alg_vrefs.resize(count);
        for (int idx = 1; idx < count + 1; idx++) {
            fmiLsDaeAlgebraicVariableHandle* h = fmiLsDae_getAlgebraicVariableByIndex(priv->fmu, idx - 1);
            priv->alg_vrefs[idx - 1] = fmiLsDae_getAlgebraicVariableValueReference(h);
        }
    }

    {
        int count = fmiLsDae_getNumberOfResiduals(priv->fmu);
        priv->res_vrefs.resize(count);
        for (int idx = 1; idx < count + 1; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getResidualByIndex(priv->fmu, idx - 1);
            priv->res_vrefs[idx - 1] = fmiLsDae_getValueReference(h);
        }
    }

    {
        int count = fmiLsDae_getNumberOfOutputs(priv->fmu);
        priv->output_vrefs.resize(count);
        for (int idx = 1; idx < count + 1; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getOutputByIndex(priv->fmu, idx - 1);
            priv->output_vrefs[idx + 1] = fmiLsDae_getValueReference(h);
        }
    }

    {
        int nVars = fmi3_getNumberOfVariables(priv->fmu);
        for (int idx = 1; idx < nVars + 1; idx++) {
            fmi3VariableHandle* var = fmi3_getVariableByIndex(priv->fmu, idx);
            if (fmi3_getVariableCausality(var) == fmi3CausalityInput) {
                priv->input_vrefs.push_back(fmi3_getVariableValueReference(var));
            }
        }
    }

    Log::info("FMI version  : {}", static_cast<int>(fmi4c_getFmiVersion(priv->fmu)));
    Log::info("Model name   : {}", fmi3_modelName(priv->fmu));
    Log::info("  nx   = {}  (states)",      nx());
    Log::info("  ndx  = {}  (derivatives)", ndx());
    Log::info("  nu   = {}  (inputs)",      nu());
    Log::info("  nw   = {}  (algebraics)",  nw());
    Log::info("  nres = {}  (residuals)",   nres());
    Log::info("  ny   = {}  (outputs)",     ny());

    // ... [previous code] ...

    Log::info("  ny   = {}  (outputs)",     ny());

    // Completing your partial loop and adding the rest
    Log::info("--- Value References ---");

    Log::info("State Value References:");
    for (auto const& vref : priv->state_vrefs) {
        Log::info("  {}", vref);
    }

    Log::info("Derivative Value References:");
    for (auto const& vref : priv->deriv_vrefs) {
        Log::info("  {}", vref);
    }

    Log::info("Input Value References:");
    for (auto const& vref : priv->input_vrefs) {
        Log::info("  {}", vref);
    }

    Log::info("Algebraic Value References:");
    for (auto const& vref : priv->alg_vrefs) {
        Log::info("  {}", vref);
    }

    Log::info("Residual Value References:");
    for (auto const& vref : priv->res_vrefs) {
        Log::info("  {}", vref);
    }

    Log::info("Output Value References:");
    for (auto const& vref : priv->output_vrefs) {
        Log::info("  {}", vref);
    }
}

void FMIData::eval_point_lfg(const f64* xu, const f64* p, f64 time, f64* out)
{
    // evaluate (L, f, g)(x, u, w, p, time)
}

void FMIData::jac_point_lfg(const f64* xu, const f64* p, f64 time, f64* out)
{
    // evaluate d(L, f, g) / d(x, u, w, p) (x, u, w, p, time)
}

void FullSweep::callback_eval(
    const f64* xu_nlp,
    const f64* p)
{
    fill_zero_eval_buffer();

    for (int i = 0; i < pc.mesh->intervals; i++) {
        for (int j = 0; j < pc.mesh->nodes[i]; j++) {
            const f64* xu_ij = get_xu_ij(xu_nlp, i, j);
            f64* eval_buf_ij = get_eval_buffer(i, j);
            fmi_data.eval_point_lfg(xu_ij, p, pc.mesh->t[i][j], eval_buf_ij);
        }
    }
}

void FullSweep::callback_jac(
    const f64* xu_nlp,
    const f64* p)
{
    fill_zero_jac_buffer();

    for (int i = 0; i < pc.mesh->intervals; i++) {
        for (int j = 0; j < pc.mesh->nodes[i]; j++) {
            const f64* xu_ij = get_xu_ij(xu_nlp, i, j);
            f64* jac_buf_ij = get_jac_buffer(i, j);
            fmi_data.eval_point_lfg(xu_ij, p, pc.mesh->t[i][j], jac_buf_ij);
        }
    }
}

void FullSweep::callback_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, const f64* lambda)
{
    // not yet implemented
}

GDOP::ProblemConstants create_problem_constants(FMIData& fmi_data)
{
    auto x_bounds = FixedVector<Bounds>(2);
    auto u_bounds = FixedVector<Bounds>(4);
    auto p_bounds = FixedVector<Bounds>(0);
    auto T_bounds = std::array<Bounds, 2>{};
    auto g_bounds = FixedVector<Bounds>(2);
    auto r_bounds = FixedVector<Bounds>(0);

    auto xu0_fixed = FixedVector<std::optional<f64>>(6);
    auto xuf_fixed = FixedVector<std::optional<f64>>(6);
    auto T_fixed = std::array<std::optional<f64>, 2>{};

    auto mesh = Mesh::create_equidistant_fixed_stages(
        /* t0 */        0.0,
        /* tf */        1.0,
        /* intervals */ 25,
        /* stages */    3,
        /* type */      MeshType::Physical
    );

    return GDOP::ProblemConstants(
        /* mayer */ false,
        /* lagrange */ false,
        std::move(x_bounds),
        std::move(u_bounds),
        std::move(p_bounds),
        std::move(T_bounds),
        std::move(xu0_fixed),
        std::move(xuf_fixed),
        std::move(T_fixed),
        std::move(r_bounds),
        std::move(g_bounds),
        *mesh
    );
};

GDOP::Problem create_gdop_problem(FMIData& fmi_data)
{
    // create ProblemConstants
    auto problem_constants = std::make_unique<GDOP::ProblemConstants>(create_problem_constants(fmi_data));
    // create bounds
    // create sparsity
    return GDOP::Problem(nullptr, nullptr, nullptr, nullptr);
}

Problem::Problem(FMIData& fmi_data)
    : GDOP::Problem(create_gdop_problem(fmi_data))
{}

void main_fmi(const char* path, const char* modelname)
{
    auto fmi_data = FMIData(path, modelname);
    auto problem = Problem(fmi_data);

    Log::info("Hi from the FMI-LS-DAE interface.");
};

} // namespace FMI
