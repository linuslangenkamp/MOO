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

    // there should be some better / a single variant:
    //     1. one big loop over all variables => split into disjoint sets we need, i.e. x, derx, u, p, alg, residuals
    //     2. immediately get their number and value references

    {
        int nVars = fmi3_getNumberOfVariables(priv->fmu);
        for (int idx = 1; idx < nVars + 1; idx++) {
            fmi3VariableHandle* var = fmi3_getVariableByIndex(priv->fmu, idx);
            int derx_ref = fmi3_getVariableValueReference(var);
            int x_ref = fmi3_getVariableDerivativeIndex(fmi3_getVariableByValueReference(priv->fmu, derx_ref));

            if (fmi3_getVariableCausality(var) == fmi3CausalityInput) {
                priv->input_vrefs.push_back(derx_ref);
            } else if (x_ref != 0) {
                priv->state_vrefs.push_back(x_ref);
                priv->deriv_vrefs.push_back(derx_ref);
            }
        }
    }

    {
        int count = fmiLsDae_getNumberOfAlgebraicVariables(priv->fmu);
        priv->alg_vrefs.reserve(count);
        for (int idx = 1; idx < count + 1; idx++) {
            fmiLsDaeAlgebraicVariableHandle* h = fmiLsDae_getAlgebraicVariableByIndex(priv->fmu, idx - 1);
            priv->alg_vrefs.push_back(fmiLsDae_getAlgebraicVariableValueReference(h));
        }
    }

    {
        int count = fmiLsDae_getNumberOfResiduals(priv->fmu);
        priv->res_vrefs.reserve(count);
        for (int idx = 1; idx < count + 1; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getResidualByIndex(priv->fmu, idx - 1);
            priv->res_vrefs.push_back(fmiLsDae_getValueReference(h));
        }
    }

    {
        int count = fmiLsDae_getNumberOfOutputs(priv->fmu);
        priv->output_vrefs.reserve(count);
        for (int idx = 1; idx < count + 1; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getOutputByIndex(priv->fmu, idx - 1);
            priv->output_vrefs.reserve(fmiLsDae_getValueReference(h));
        }
    }

    Log::info("FMI version  : {}", static_cast<int>(fmi4c_getFmiVersion(priv->fmu)));
    Log::info("Model name   : {}", fmi3_modelName(priv->fmu));
    Log::info("  #x      = {}  (states)",      priv->state_vrefs.size());
    Log::info("  #der(x) = {}  (derivatives)", priv->deriv_vrefs.size());
    Log::info("  #u      = {}  (inputs)",      priv->input_vrefs.size());
    Log::info("  #z      = {}  (algebraics)",  priv->alg_vrefs.size());
    Log::info("  #g      = {}  (residuals)",   priv->res_vrefs.size());
    Log::info("  #y      = {}  (outputs)",     priv->output_vrefs.size());

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

GDOP::ProblemConstants create_problem_constants(FMIData& fmi_data)
{
    auto x_bounds = FixedVector<Bounds>(fmi_data.priv->state_vrefs.size());
    auto u_bounds = FixedVector<Bounds>(fmi_data.priv->input_vrefs.size() + fmi_data.priv->alg_vrefs.size());
    auto p_bounds = FixedVector<Bounds>(0);
    auto T_bounds = std::array<Bounds, 2>{};
    auto g_bounds = FixedVector<Bounds>(fmi_data.priv->res_vrefs.size());
    for (auto i = 0; i < g_bounds.int_size(); i++)
    {
        g_bounds[i].lb = 0.0;
        g_bounds[i].ub = 0.0;
    }

    auto r_bounds = FixedVector<Bounds>(0);

    auto xu0_fixed = FixedVector<std::optional<f64>>(x_bounds.size() + u_bounds.size());
    auto xuf_fixed = FixedVector<std::optional<f64>>(x_bounds.size() + u_bounds.size());
    auto T_fixed = std::array<std::optional<f64>, 2>{};

    Log::info("=== Problem Bounds ===");

    Log::info("State bounds (x):");
    for (size_t i = 0; i < x_bounds.size(); ++i) {
        const auto& b = x_bounds[i];
        Log::info("  x[{}]: lb = {}, ub = {}", i, b.lb, b.ub);
    }

    Log::info("Input bounds (u):");
    size_t nu = fmi_data.priv->input_vrefs.size();

    for (size_t i = 0; i < nu; ++i) {
        const auto& b = u_bounds[i];
        Log::info("  u[{}]: lb = {}, ub = {}", i, b.lb, b.ub);
    }

    Log::info("Algebraic bounds (z):");
    for (size_t i = nu; i < u_bounds.size(); ++i) {
        const auto& b = u_bounds[i];
        Log::info("  z[{}]: lb = {}, ub = {}", i - nu, b.lb, b.ub);
    }

    Log::info("Parameter bounds (p):");
    for (size_t i = 0; i < p_bounds.size(); ++i) {
        const auto& b = p_bounds[i];
        Log::info("  p[{}]: lb = {}, ub = {}", i, b.lb, b.ub);
    }

    Log::info("Residual bounds (g):");
    for (size_t i = 0; i < g_bounds.size(); ++i) {
        const auto& b = g_bounds[i];
        Log::info("  g[{}]: lb = {}, ub = {}", i, b.lb, b.ub);
    }

    Log::info("================================");

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
}

FullSweep::FullSweep(GDOP::FullSweepLayout&& layout_in,
                     const GDOP::ProblemConstants& pc,
                     FMIData& fmi_data)
    : GDOP::FullSweep(std::move(layout_in), pc),
      fmi_data(fmi_data)
{}

BoundarySweep::BoundarySweep(GDOP::BoundarySweepLayout&& layout_in,
                             const GDOP::ProblemConstants& pc,
                             FMIData& fmi_data)
    : GDOP::BoundarySweep(std::move(layout_in), pc),
      fmi_data(fmi_data)
{}


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

void BoundarySweep::callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, f64 t0, f64 tf)
{
    // not yet implemented
}

void BoundarySweep::callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, f64 t0, f64 tf)
{
    // not yet implemented
}

void BoundarySweep::callback_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, f64 t0, f64 tf, const f64 mayer_factor, const f64* lambda)
{
    // not yet implemented
}

// we just fill the buffer in order (L, f1, ..., fn, g1, ..., gn) <- buffer layout
void layout_lfg_init_eval(GDOP::FullSweepLayout& layout, FMIData& fmi_data, GDOP::ProblemConstants& pc)
{
    if (layout.L) {
        layout.L->buf_index = 0;
    }

    int offset = (layout.L ? 1 : 0);
    for (int f_idx = 0; f_idx < pc.x_size; f_idx++) {
        layout.f[f_idx].buf_index = offset + f_idx;
    }

    offset += pc.x_size;
    for (int g_idx = 0; g_idx < pc.g_size; g_idx++) {
        layout.g[g_idx].buf_index = offset + g_idx;
    }
}

void layout_lfg_init_jac(GDOP::FullSweepLayout& layout, FMIData& fmi_data, GDOP::ProblemConstants& pc)
{

}

GDOP::FullSweepLayout create_fullsweep_layout(FMIData& fmi_data, GDOP::ProblemConstants& pc) {
    auto layout_lfg = GDOP::FullSweepLayout(pc.has_lagrange, pc.x_size, pc.g_size);

    layout_lfg_init_eval(layout_lfg, fmi_data, pc);
    layout_lfg_init_jac(layout_lfg, fmi_data, pc);

    return layout_lfg;
}

GDOP::Problem create_gdop_problem(FMIData& fmi_data)
{
    auto pc = std::make_unique<GDOP::ProblemConstants>(create_problem_constants(fmi_data));
    auto fs = std::make_unique<FullSweep>(create_fullsweep_layout(fmi_data, *pc), *pc, fmi_data); 
    auto empty_bs = std::make_unique<BoundarySweep>(GDOP::BoundarySweepLayout(false, 0), *pc, fmi_data);

    return GDOP::Problem(std::move(fs), std::move(empty_bs), std::move(pc), nullptr);
}

Problem::Problem(FMIData& fmi_data)
    : GDOP::Problem(create_gdop_problem(fmi_data))
{}

void main_fmi(const char* path, const char* modelname)
{
    Log::info("Hi from the FMI-LS-DAE interface.");
    Log::info("Path {}, Name {}", path, modelname);

    auto fmi_data = FMIData(path, modelname);
    auto problem = Problem(fmi_data);
}

} // namespace FMI
