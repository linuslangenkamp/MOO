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

#include <nlp/instances/gdop/orchestrator.h>
#include <nlp/solvers/ipopt/solver.h>
#include <nlp/instances/gdop/strategies.h>

#include <fmi4c.h>
#include <fmi4c_lsdae.h>

namespace FMI {

enum class VarCategory {
    State,      // x
    Input,      // u
    Algebraic,  // z
    Parameter,  // p
    Unknown
};

struct VarIndex {
    VarCategory category;
    int local_index;  // index inside strictly x, u, or z
    int semi_local;   // index inside x or concatenated (u+z)
    int xu_index;     // index inside fully concatenated (x, u+z)
};

struct Dependency {
    fmi3ValueReference vref;
    VarIndex var;
};

struct FMIData_priv {
    fmuHandle* fmu;
    fmi3InstanceHandle* instance;

    std::vector<fmi3ValueReference> lagrange_vrefs;
    std::vector<fmi3ValueReference> state_vrefs;
    std::vector<fmi3ValueReference> deriv_vrefs;
    std::vector<fmi3ValueReference> input_vrefs;
    std::vector<fmi3ValueReference> alg_vrefs;
    std::vector<fmi3ValueReference> res_vrefs;
    std::vector<fmi3ValueReference> output_vrefs;
    std::vector<fmi3ValueReference> parameter_vrefs;

    std::unordered_map<fmi3ValueReference, VarIndex> vref_map;
    std::vector<std::vector<Dependency>> lagr_dependencies;
    std::vector<std::vector<Dependency>> deriv_dependencies;
    std::vector<std::vector<Dependency>> res_dependencies;

    std::vector<fmi3ValueReference> Lfg_vrefs; // value references of the (L, f1, ..., fn, g1, ..., gm) vector
    std::vector<fmi3ValueReference> xuz_vrefs; // value references of the (x1, ..., xn, u1, ..., uk, z1, ..., zl) vector

    FMIData_priv(const char* path, const char* modelname);
};

FMIData_priv::FMIData_priv(const char* path, const char* modelname)
    : fmu(fmi4c_loadFmu(path, modelname)),
      instance(fmi3_instantiateModelExchange(fmu, true, true, nullptr, nullptr))
{}

void FMIData::print() {
    Log::info("\n=== FMI Model Metadata ===");
    Log::info("FMI version  : {}", static_cast<int>(fmi4c_getFmiVersion(priv->fmu)));
    Log::info("Model name   : {}", fmi3_modelName(priv->fmu));

    FixedTableFormat<4> summary_fmt = {{15, 8, 25, 10}, {Align::Center, Align::Center, Align::Center, Align::Center}};
    Log::start_module(summary_fmt, "Variable Counts & Preview");
    Log::row(summary_fmt, "Category", "Count", "First 2 VRefs", "Symbol");
    Log::dashes(summary_fmt);

    auto get_preview = [](const std::vector<fmi3ValueReference>& vrefs) -> std::string {
        if (vrefs.empty()) return "-";

        std::string s;
        size_t take = std::min<size_t>(vrefs.size(), 2);
        for (size_t i = 0; i < take; ++i) {
            s += std::to_string(vrefs[i]) + (i < take - 1 ? ", " : "");
        }
        if (vrefs.size() > 2) s += "...";
        return s;
    };

    Log::row(summary_fmt, "States",      priv->state_vrefs.size(),  get_preview(priv->state_vrefs),  "x");
    Log::row(summary_fmt, "Derivatives", priv->deriv_vrefs.size(),  get_preview(priv->deriv_vrefs),  "der(x)");
    Log::row(summary_fmt, "Inputs",      priv->input_vrefs.size(),  get_preview(priv->input_vrefs),  "u");
    Log::row(summary_fmt, "Algebraics",  priv->alg_vrefs.size(),   get_preview(priv->alg_vrefs),    "z");
    Log::row(summary_fmt, "Residuals",   priv->res_vrefs.size(),    get_preview(priv->res_vrefs),    "g");
    Log::row(summary_fmt, "Outputs",     priv->output_vrefs.size(), get_preview(priv->output_vrefs), "y");
    Log::row(summary_fmt, "Parameter",   priv->parameter_vrefs.size(), get_preview(priv->parameter_vrefs), "p");
    Log::dashes(summary_fmt);

    FixedTableFormat<5> map_fmt = {{15, 10, 10, 12, 12}, {Align::Center, Align::Center, Align::Center, Align::Center, Align::Center}};
    Log::start_module(map_fmt, "Index Mapping");
    Log::row(map_fmt, "Category", "VRef", "Local (L)", "Semi-Loc (SL)", "Global (G)");
    Log::dashes(map_fmt);

    auto print_map_rows = [&](const std::string& cat, const std::vector<fmi3ValueReference>& vrefs) {
        for (auto vref : vrefs) {
            auto it = priv->vref_map.find(vref);
            if (it != priv->vref_map.end()) {
                Log::row(map_fmt, cat, vref, it->second.local_index, it->second.semi_local, it->second.xu_index);
            }
        }
    };

    print_map_rows("State (x)", priv->state_vrefs);
    print_map_rows("Input (u)", priv->input_vrefs);
    print_map_rows("Algebraic (z)", priv->alg_vrefs);
    Log::dashes(map_fmt);

    FixedTableFormat<3> dep_fmt = {{20, 15, 45}, {Align::Center, Align::Center, Align::Center}};
    Log::start_module(dep_fmt, "System Sparsity (Global Jacobian Indices)");
    Log::row(dep_fmt, "Equation", "VRef", "Global Column Indices (G)");
    Log::dashes(dep_fmt);

    auto print_deps = [&](const std::string& prefix,
                          const std::vector<fmi3ValueReference>& vrefs,
                          const std::vector<std::vector<Dependency>>& deps_list) {
        for (size_t i = 0; i < vrefs.size(); ++i) {
            std::string indices_str;
            for (const auto& d : deps_list[i]) {
                if (d.var.category == VarCategory::Unknown) {
                    indices_str += "??, ";
                } else {
                    indices_str += std::to_string(d.var.xu_index) + ", ";
                }
            }

            if (!indices_str.empty()) {
                indices_str.erase(indices_str.length() - 2);
            }

            Log::row(dep_fmt, fmt::format("{}[{}]", prefix, i), i, indices_str.empty() ? "none" : indices_str);
        }
    };

    print_deps("L", priv->lagrange_vrefs, priv->lagr_dependencies);
    Log::dashes(dep_fmt);
    print_deps("der(x)", priv->deriv_vrefs, priv->deriv_dependencies);
    Log::dashes(dep_fmt);
    print_deps("Residual g", priv->res_vrefs, priv->res_dependencies);
    Log::dashes(dep_fmt);
}

FMIData::FMIData(FMISettings& settings)
    : priv(std::make_unique<FMIData_priv>(settings.path, settings.modelname))
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

    int states = fmiLsDae_getNumberOfContinuousStateDerivatives(priv->fmu);
    for (int s = 0; s < states; s++)
    {
        auto structure = fmiLsDae_getContinuousStateDerivativeByIndex(priv->fmu, s);
        auto der_ref = fmiLsDae_getValueReference(structure);
        auto x_ref = fmi3_getVariableDerivativeIndex(fmi3_getVariableByValueReference(priv->fmu, der_ref));
        priv->state_vrefs.push_back(x_ref);
        priv->deriv_vrefs.push_back(der_ref);
    }

    {
        int nVars = fmi3_getNumberOfVariables(priv->fmu);
        for (int idx = 1; idx < nVars + 1; idx++) {
            fmi3VariableHandle* var = fmi3_getVariableByIndex(priv->fmu, idx);
            int ref = fmi3_getVariableValueReference(var);

            if (fmi3_getVariableCausality(var) == fmi3CausalityInput) {
                priv->input_vrefs.push_back(ref);
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
            priv->output_vrefs.push_back(fmiLsDae_getValueReference(h));
        }
    }

    for (auto& ref : settings.parameter_vrefs)
    {
        auto handle = fmi3_getVariableByValueReference(priv->fmu, ref);
        auto variability = fmi3_getVariableVariability(handle);
        if (variability == fmi3VariabilityTunable) {
            priv->parameter_vrefs.push_back(ref);
        }
        else {
            Log::error("Provided parameter value reference {} is not tunable.", ref);
            throw;
        }
    }

    if (settings.lagrange_vref != nullptr)
    {
        auto it = std::find(priv->output_vrefs.begin(), priv->output_vrefs.end(), static_cast<fmi3ValueReference>(*settings.lagrange_vref));
        if (it == priv->output_vrefs.end()) {
            Log::error("Provided Lagrange term value reference {} does not exist.", *settings.lagrange_vref);
            throw;
        }
        priv->lagrange_vrefs.push_back(*settings.lagrange_vref);
    }

    // build the Variable Index Map (Concatenated xu = [x, u, z])
    {
        int xu_offset = 0;
        int u_size = static_cast<int>(priv->input_vrefs.size());

        for (size_t i = 0; i < priv->state_vrefs.size(); ++i) {
            priv->vref_map[priv->state_vrefs[i]] = {
                VarCategory::State, static_cast<int>(i), static_cast<int>(i), xu_offset++
            };
        }

        for (size_t i = 0; i < priv->input_vrefs.size(); ++i) {
            priv->vref_map[priv->input_vrefs[i]] = {
                VarCategory::Input, static_cast<int>(i), static_cast<int>(i), xu_offset++
            };
        }

        for (size_t i = 0; i < priv->alg_vrefs.size(); ++i) {
            priv->vref_map[priv->alg_vrefs[i]] = {
                VarCategory::Algebraic, static_cast<int>(i), u_size + static_cast<int>(i), xu_offset++
            };
        }

        for (size_t i = 0; i < priv->parameter_vrefs.size(); ++i) {
            priv->vref_map[priv->parameter_vrefs[i]] = {
                VarCategory::Parameter, static_cast<int>(i), static_cast<int>(i), static_cast<int>(i)
            };
        }
    }

    // build dependency structure

    if (priv->lagrange_vrefs.size() != 0)
    {
        int count = static_cast<int>(priv->output_vrefs.size());
        priv->lagr_dependencies.reserve(1);

        for (int idx = 0; idx < count; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getOutputByIndex(priv->fmu, idx);
            if (h->valueReference != static_cast<fmi3ValueReference>(priv->lagrange_vrefs[0])) continue;

            int n_deps = fmiLsDae_getNumberOfDependencies(h);
            std::vector<fmi3ValueReference> raw_deps(n_deps);
            std::vector<Dependency> deps;
            deps.reserve(n_deps);

            if (n_deps > 0) {
                fmiLsDae_getDependencies(h, raw_deps.data(), n_deps);
                for (auto vref : raw_deps) {
                    auto it = priv->vref_map.find(vref);
                    if (it != priv->vref_map.end()) {
                        deps.push_back({vref, it->second});
                    } else {
                        // unmapped e.g. parameters?
                        deps.push_back({vref, {VarCategory::Unknown, -1, -1, -1}});
                    }
                }
            }
            priv->lagr_dependencies.push_back(std::move(deps));
        }
    }

    {
        int count = static_cast<int>(priv->state_vrefs.size());
        priv->deriv_dependencies.reserve(count);

        for (int idx = 0; idx < count; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getContinuousStateDerivativeByIndex(priv->fmu, idx);

            int n_deps = (h == nullptr) ? 0 : fmiLsDae_getNumberOfDependencies(h);
            std::vector<fmi3ValueReference> raw_deps(n_deps);
            std::vector<Dependency> deps;
            deps.reserve(n_deps);

            if (n_deps > 0) {
                fmiLsDae_getDependencies(h, raw_deps.data(), n_deps);
                for (auto vref : raw_deps) {
                    auto it = priv->vref_map.find(vref);
                    if (it != priv->vref_map.end()) {
                        deps.push_back({vref, it->second});
                    } else {
                        // unmapped e.g. parameters?
                        deps.push_back({vref, {VarCategory::Unknown, -1, -1, -1}});
                    }
                }
            }
            priv->deriv_dependencies.push_back(std::move(deps));
        }
    }

    {
        int count = static_cast<int>(priv->res_vrefs.size());
        priv->res_dependencies.reserve(count);

        for (int idx = 0; idx < count; idx++) {
            fmiLsDaeModelStructureHandle* h = fmiLsDae_getResidualByIndex(priv->fmu, idx);

            int n_deps = fmiLsDae_getNumberOfDependencies(h);
            std::vector<fmi3ValueReference> raw_deps(n_deps);
            std::vector<Dependency> deps;
            deps.reserve(n_deps);

            if (n_deps > 0) {
                fmiLsDae_getDependencies(h, raw_deps.data(), n_deps);
                for (auto vref : raw_deps) {
                    auto it = priv->vref_map.find(vref);
                    if (it != priv->vref_map.end()) {
                        deps.push_back({vref, it->second});
                    } else {
                        // unmapped e.g. parameters?
                        deps.push_back({vref, {VarCategory::Unknown, -1, -1, -1}});
                    }
                }
            }
            priv->res_dependencies.push_back(std::move(deps));
        }
    }

    // create vectors of variable references as in our sorted MOO structures

    if (settings.lagrange_vref != nullptr)
    {
        priv->Lfg_vrefs.push_back(priv->lagrange_vrefs[0]);
    }

    for (auto f : priv->deriv_vrefs) {
        priv->Lfg_vrefs.push_back(f);
    }

    for (auto g : priv->res_vrefs) {
        priv->Lfg_vrefs.push_back(g);
    }

    for (auto x : priv->state_vrefs) {
        priv->xuz_vrefs.push_back(x);
    }

    for (auto u : priv->input_vrefs) {
        priv->xuz_vrefs.push_back(u);
    }

    for (auto z : priv->alg_vrefs) {
        priv->xuz_vrefs.push_back(z);
    }

    print();
}

void FMIData::initialize(f64 t_start, f64 t_stop) {
    auto* inst = priv->instance;

    fmi3Status status = fmi3_enterInitializationMode(
        inst,
        fmi3False, 0.0,     // toleranceDefined, tolerance
        t_start,            // startTime
        fmi3True, t_stop    // stopTimeDefined, stopTime
    );

    if (status != fmi3OK) {
        throw std::runtime_error("FMU: Failed to enter initialization mode");
    }

    status = fmi3_exitInitializationMode(inst);
    if (status != fmi3OK) {
        throw std::runtime_error("FMU: Failed to exit initialization mode");
    }

    /* TODO: What do we need to do here?

    status = fmi3_enterEventMode(inst);
    if (status != fmi3OK) {
        throw std::runtime_error("FMU: Failed to enter event mode");
    }

    status = fmi3_enterStepMode(inst);
    if (status != fmi3OK) {
        throw std::runtime_error("FMU: Failed to enter step time mode");
    }

    */

    Log::info("FMU is now in Continuous-Time Mode (ready for simulation).");
}

std::vector<f64> FMIData::get_initial_states() {
    if (priv->state_vrefs.empty()) return {};

    std::vector<f64> states(priv->state_vrefs.size());

    fmi3Status status = fmi3_getFloat64(priv->instance,
                                        priv->state_vrefs.data(),
                                        priv->state_vrefs.size(),
                                        states.data(),
                                        states.size());

    if (status != fmi3OK) {
        Log::error("Failed to retrieve initial states from FMU");
    }

    return states;
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

    // get and fix start values for state variables
    auto x0 = fmi_data.get_initial_states();

    for (size_t i = 0; i < x0.size(); ++i) {
        xu0_fixed[i] = x0[i];
    }

    for (size_t i = 0; i < fmi_data.priv->input_vrefs.size(); ++i) {
        u_bounds[i].lb = -1;
        u_bounds[i].ub = 1;
    }


    Log::info("\n=== Start values for States x(0.0) ===");
    for (size_t i = 0; i < x0.size(); ++i) {
        Log::info("  x[{}](0.0) = {}", i, x0[i]);
    }

    Log::info("\n=== Problem Bounds ===");

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
        /* tf */        0.5,
        /* intervals */ 15,
        /* stages */    3,
        /* type */      MeshType::Physical
    );

    T_fixed[0] = mesh->t0;
    T_fixed[1] = mesh->tf;

    return GDOP::ProblemConstants(
        /* mayer */ false,
        /* lagrange */ !fmi_data.priv->lagrange_vrefs.empty(),
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

// evaluate (L, f, g)(x, u, w, p, time)
void FMIData::eval_point_lfg(const f64* xu, const f64* p, f64 time, f64* out)
{
    fmi3_setTime(priv->instance, time);

    fmi3_setFloat64(priv->instance,
                    priv->xuz_vrefs.data(),
                    priv->xuz_vrefs.size(),
                    xu,
                    priv->xuz_vrefs.size());

    fmi3_getContinuousStateDerivatives(priv->instance, out, priv->state_vrefs.size());

    fmi3_getFloat64(priv->instance,
                    priv->Lfg_vrefs.data(),
                    priv->Lfg_vrefs.size(),
                    out,
                    priv->Lfg_vrefs.size());
}

// evaluate d(L, f, g) / d(x, u, z, p) (x, u, z, p, time)
void FMIData::jac_point_lfg(const f64* xu, const f64* p, f64 time, f64* out)
{
    fmi3InstanceHandle* inst = priv->instance;

    fmi3_setTime(inst, time);
    fmi3_setFloat64(inst,
                    priv->xuz_vrefs.data(),
                    priv->xuz_vrefs.size(),
                    xu,
                    priv->xuz_vrefs.size());

    int out_idx = 0;

    auto get_partial = [&](fmi3ValueReference value_vref, fmi3ValueReference independent_vref) -> f64 {
        f64 delta = 1.0;
        f64 partial_derivative = 0.0;

        fmi3Status status = fmi3_getDirectionalDerivative(inst,
                                                          &value_vref, 1,
                                                          &independent_vref, 1,
                                                          &delta, 1,
                                                          &partial_derivative, 1);

        return (status == fmi3OK) ? partial_derivative : 0.0;
    };

    for (size_t i = 0; i < priv->lagrange_vrefs.size(); ++i) {
        fmi3ValueReference L_vref = priv->lagrange_vrefs[i];
        for (const auto& dep : priv->lagr_dependencies[i]) {
            out[out_idx++] = get_partial(L_vref, dep.vref);
        }
    }

    for (size_t i = 0; i < priv->deriv_vrefs.size(); ++i) {
        fmi3ValueReference f_vref = priv->deriv_vrefs[i];
        for (const auto& dep : priv->deriv_dependencies[i]) {
            out[out_idx++] = get_partial(f_vref, dep.vref);
        }
    }

    for (size_t i = 0; i < priv->res_vrefs.size(); ++i) {
        fmi3ValueReference g_vref = priv->res_vrefs[i];
        for (const auto& dep : priv->res_dependencies[i]) {
            out[out_idx++] = get_partial(g_vref, dep.vref);
        }
    }
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
            fmi_data.jac_point_lfg(xu_ij, p, pc.mesh->t[i][j], jac_buf_ij);
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
    int jac_offset = 0;

    if (layout.L) {
        auto& L_jac = layout.L->jac;
        const auto& fmi_deps = fmi_data.priv->lagr_dependencies[0];

        for (const auto& dep : fmi_deps) {
            if (dep.var.category == VarCategory::State) {
                L_jac.dx.push_back({dep.var.semi_local, jac_offset++});
            }
            else if (dep.var.category == VarCategory::Input || dep.var.category == VarCategory::Algebraic) {
                L_jac.du.push_back({dep.var.semi_local, jac_offset++});
            }
            else if (dep.var.category == VarCategory::Parameter) {
                L_jac.dp.push_back({dep.var.local_index, jac_offset++});
            }
        }
    }

    for (int i = 0; i < pc.x_size; i++) {
        auto& current_f_jac = layout.f[i].jac;
        const auto& fmi_deps = fmi_data.priv->deriv_dependencies[i];

        for (const auto& dep : fmi_deps) {
            if (dep.var.category == VarCategory::State) {
                current_f_jac.dx.push_back({dep.var.semi_local, jac_offset++});
            }
            else if (dep.var.category == VarCategory::Input || dep.var.category == VarCategory::Algebraic) {
                current_f_jac.du.push_back({dep.var.semi_local, jac_offset++});
            }
            else if (dep.var.category == VarCategory::Parameter) {
                current_f_jac.dp.push_back({dep.var.local_index, jac_offset++});
            }
        }
    }

    for (int i = 0; i < pc.g_size; i++) {
        auto& current_g_jac = layout.g[i].jac;
        const auto& fmi_deps = fmi_data.priv->res_dependencies[i];

        for (const auto& dep : fmi_deps) {
            if (dep.var.category == VarCategory::State) {
                current_g_jac.dx.push_back({dep.var.semi_local, jac_offset++});
            }
            else if (dep.var.category == VarCategory::Input || dep.var.category == VarCategory::Algebraic) {
                current_g_jac.du.push_back({dep.var.semi_local, jac_offset++});
            }
            else if (dep.var.category == VarCategory::Parameter) {
                current_g_jac.dp.push_back({dep.var.local_index, jac_offset++});
            }
        }
    }

    Log::info("Jacobian sparsity initialized. Total NNZ per node: {}", jac_offset);
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

    // test sparsity of the fullsweep
    fs->print_jacobian_sparsity_pattern();

    return GDOP::Problem(std::move(fs), std::move(empty_bs), std::move(pc), nullptr);
}

Problem::Problem(FMIData& fmi_data)
    : GDOP::Problem(create_gdop_problem(fmi_data))
{}

void main_fmi(FMISettings& settings) {
    Log::info("Hi from the FMI-LS-DAE interface.");
    Log::info("Path {}, Name {}", settings.path, settings.modelname);

    auto fmi_data = FMIData(settings);
    fmi_data.initialize(0.0, 0.0);
    auto problem = Problem(fmi_data);

    auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    strategies->emitter = std::make_shared<GDOP::CSVEmitter>("optimal_solution.csv", false);

    auto gdop = GDOP::GDOP(problem);

    auto nlp_solver_settings = NLP::NLPSolverSettings(0, nullptr);
    nlp_solver_settings.set(NLP::Option::Hessian, NLP::HessianOption::LBFGS);
    nlp_solver_settings.set(NLP::Option::Jacobian, NLP::JacobianOption::Exact);
    nlp_solver_settings.set(NLP::Option::Gradient, NLP::GradientOption::Exact);
    // nlp_solver_settings.set(NLP::Option::IpoptDerivativeTest, true);
    nlp_solver_settings.set(NLP::Option::Tolerance, 1e-8);
    nlp_solver_settings.print();

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_settings);

    auto orchestrator = GDOP::MeshRefinementOrchestrator(gdop, std::move(strategies), ipopt_solver);

    orchestrator.optimize();

    TimingTree::instance().print_tree_table();
}

} // namespace FMI
