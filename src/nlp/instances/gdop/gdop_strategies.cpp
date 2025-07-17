#include "gdop_strategies.h"
#include "gdop.h"

// TODO: add doxygen everywhere
// TODO: replace couts with format error_log() and log() prints

namespace GDOP {

// ==================== no-op strategies ====================

// no simulation available
std::unique_ptr<Trajectory> DefaultNoSimulation::operator()(const GDOP& gdop, const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
    LOG_ERROR("No Simulation strategy set: returning nullptr.");
    return nullptr;
}

// no simulation step available
std::unique_ptr<Trajectory> DefaultNoSimulationStep::operator()(const GDOP& gdop, const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) {
    LOG_ERROR("No SimulationStep strategy set: returning nullptr.");
    return nullptr;
}

// no mesh refinement available
std::unique_ptr<Trajectory> DefaultNoMeshRefinement::operator()(const GDOP& gdop) {
    LOG_ERROR("No MeshRefinement strategy set: returning nullptr.");
    return nullptr;
}

// no emitter
int DefaultNoEmitter::operator()(const GDOP& gdop, const Trajectory& trajectory) {
    LOG_ERROR("No Emitter strategy set: returning -1.");
    return -1;
}

// no verifier
bool DefaultNoVerifier::operator()(const GDOP& gdop, const Trajectory& trajectory) {
    LOG_ERROR("No Verifier strategy set: returning false.");
    return false;
}

// no scaling
std::shared_ptr<NLP::Scaling> DefaultNoScalingFactory::operator()(const GDOP& gdop) {
    return std::make_shared<NLP::NoScaling>(NLP::NoScaling());
};

// ==================== non no-op strategies ====================

// default initialization (is not really proper, but an implementation)
std::unique_ptr<Trajectory> DefaultConstantInitialization::operator()(const GDOP& gdop) {
    const auto& problem = gdop.problem;

    const size_t x_size = problem.x_bounds.size();
    const size_t u_size = problem.u_bounds.size();
    const size_t p_size = problem.p_bounds.size();

    // time vector with start and end times
    std::vector<f64> t = { 0.0, gdop.mesh.tf };

    std::vector<std::vector<f64>> x_guess(x_size);
    std::vector<std::vector<f64>> u_guess(u_size);
    std::vector<f64>              p_guess(p_size, 0.0);  // TODO: implement

    InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    // fill state guesses
    for (size_t x = 0; x < x_size; x++) {
        auto x0_opt = problem.x0_fixed[x];
        auto xf_opt = problem.xf_fixed[x];

        if (x0_opt && xf_opt) {
            // initial and final fixed: linear interpolation
            x_guess[x] = { *x0_opt, *xf_opt };
        } else if (x0_opt) {
            // initial fixed: constant at initial
            x_guess[x] = { *x0_opt, *x0_opt };
        } else if (xf_opt) {
            // final fixed: constant at final
            x_guess[x] = { *xf_opt, *xf_opt };
        } else {
            // nothing fixed: use midpoint of bounds or zero if unbounded
            double val = 0.0;
            if (problem.x_bounds[x].has_lower() && problem.x_bounds[x].has_upper()) {
                val = 0.5 * (problem.x_bounds[x].lb + problem.x_bounds[x].ub);
            }
            x_guess[x] = { val, val };
        }
    }

    // fill control guesses: use midpoint of bounds or zero
    for (size_t u = 0; u < u_size; u++) {
        double val = 0.0;
        if (problem.u_bounds[u].has_lower() && problem.u_bounds[u].has_upper()) {
            val = 0.5 * (problem.u_bounds[u].lb + problem.u_bounds[u].ub);
        }
        u_guess[u] = { val, val };
    }

    // fill parameter guesses: use midpoint of bounds or zero
    for (size_t p = 0; p < p_size; p++) {
        double val = 0.0;
        if (problem.p_bounds[p].has_lower() && problem.p_bounds[p].has_upper()) {
            val = 0.5 * (problem.p_bounds[p].lb + problem.p_bounds[p].ub);
        }
        p_guess[p] = val;
    }

    return std::make_unique<Trajectory>(Trajectory{ t, x_guess, u_guess, p_guess, interpolation });
}

// proper simulation-based initialization strategy
SimulationInitialization::SimulationInitialization(std::shared_ptr<Initialization> initialization,
                                                   std::shared_ptr<Simulation> simulation)
  : initialization(initialization), simulation(simulation) {}

std::unique_ptr<Trajectory> SimulationInitialization::operator()(const GDOP& gdop) {
    auto simple_guess       = (*initialization)(gdop);                                           // call simple, e.g. constant guess
    auto extracted_controls = simple_guess->copy_extract_controls();                             // extract controls from the guess
    auto exctracted_x0      = simple_guess->extract_initial_states();                            // extract x(t_0) from the guess
    auto simulated_guess    = (*simulation)(gdop, extracted_controls, gdop.mesh.node_count, 0.0, // perform simulation using the controls and gdop config
                                            gdop.mesh.tf, exctracted_x0.raw());
    return simulated_guess;
}

// csv emit
CSVEmitter::CSVEmitter(std::string filename) : filename(filename) {}

int CSVEmitter::operator()(const GDOP& gdop, const Trajectory& trajectory) { return trajectory.to_csv(filename); }

// simulation-based verification
SimulationVerifier::SimulationVerifier(std::shared_ptr<Simulation> simulation,
                                       Linalg::Norm norm,
                                       FixedVector<f64>&& tolerances)
    : simulation(simulation), norm(norm), tolerances(std::move(tolerances)) {}

bool SimulationVerifier::operator()(const GDOP& gdop, const Trajectory& trajectory) {
    auto extracted_controls = trajectory.copy_extract_controls();   // extract controls from the trajectory
    auto exctracted_x0      = trajectory.extract_initial_states();  // extract x(t_0) from the trajectory

    // perform simulation using the controls, gdop config and a very high number of nodes
    int  high_node_count    = 10 * gdop.mesh.node_count;
    auto simulation_result  = (*simulation)(gdop, extracted_controls, high_node_count,
                                             0.0, gdop.mesh.tf, exctracted_x0.raw());

    // result of high resolution simulation is interpolated onto lower resolution mesh
    auto interpolated_sim   = simulation_result->interpolate_onto_mesh(gdop.mesh, gdop.collocation);

    // calculate errors for each state in given norm (between provided and simulated states)
    auto errors             = trajectory.state_errors(interpolated_sim, norm);

    bool is_valid = true;

    LOG_START_MODULE("Simulation-Based Verification");

    FixedTableFormat<4> fmt = {{7, 12, 12, 4}};

    LOG_HEADER(fmt, "State", "Error", "Tolerance", "Pass");
    LOG_DASHES(fmt);

    for (size_t x_idx = 0; x_idx < trajectory.x.size(); x_idx++) {
        f64 tol = tolerances[x_idx];
        f64 err = errors[x_idx];
        bool pass = (err <= tol);

        LOG_ROW(fmt,
            fmt::format("x[{}]", x_idx),
            fmt::format("{:.3e}", err),
            fmt::format("{:.3e}", tol),
            pass ? "" : "x");

        if (!pass) {
            is_valid = false;
        }
    }

    LOG_DASHES(fmt);

    if (is_valid) {
        LOG_SUCCESS("All state errors within tolerances.");
    } else {
        LOG_WARNING("One or more state errors exceeded tolerances.");
    }

    return is_valid;
}

// default strategy collection
Strategies Strategies::default_strategies() {
    Strategies strategy;
    strategy.initialization  = std::make_shared<DefaultConstantInitialization>();
    strategy.simulation      = std::make_shared<DefaultNoSimulation>();
    strategy.simulation_step = std::make_shared<DefaultNoSimulationStep>();
    strategy.mesh_refinement = std::make_shared<DefaultNoMeshRefinement>();
    strategy.emitter         = std::make_shared<DefaultNoEmitter>();
    strategy.verifier        = std::make_shared<DefaultNoVerifier>();
    strategy.scaling_factory = std::make_shared<DefaultNoScalingFactory>();
    return strategy;
};

} // namespace GDOP
