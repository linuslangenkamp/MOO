#include "gdop_strategies.h"
#include "gdop.h"

// TODO: add doxygen everywhere
// TODO: replace couts with format error_log() and log() prints
// TODO: add namespaces for every "major folder" :: NLP, Ipopt, OM, Base

namespace GDOP {

std::unique_ptr<Trajectory> DefaultConstantInitialization::operator()(GDOP& gdop) {
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

        // initial and final fixed: linear interpolation
        if (x0_opt && xf_opt) {
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

// no simulation available
std::unique_ptr<Trajectory> DefaultNoSimulation::operator()(GDOP& gdop, const ControlTrajectory& u, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
    std::cerr << "[Warning] No Simulation strategy set: returning nullptr." << std::endl;
    return nullptr;
}

// no simulation step available
std::unique_ptr<Trajectory> DefaultNoSimulationStep::operator()(GDOP& gdop, const ControlTrajectory& u, f64 start_time, f64 stop_time, f64* x_start_values) {
    std::cerr << "[Warning] No SimulationStep strategy set: returning nullptr." << std::endl;
    return nullptr;
}

// no mesh refinement available
std::unique_ptr<Trajectory> DefaultNoMeshRefinement::operator()(GDOP& gdop) {
    std::cerr << "[Warning] No MeshRefinement strategy set: returning nullptr." << std::endl;
    return nullptr;
}

// proper strategies

SimulationInitialization::SimulationInitialization(std::shared_ptr<Initialization> initialization,
                                                   std::shared_ptr<Simulation> simulation)
  : initialization(initialization), simulation(simulation) {}

std::unique_ptr<Trajectory> SimulationInitialization::operator()(GDOP& gdop) {
    auto simple_guess       = (*initialization)(gdop);                 // call simple, e.g. constant guess
    auto extracted_controls = simple_guess->copy_extract_controls();   // extract controls of the guess

    size_t x_size = simple_guess->x.size();
    FixedVector<f64> x0(x_size);
    for (size_t i = 0; i < x_size; i++) { x0[i] = simple_guess->x[i][0]; }

    auto simulated_guess    = (*simulation)(gdop, extracted_controls, gdop.mesh.node_count, // perform simulation using the controls and config
                                            0.0, gdop.mesh.tf, x0.raw());
    return simulated_guess;
}

// default strategy collection

Strategies Strategies::default_strategies() {
    return {std::make_shared<DefaultConstantInitialization>(),
            std::make_shared<DefaultNoSimulation>(),
            std::make_shared<DefaultNoSimulationStep>(),
            std::make_shared<DefaultNoMeshRefinement>()};
};

} // namespace GDOP
