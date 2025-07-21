#include "gdop_strategies.h"
#include "gdop.h"

// TODO: add doxygen everywhere
// TODO: replace couts with format error_log() and log() prints

namespace GDOP {

// ==================== no-op strategies ====================

// no simulation available
std::unique_ptr<Trajectory> DefaultNoSimulation::operator()(const GDOP& gdop, const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
    LOG_WARNING("No Simulation strategy set: returning nullptr.");
    return nullptr;
}

// no simulation step available
std::unique_ptr<Trajectory> DefaultNoSimulationStep::operator()(const GDOP& gdop, const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) {
    LOG_WARNING("No SimulationStep strategy set: returning nullptr.");
    return nullptr;
}

// no mesh refinement available
void DefaultNoMeshRefinement::reinit(const GDOP& gdop) {}

std::unique_ptr<MeshUpdate> DefaultNoMeshRefinement::operator()(const GDOP& gdop) {
    LOG_WARNING("No MeshRefinement strategy set: returning nullptr.");
    return nullptr;
}

// no emitter
int DefaultNoEmitter::operator()(const GDOP& gdop, const Trajectory& trajectory) {
    LOG_WARNING("No Emitter strategy set: returning -1.");
    return -1;
}

// no verifier
bool DefaultNoVerifier::operator()(const GDOP& gdop, const Trajectory& trajectory) {
    LOG_WARNING("No Verifier strategy set: returning false.");
    return false;
}

// no scaling
std::shared_ptr<NLP::Scaling> DefaultNoScalingFactory::operator()(const GDOP& gdop) {
    LOG_WARNING("No ScalingFactory strategy set: fallback to DefaultNoScalingFactory.");
    return std::make_shared<NLP::NoScaling>(NLP::NoScaling());
};

// ==================== non no-op strategies ====================

// default initialization (is not really proper, but an implementation)
std::unique_ptr<Trajectory> DefaultConstantInitialization::operator()(const GDOP& gdop) {
    LOG_WARNING("No Initialization strategy set: fallback to DefaultConstantInitialization.");

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

// interpolate trajectory to new mesh with simple linear interpolation
std::unique_ptr<Trajectory> DefaultLinearInterpolation::operator()(const GDOP& gdop,
                                                                   const Mesh& new_mesh,
                                                                   const Trajectory& trajectory) {
    return std::make_unique<Trajectory>(trajectory.interpolate_onto_mesh_linear(new_mesh, gdop.collocation));
}

// proper simulation-based initialization strategy
SimulationInitialization::SimulationInitialization(std::shared_ptr<Initialization> initialization,
                                                   std::shared_ptr<Simulation> simulation)
  : initialization(initialization), simulation(simulation) {}

std::unique_ptr<Trajectory> SimulationInitialization::operator()(const GDOP& gdop) {
    auto simple_guess       = (*initialization)(gdop);                                             // call simple, e.g. constant guess
    auto extracted_controls = simple_guess->copy_extract_controls();                               // extract controls from the guess
    auto exctracted_x0      = simple_guess->extract_initial_states();                              // extract x(t_0) from the guess
    auto simulated_guess    = (*simulation)(gdop, extracted_controls, gdop.mesh.node_count, 0.0,   // perform simulation using the controls and gdop config
                                            gdop.mesh.tf, exctracted_x0.raw());
    auto interpolated_sim   = simulated_guess->interpolate_onto_mesh(gdop.mesh, gdop.collocation); // interpolate simulation to current mesh + collocation
    return std::make_unique<Trajectory>(interpolated_sim);
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

    // perform simulation using the controls, gdop config and a high number of nodes
    int  high_node_count    = 10 * gdop.mesh.node_count;
    auto simulation_result  = (*simulation)(gdop, extracted_controls, high_node_count,
                                             0.0, gdop.mesh.tf, exctracted_x0.raw());

    // result of high resolution simulation is interpolated onto lower resolution mesh
    auto interpolated_sim   = simulation_result->interpolate_onto_mesh(gdop.mesh, gdop.collocation);

    // calculate errors for each state in given norm (between provided and simulated states)
    auto errors             = trajectory.state_errors(interpolated_sim, norm);

    bool is_valid = true;

    FixedTableFormat<4> ftf = {{7,             13,            11,            4},
                               {Align::Center, Align::Center, Align::Center, Align::Center}};

    LOG_START_MODULE(ftf, "Simulation-Based Verification");

    LOG_HEADER(ftf, "State", fmt::format("Error [{}]", Linalg::norm_to_string(norm)), "Tolerance", "Pass");
    LOG_DASHES(ftf);

    for (size_t x_idx = 0; x_idx < trajectory.x.size(); x_idx++) {
        f64 tol = tolerances[x_idx];
        f64 err = errors[x_idx];
        bool pass = (err <= tol);

        LOG_ROW(ftf,
            fmt::format("x[{}]", x_idx),
            fmt::format("{:.3e}", err),
            fmt::format("{:.3e}", tol),
            pass ? "PASS" : "FAIL");

        if (!pass) {
            is_valid = false;
        }
    }

    LOG_DASHES(ftf);

    if (is_valid) {
        LOG_SUCCESS("All state errors within tolerances.");
    } else {
        LOG_WARNING("One or more state errors exceeded tolerances.");
    }

    LOG_DASHES_LN(ftf);

    return is_valid;
}

// interpolate trajectory to new mesh with collocation scheme - polynomial interpolation
std::unique_ptr<Trajectory> PolynomialInterpolation::operator()(const GDOP& gdop,
                                                                 const Mesh& new_mesh,
                                                                 const Trajectory& trajectory) {
    const auto& old_mesh = gdop.mesh;
    const auto& colloc   = gdop.collocation;

    const auto& old_x = trajectory.x;
    const auto& old_u = trajectory.u;
    const auto& old_p = trajectory.p;

    const int x_size = static_cast<int>(old_x.size());
    const int u_size = static_cast<int>(old_u.size());
    const int grid_size = new_mesh.node_count + 1;

    std::vector<f64> new_t;
    std::vector<std::vector<f64>> new_x(x_size, std::vector<f64>(grid_size));
    std::vector<std::vector<f64>> new_u(u_size, std::vector<f64>(grid_size));

    new_t.push_back(0.0);

    for (int x_idx = 0; x_idx < x_size; x_idx++) {
        new_x[x_idx][0] = old_x[x_idx][0]; // x(t_0)
    }

    for (int u_idx = 0; u_idx < u_size; u_idx++) {
        new_u[u_idx][0] = old_u[u_idx][0]; // u(t_0)
    }

    int global_grid_index = 1;
    int current_old_interval = 0;
    int offset = 0;

    for (int i = 0; i < new_mesh.intervals; i++) {
        int scheme_new = new_mesh.nodes[i];

        for (int j = 0; j < scheme_new; j++) {
            f64 t_query = new_mesh.t[i][j];
            new_t.push_back(t_query);

            while (current_old_interval + 1 < old_mesh.intervals &&
                old_mesh.grid[current_old_interval + 1] < t_query) {
                offset += old_mesh.nodes[current_old_interval];
                current_old_interval++;
            }

            const int old_p_order = old_mesh.nodes[current_old_interval];
            const bool contains_zero = true;

            f64 t_start = old_mesh.grid[current_old_interval];
            f64 t_end   = t_start + old_mesh.delta_t[current_old_interval];

            for (int x_idx = 0; x_idx < x_size; x_idx++) {
                const f64* x_vals = &old_x[x_idx][offset];
                new_x[x_idx][global_grid_index] = colloc.interpolate(
                    old_p_order, contains_zero, x_vals, 1,
                    t_start, t_end, t_query
                );
            }

            for (int u_idx = 0; u_idx < u_size; u_idx++) {
                const f64* u_vals = &old_u[u_idx][offset];
                new_u[u_idx][global_grid_index] = colloc.interpolate(
                    old_p_order, contains_zero, u_vals, 1,
                    t_start, t_end, t_query
                );
            }

            global_grid_index++;
        }
    }

    assert(global_grid_index == grid_size);
    return std::make_unique<Trajectory>(
        std::move(new_t), std::move(new_x), std::move(new_u),
        old_p, InterpolationMethod::LINEAR);
}


// L2BoundaryNorm reinit for next call
void L2BoundaryNorm::reinit(const GDOP& gdop) {
    phase_one_iteration = 0;
    phase_two_iteration = 0;
    max_phase_one_iterations = 0;
    max_phase_two_iterations = 15;

    // on-interval
    lambda         = 0.0;
    mesh_size_zero = gdop.mesh.intervals;

    // corner
    CTOL_1 = FixedVector<f64>(gdop.off_u);
    CTOL_2 = FixedVector<f64>(gdop.off_u);
    for (size_t i = 0; i < CTOL_1.size(); i++) { CTOL_1[i] = 0.1; }
    for (size_t i = 0; i < CTOL_2.size(); i++) { CTOL_2[i] = 0.1; }
}

// L2BoundaryNorm mesh refinement algorithm
std::unique_ptr<MeshUpdate> L2BoundaryNorm::operator()(const GDOP& gdop) {
    const auto& mesh = gdop.mesh;
    const auto& collocation = gdop.collocation;

    // constant degree for all intervals
    const int p = mesh.nodes[0];

    bool terminated = true;
    FixedVector<f64> new_grid;

    // L2BoundaryNorm
    if (phase_one_iteration < max_phase_one_iterations) {
        // phase I - full bisection case
        terminated = false;
        new_grid = FixedVector<f64>(2 * mesh.grid.size() - 1);
        for (int i = 0; i < mesh.grid.int_size() - 1; i++) {
            new_grid[2 * i] = mesh.grid[i];
            new_grid[2 * i + 1] = 0.5 * (mesh.grid[i] + mesh.grid[i + 1]);
        }
        new_grid[2 * mesh.grid.size() - 2] = mesh.tf;
        phase_one_iteration++;
    }
    else if (phase_two_iteration < max_phase_two_iterations) {
        // phase II - non-smoothness detection
        std::set<f64> set_new_grid;

        // p', p'' in `Enhancing Collocation-Based Dynamic Optimization through Adaptive Mesh Refinement` (Algorithm 3.2)
        FixedVector<f64> p_1(collocation.get_max_scheme() + 1);
        FixedVector<f64> p_2(collocation.get_max_scheme() + 1);

        // p'_i(t_{i+1}) and p''_i(t_{i+1}) for boundary condition
        bool has_last_boundary = false;
        f64 p_boundary_1_last_end;
        f64 p_boundary_2_last_end;
        f64 p_boundary_1_this_end;
        f64 p_boundary_2_this_end;

        for (auto const& u_vec : gdop.optimal_solution->u) {
            // compute range of u
            auto [min_it, max_it] = std::minmax_element(u_vec.begin(), u_vec.end());
            f64 u_range = *max_it - *min_it;

            // TODO: what about nominals? What about u == const. ?
            f64 TOL_1 = u_range * pow(10, -lambda) / mesh_size_zero;
            f64 TOL_2 = TOL_1 / 2;

            for (int i = 0; i < mesh.grid.int_size() - 1; i++) {
                int start_index = mesh.acc_nodes[i][0] + (int)(i != 0); // u(t_{i}), offset stems from 0 not contained in the mesh

                // apply D to the controls on interval i
                collocation.diff_matrix_multiply(mesh.nodes[i], &u_vec[start_index], p_1.raw());   // p'  = D^(1) * u_hat
                collocation.diff_matrix_multiply(mesh.nodes[i], p_1.raw(), p_2.raw());             // p'' = D^(2) * u_hat = D^(1) * p'

                // remember last value for next boundary condition
                p_boundary_1_this_end = p_1.back();
                p_boundary_2_this_end = p_2.back();

                // retrieve pointers at index 1 (ignore index 0 as its not needed for the condition) - note index 0 untouched!
                f64 *q_1 = &p_1[1];
                f64 *q_2 = &p_2[1];

                // element-wise square (inplace)
                Linalg::square(mesh.nodes[i], q_1); // q'  = ((p_1')^2,  (p_2')^2,  ..., (p_m')^2)
                Linalg::square(mesh.nodes[i], q_2); // q'' = ((p_1'')^2, (p_2'')^2, ..., (p_m'')^2)

                // calculate the L_2 norm with exact quadrature
                f64 norm_p_1 = sqrt(collocation.integrate(mesh.nodes[i], q_1)); // norm_p_1 := sqrt(b^T * q')
                f64 norm_p_2 = sqrt(collocation.integrate(mesh.nodes[i], q_2)); // norm_p_2 := sqrt(b^T * q'')

                // === on-interval condition ===
                if (norm_p_1 > TOL_1 || norm_p_2 > TOL_2) {
                    terminated = false;

                    // left interval
                    if (i > 0) {
                        set_new_grid.insert(0.5 * (mesh.grid[i - 1] + mesh.grid[i]));
                    }

                    // this / center interval
                    set_new_grid.insert(0.5 * (mesh.grid[i] + mesh.grid[i + 1]));

                    // right interval
                    if (i < mesh.grid.int_size() - 2) {
                        set_new_grid.insert(0.5 * (mesh.grid[i + 1] + mesh.grid[i + 2]));
                    }
                }

                // === boundary condition: Compare end of i-1 with start of i ===
                if (has_last_boundary) {
                    f64 p1_i1 = p_1[0];  // this interval, start i - index was untouched in on-interval computation
                    f64 p2_i1 = p_2[0];

                    // boundary condition plus-1-error
                    f64 ERR_1 = std::abs(p1_i1 - p_boundary_1_last_end) / (1.0 + std::min(std::abs(p1_i1), std::abs(p_boundary_1_last_end)));
                    f64 ERR_2 = std::abs(p2_i1 - p_boundary_2_last_end) / (1.0 + std::min(std::abs(p2_i1), std::abs(p_boundary_2_last_end)));

                    if (ERR_1 > TOL_1 || ERR_2 > TOL_2) {
                        terminated = false;

                        // insert this / center midpoint + left / previous midpoint
                        set_new_grid.insert(0.5 * (mesh.grid[i] + mesh.grid[i + 1]));
                        set_new_grid.insert(0.5 * (mesh.grid[i - 1] + mesh.grid[i]));
                    }
                }

                // Save last p_k for next round (last = end of this interval)
                p_boundary_1_last_end = p_boundary_1_this_end;
                p_boundary_2_last_end = p_boundary_2_this_end;
                has_last_boundary = true;

                set_new_grid.insert(mesh.grid[i]);
            }
        }

        set_new_grid.insert(mesh.grid.back());

        // create new_grid if not terminated
        if (!terminated) {
            new_grid = FixedVector<f64>(set_new_grid.begin(), set_new_grid.end());
        }

        phase_two_iteration++;
    }

    if (terminated) {
        return nullptr;
    }

    // set constant polynomial degree
    FixedVector<int> new_nodes(new_grid.size() - 1);
    for (int i = 0; i < new_nodes.int_size(); i++) {
        new_nodes[i] = p;
    }

    return std::make_unique<MeshUpdate>(std::move(new_grid),
                                        std::move(new_nodes));
}

// default strategy collection
Strategies Strategies::default_strategies() {
    Strategies strategy;
    strategy.initialization  = std::make_shared<DefaultConstantInitialization>();
    strategy.simulation      = std::make_shared<DefaultNoSimulation>();
    strategy.simulation_step = std::make_shared<DefaultNoSimulationStep>();
    strategy.mesh_refinement = std::make_shared<DefaultNoMeshRefinement>();
    strategy.interpolation   = std::make_shared<DefaultLinearInterpolation>();
    strategy.emitter         = std::make_shared<DefaultNoEmitter>();
    strategy.verifier        = std::make_shared<DefaultNoVerifier>();
    strategy.scaling_factory = std::make_shared<DefaultNoScalingFactory>();
    return strategy;
};

} // namespace GDOP
