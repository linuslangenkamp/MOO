#ifndef OPT_GDOP_STRATEGIES_H
#define OPT_GDOP_STRATEGIES_H

#include <functional>
#include <memory>
#include <src/base/trajectory.h>
#include <src/nlp/nlp.h>
#include <src/base/log.h>

// Strategies define interchangeable behaviors for key stages such as initialization, simulation,
// mesh refinement, result emission, and optimality verification in the GDOP optimization process.
// This file offers simple default and generic advanced strategy implementations that may be used.

// Each strategy can optionally modify the GDOP state when necessary (e.g. interpolation in mesh refinements),
// otherwise, GDOP is passed as const & to ensure safe, read-only access by default.
// This modular approach promotes flexibility, allowing strategies to be swapped or combined easily,
// simplifies testing and maintenance, and supports extensibility for custom behaviors.

// -- Base Strategy interfaces --

namespace GDOP {

class GDOP;


/**
 * @brief Strategy for initializing the GDOP.
 *
 * This strategy is responsible for creating an initial guess for the variables.
 *
 * Implementations may use bounds, analytical guesses, or results of simulation.
 *
 * @param gdop Optimization problem (read-only).
 * @return A unique_ptr to a PrimalDualTrajectory object representing the initial state guess.
 */
class Initialization {
public:
    virtual std::unique_ptr<PrimalDualTrajectory> operator()(const GDOP& gdop) = 0; 
};

/**
 * @brief Strategy for re-initializing the GDOP.
 *
 * This strategy is responsible for creating an initial guess for the variables,
 * based on an already existing trajectory, that is optimal on old_mesh.
 *
 * Implementations may:
 * - interpolate the provided solution in some clever way onto new_mesh.
 * - leave lower, upper and standard costates as nullptr, so only the old x is used.
 *
 * @param old_mesh Old Mesh of the GDOP
 * @param new_mesh New Mesh of the GDOP (interpolate onto this)
 * @param collocation Given Collocation scheme
 * @param trajectory Optimal trajectory on old_mesh
 * @return A unique_ptr to a PrimalDualTrajectory object representing the new, interpolated guess.
 */
class RefinedInitialization {
public:
    virtual std::unique_ptr<PrimalDualTrajectory> operator()(const Mesh& old_mesh,
                                                             const Mesh& new_mesh,
                                                             const Collocation& collocation,
                                                             const PrimalDualTrajectory& trajectory) = 0; 
};

/**
 * @brief Strategy for simulating the full system over the entire time horizon.
 *
 * @param controls Control inputs to apply.
 * @param num_steps Number of integration steps.
 * @param start_time Starting time of the simulation.
 * @param stop_time Final time of the simulation.
 * @param x_start_values Initial state values at start_time.
 * @return A unique_ptr to a Trajectory representing result file of the simulation.
 */
class Simulation {
public:
    virtual std::unique_ptr<Trajectory> operator()(const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) = 0;
};

/**
 * @brief Strategy for simulating a short segment, i.e. one step.
 *
 * Useful for finer-grained validation or model checking.
 *
 * @param controls Control input to apply (must include the time interval).
 * @param start_time Start of the step.
 * @param stop_time End of the step.
 * @param x_start_values Initial state at start_time.
 * @return A unique_ptr to the resulting trajectory segment.
 */
class SimulationStep {
public:
    virtual std::unique_ptr<Trajectory> operator()(const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) = 0;
};

/**
 * @brief Strategy for refining the mesh points and polynomial degrees.
 */
class MeshRefinement {
public:
    virtual void reset(const GDOP& gdop) = 0;
    virtual std::unique_ptr<MeshUpdate> operator()(const Mesh& mesh,
                                                   const Collocation& collocation,
                                                   const PrimalDualTrajectory& trajectory) = 0;
};

/**
 * @brief Strategy for interpolating a Trajectory onto a new Mesh
 *
 * @param gdop Optimization problem with old Mesh.
 * @param new_mesh New Mesh to interpolate onto.
 * @param trajectory Trajectory to interpolate.
 * @return A unique_ptr to the resulting interpolated trajectory.
 */
class Interpolation {
public:
    virtual std::vector<f64> operator()(const Mesh& old_mesh,
                                        const Mesh& new_mesh,
                                        const Collocation& collocation,
                                        const std::vector<f64>& values,
                                        bool contains_zero) = 0;
};

/**
 * @brief Strategy for emitting output, such as writing CSV, MAT files or logging.
 *
 * Called after optimization finishes to output the resulting trajectory.
 *
 * @param gdop Optimization problem.
 * @param trajectory Final trajectory.
 * @return 0 on success, nonzero on failure.
 */
class Emitter {
public:
    virtual int operator()(const Trajectory& trajectory) = 0;
};


// TODO: this should accept GDOP&, Trajectory& (optimal solution), CostateTrajectory& (dual optimum)
/**
 * @brief Strategy for verifying optimality / quality post-optimization.
 *
 * @param gdop Optimization problem.
 * @param trajectory Final trajectory.
 * @return true if verified successfully, false otherwise.
 */
class Verifier {
public:
    virtual bool operator()(const GDOP& gdop, const PrimalDualTrajectory& trajectory) = 0;
};

/**
 * @brief Strategy for creating and injecting NLP variable/function scaling.
 *
 * This factory creates a `NLP::Scaling` object based on the GDOP model.
 *
 * Once returned, the created scaling object will be **set into the `NLP` instance** 
 * (which is a parent of `GDOP`) and automatically used in solver routines.
 *
 * @param gdop Optimization problem.
 * @return A shared pointer to the created `NLP::Scaling` object.
 */
class ScalingFactory {
public:
    virtual std::shared_ptr<NLP::Scaling> operator()(const GDOP& gdop) = 0;
};

// ==================== Default Strategy implementations ====================

class DefaultNoSimulation : public Simulation {
public:
    std::unique_ptr<Trajectory> operator()(const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) override;
};

class DefaultNoSimulationStep : public SimulationStep {
public:
    std::unique_ptr<Trajectory> operator()(const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) override;
};

class DefaultNoMeshRefinement : public MeshRefinement {
public:
    void reset(const GDOP& gdop) override;
    std::unique_ptr<MeshUpdate> operator()(const Mesh& mesh, const Collocation& collocation, const PrimalDualTrajectory& trajectory) override;
};

class DefaultLinearInterpolation : public Interpolation {
public:
    std::vector<f64> operator()(const Mesh& old_mesh,
                                const Mesh& new_mesh,
                                const Collocation& collocation,
                                const std::vector<f64>& values,
                                bool contains_zero) override;
};

class DefaultInterpolationRefinedInitialization : public RefinedInitialization {
public:
    std::shared_ptr<Interpolation> interpolation;
    bool interpolate_primals;
    bool interpolate_costates_constraints;
    bool interpolate_costates_bounds;

    DefaultInterpolationRefinedInitialization(std::shared_ptr<Interpolation> interpolation_,
                                              bool interpolate_primals_,
                                              bool interpolate_costates_constraints_,
                                              bool interpolate_costates_bounds_);

    std::unique_ptr<PrimalDualTrajectory> operator()(const Mesh& old_mesh,
                                                     const Mesh& new_mesh,
                                                     const Collocation& collocation,
                                                     const PrimalDualTrajectory& trajectory) override;
};

class DefaultNoEmitter : public Emitter {
public:
    int operator()(const Trajectory& trajectory) override;
};

class DefaultNoVerifier : public Verifier {
public:
    bool operator()(const GDOP& gdop, const PrimalDualTrajectory& trajectory) override;
};

// -- simple default scaling (no scaling) --
class DefaultNoScalingFactory : public ScalingFactory {
public:
    std::shared_ptr<NLP::Scaling> operator()(const GDOP& gdop) override;
};

// -- simple default initialization (checks bounds and chooses initial value depending on that) --
class DefaultConstantInitialization : public Initialization {
public:
    std::unique_ptr<PrimalDualTrajectory> operator()(const GDOP& gdop) override;
};

// ==================== more advanced Strategies ====================

// -- uses Collocation scheme to interpolate States and Controls --
class PolynomialInterpolation : public Interpolation {
public:
    std::vector<f64> operator()(const Mesh& old_mesh,
                                const Mesh& new_mesh,
                                const Collocation& collocation,
                                const std::vector<f64>& values,
                                bool contains_zero) override;
};

// -- combined Strategy (simple Initialization, extract Controls, simulate) --
class SimulationInitialization : public Initialization {
public:
    std::shared_ptr<Initialization> initialization;
    std::shared_ptr<Simulation>     simulation;

    SimulationInitialization(std::shared_ptr<Initialization> initialization, std::shared_ptr<Simulation> simulation);

    std::unique_ptr<PrimalDualTrajectory> operator()(const GDOP& gdop) override;
};

// -- L2-Boundary-Norm Mesh Refinement Strategy --
class L2BoundaryNorm : public MeshRefinement {
public:
    int phase_one_iteration;
    int phase_two_iteration;
    int max_phase_one_iterations;
    int max_phase_two_iterations;

    // on-interval
    f64 mesh_lambda;
    f64 mesh_size_zero;

    // corner
    FixedVector<f64> CTOL_1;
    FixedVector<f64> CTOL_2;

    void reset(const GDOP& gdop) override;

    std::unique_ptr<MeshUpdate> operator()(const Mesh& mesh, const Collocation& collocation, const PrimalDualTrajectory& trajectory) override;
};

// -- emit optimal solution to csv --
class CSVEmitter : public Emitter {
public:
    std::string filename;

    CSVEmitter(std::string filename);

    int operator()(const Trajectory& trajectory) override;
};

// -- verify optimality by full simulation and state comparison with given norm --
class SimulationVerifier : public Verifier {
public:
    std::shared_ptr<Simulation> simulation;
    Linalg::Norm norm;
    FixedVector<f64> tolerances;

    SimulationVerifier(std::shared_ptr<Simulation> simulation, Linalg::Norm norm, FixedVector<f64>&& tolerances);

    bool operator()(const GDOP& gdop, const PrimalDualTrajectory& trajectory) override;
};

// TODO: add costates verifier

// ==================== Strategies Object ====================

/**
 * @brief Aggregates all strategy components into a single object.
 *
 * This class holds shared pointers to each pluggable strategy interface:
 * initialization, simulation, mesh refinement, emission, verification, scaling ...
 *
 * You can provide your own strategy objects or use the defaults via `default_strategies()`.
 */
class Strategies {
public:
    std::shared_ptr<Initialization>        initialization;
    std::shared_ptr<RefinedInitialization> refined_initialization;
    std::shared_ptr<Simulation>            simulation;
    std::shared_ptr<SimulationStep>        simulation_step;
    std::shared_ptr<MeshRefinement>        mesh_refinement;
    std::shared_ptr<Interpolation>         interpolation;
    std::shared_ptr<Emitter>               emitter;
    std::shared_ptr<Verifier>              verifier;
    std::shared_ptr<ScalingFactory>        scaling_factory;

    static Strategies default_strategies();

    auto get_initial_guess(const GDOP& gdop) {
        return (*initialization)(gdop);
    }

    auto get_refined_initial_guess(const Mesh& old_mesh, const Mesh& new_mesh, const Collocation& collocation, const PrimalDualTrajectory& trajectory) {
        return (*refined_initialization)(old_mesh, new_mesh, collocation, trajectory);
    }

    auto simulate(const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
        return (*simulation)(controls, num_steps, start_time, stop_time, x_start_values);
    }

    auto simulate_step(const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) {
        return (*simulation_step)(controls, start_time, stop_time, x_start_values);
    }

    auto detect(const Mesh& mesh, const Collocation& collocation, const PrimalDualTrajectory& trajectory) {
        return (*mesh_refinement)(mesh, collocation, trajectory);
    }

    auto interpolate(const Mesh& old_mesh, const Mesh& new_mesh, const Collocation& collocation, const std::vector<f64>& values, bool contains_zero) {
        return (*interpolation)(old_mesh, new_mesh, collocation, values, contains_zero);
    }

    auto emit(const Trajectory& trajectory) {
        return (*emitter)(trajectory);
    }

    auto verify(const GDOP& gdop, const PrimalDualTrajectory& trajectory) {
        return (*verifier)(gdop, trajectory);
    }

    auto create_scaling(const GDOP& gdop) {
        return (*scaling_factory)(gdop);
    }

    void reset(const GDOP& gdop) {
        // add others if we have an internal state that changes during optimization (e.g. mesh refinement iteration count)
        mesh_refinement->reset(gdop);
    }
};

} // namespace GDOP

#endif // OPT_GDOP_STRATEGIES_H
