#ifndef OPT_GDOP_STRATEGIES_H
#define OPT_GDOP_STRATEGIES_H

#include <functional>
#include <memory>
#include <src/base/trajectory.h>
#include <src/base/log.h>

#include <src/nlp/nlp.h>

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

struct Initialization {
    virtual std::unique_ptr<Trajectory> operator()(const GDOP& gdop) = 0; 
};

struct Simulation {
    virtual std::unique_ptr<Trajectory> operator()(const GDOP& gdop, const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) = 0;
};

struct SimulationStep {
    virtual std::unique_ptr<Trajectory> operator()(const GDOP& gdop, const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) = 0;
};

struct MeshRefinement {
    virtual std::unique_ptr<Trajectory> operator()(const GDOP& gdop) = 0;
};

struct Emitter {
    virtual int operator()(const GDOP& gdop, const Trajectory& trajectory) = 0;
};

struct Verifier {
    virtual bool operator()(const GDOP& gdop, const Trajectory& trajectory) = 0;
};

// -- Default Strategy implementations --

struct DefaultNoSimulation : public Simulation {
    std::unique_ptr<Trajectory> operator()(const GDOP& gdop, const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) override;
};

struct DefaultNoSimulationStep : public SimulationStep {
    std::unique_ptr<Trajectory> operator()(const GDOP& gdop, const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) override;
};

struct DefaultNoMeshRefinement : public MeshRefinement {
    std::unique_ptr<Trajectory> operator()(const GDOP& gdop) override;
};

struct DefaultNoEmitter : public Emitter {
    int operator()(const GDOP& gdop, const Trajectory& trajectory) override;
};

struct DefaultNoVerifier : public Verifier {
    bool operator()(const GDOP& gdop, const Trajectory& trajectory) override;
};


// -- actual implementations --

// -- simple default initialization (checks bounds and chooses initial value depending on that) --
struct DefaultConstantInitialization : public Initialization {
    std::unique_ptr<Trajectory> operator()(const GDOP& gdop) override;
};

// -- combined Strategy (simple Initialization, extract Controls, simulate) --
struct SimulationInitialization : public Initialization {
    std::shared_ptr<Initialization> initialization;
    std::shared_ptr<Simulation>     simulation;

    SimulationInitialization(std::shared_ptr<Initialization> initialization, std::shared_ptr<Simulation> simulation);

    std::unique_ptr<Trajectory> operator()(const GDOP& gdop) override;
};

// -- emit optimal solution to csv --
struct CSVEmitter : public Emitter {
    std::string filename;

    CSVEmitter(std::string filename);

    int operator()(const GDOP& gdop, const Trajectory& trajectory) override;
};

// -- verify optimality by full simulation and state comparison with given norm --
struct SimulationVerifier : public Verifier {
    std::shared_ptr<Simulation> simulation;
    Linalg::Norm norm;
    FixedVector<f64> tolerances;

    SimulationVerifier(std::shared_ptr<Simulation> simulation, Linalg::Norm norm, FixedVector<f64>&& tolerances);

    bool operator()(const GDOP& gdop, const Trajectory& trajectory) override;
};

struct ScalingFactory {
    virtual std::shared_ptr<NLP::Scaling> operator()(const GDOP& gdop) = 0;
};

// TODO: add default scaling here also

// -- Strategies --

struct Strategies {
    std::shared_ptr<Initialization> initialization;
    std::shared_ptr<Simulation>     simulation;
    std::shared_ptr<SimulationStep> simulation_step;
    std::shared_ptr<MeshRefinement> mesh_refinement;
    std::shared_ptr<Emitter>        emitter;
    std::shared_ptr<Verifier>       verifier;
    std::shared_ptr<ScalingFactory> scaling_factory;

    static Strategies default_strategies();

    auto initialize(const GDOP& gdop) {
        return (*initialization)(gdop);
    }

    auto simulate(const GDOP& gdop, const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
        return (*simulation)(gdop, controls, num_steps, start_time, stop_time, x_start_values);
    }

    auto simulate_step(const GDOP& gdop, const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values) {
        return (*simulation_step)(gdop, controls, start_time, stop_time, x_start_values);
    }

    auto refine(const GDOP& gdop) {
        return (*mesh_refinement)(gdop);
    }

    auto emit(const GDOP& gdop, const Trajectory& trajectory) {
        return (*emitter)(gdop, trajectory);
    }

    auto verify(const GDOP& gdop, const Trajectory& trajectory) {
        return (*verifier)(gdop, trajectory);
    }

    auto create_scaling(const GDOP& gdop) {
        return (*scaling_factory)(gdop);
    }
};

} // namespace GDOP

#endif // OPT_GDOP_STRATEGIES_H
