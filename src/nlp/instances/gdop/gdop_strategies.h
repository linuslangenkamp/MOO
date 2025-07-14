#ifndef OPT_GDOP_STRATEGIES_H
#define OPT_GDOP_STRATEGIES_H

#include <functional>
#include <memory>
#include <src/base/mesh.h>

// -- Base Strategy interfaces --

namespace GDOP {

class GDOP;

struct Initialization {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop) = 0; 
};

struct Simulation {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop, const ControlTrajectory& u, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) = 0;
};

struct SimulationStep {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop, const ControlTrajectory& u, f64 start_time, f64 stop_time, f64* x_start_values) = 0;
};

struct MeshRefinement {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop) = 0;
};

// TODO: Emit to file?!

// TODO: Verify?!

// -- Default Strategy implementations --

struct DefaultConstantInitialization : public Initialization {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop) override;
};

struct DefaultNoSimulation : public Simulation {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop, const ControlTrajectory& u, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) override;
};

struct DefaultNoSimulationStep : public SimulationStep {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop, const ControlTrajectory& u, f64 start_time, f64 stop_time, f64* x_start_values) override;
};

struct DefaultNoMeshRefinement : public MeshRefinement {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop) override;
};

// -- Combined Strategy (Simple Initialization, Extract Controls, Simulate)
struct SimulationInitialization : public Initialization {
    std::shared_ptr<Initialization> initialization;
    std::shared_ptr<Simulation>     simulation;

    SimulationInitialization(std::shared_ptr<Initialization> initialization, std::shared_ptr<Simulation> simulation);

    std::unique_ptr<Trajectory> operator()(GDOP& gdop) override;
};

// -- Strategies --

struct Strategies {
    std::shared_ptr<Initialization> initialization;
    std::shared_ptr<Simulation>     simulation;
    std::shared_ptr<SimulationStep> simulation_step;
    std::shared_ptr<MeshRefinement> mesh_refinement;

    static Strategies default_strategies();

    auto initialize(GDOP& gdop) {
        return (*initialization)(gdop);
    }

    auto simulate(GDOP& gdop, const ControlTrajectory& u, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values) {
        return (*simulation)(gdop, u, num_steps, start_time, stop_time, x_start_values);
    }

    auto simulate_step(GDOP& gdop, const ControlTrajectory& u, f64 start_time, f64 stop_time, f64* x_start_values) {
        return (*simulation_step)(gdop, u, start_time, stop_time, x_start_values);
    }

    auto refine(GDOP& gdop) {
        return (*mesh_refinement)(gdop);
    }
};

} // namespace GDOP

#endif // OPT_GDOP_STRATEGIES_H
