#ifndef OPT_OM_STRATEGIES_H
#define OPT_OM_STRATEGIES_H

#include "simulation/simulation_runtime.h"

#include <src/nlp/instances/gdop/gdop.h>

#include "info_gdop.h"
#include "evaluations.h"

namespace OpenModelica {


struct AuxiliaryTrajectory {
    Trajectory& trajectory;
    InfoGDOP& info;
    SOLVER_INFO* solver_info;
};

struct AuxiliaryControls {
    const ControlTrajectory& controls;
    InfoGDOP& info;
    FixedVector<f64> u_interpolation;
};

void initialize_model(InfoGDOP& info);
void emit_to_result_file(Trajectory& trajectory, InfoGDOP& info);

struct ConstantInitialization : public GDOP::Initialization {
    InfoGDOP& info;
    ConstantInitialization(InfoGDOP& info);

    std::unique_ptr<Trajectory> operator()(GDOP::GDOP& gdop) override;
};

// TODO: maybe think about how we provide these configurations to the Strategy?
// Do we want to pass a grant SimulationInfo { num_steps, start_time, stop_time } ?
struct Simulation : public GDOP::Simulation {
    InfoGDOP& info;
    SOLVER_METHOD solver;

    Simulation(InfoGDOP& info, SOLVER_METHOD solver);

    std::unique_ptr<Trajectory> operator()(GDOP::GDOP& gdop, const ControlTrajectory& u, int num_steps,
                                           f64 start_time, f64 stop_time, f64* x_start_values) override;
};

struct SimulationStep : public GDOP::SimulationStep {
    std::shared_ptr<Simulation> simulation;

    SimulationStep(std::shared_ptr<Simulation> simulation);

    std::unique_ptr<Trajectory> operator()(GDOP::GDOP& gdop, const ControlTrajectory& u,
                                           f64 start_time, f64 stop_time, f64* x_start_values) override;
};

// Strategies object
GDOP::Strategies default_strategies(InfoGDOP& info, SOLVER_METHOD solver);

} // namespace OpenModelica

#endif // OPT_OM_STRATEGIES_H
