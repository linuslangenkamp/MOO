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

class ConstantInitialization : public GDOP::Initialization {
public:
    InfoGDOP& info;
    ConstantInitialization(InfoGDOP& info);

    std::unique_ptr<PrimalDualTrajectory> operator()(const GDOP::GDOP& gdop) override;
};

// TODO: maybe think about how we provide these configurations to the Strategy?
// Do we want to pass a grant SimulationInfo { num_steps, start_time, stop_time } ?
class Simulation : public GDOP::Simulation {
public:
    InfoGDOP& info;
    SOLVER_METHOD solver;

    Simulation(InfoGDOP& info, SOLVER_METHOD solver);

    std::unique_ptr<Trajectory> operator()(const ControlTrajectory& controls, int num_steps,
                                           f64 start_time, f64 stop_time, f64* x_start_values) override;
};

class SimulationStep : public GDOP::SimulationStep {
public:
    std::shared_ptr<Simulation> simulation;

    SimulationStep(std::shared_ptr<Simulation> simulation);

    std::unique_ptr<Trajectory> operator()(const ControlTrajectory& controls,
                                           f64 start_time, f64 stop_time, f64* x_start_values) override;
};

class MatEmitter : public GDOP::Emitter {
public:
    InfoGDOP& info;

    MatEmitter(InfoGDOP& info);

    int operator()(const Trajectory& trajectory) override;
};

class NominalScalingFactory : public GDOP::ScalingFactory {
public:
    InfoGDOP& info;

    NominalScalingFactory(InfoGDOP& info) : info{info} {}

    std::shared_ptr<NLP::Scaling> operator()(const GDOP::GDOP& gdop) override;
};

GDOP::Strategies default_strategies(InfoGDOP& info);

} // namespace OpenModelica

#endif // OPT_OM_STRATEGIES_H
