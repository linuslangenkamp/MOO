#ifndef OPT_GDOP_STRATEGIES_H
#define OPT_GDOP_STRATEGIES_H

#include <functional>
#include <memory>
#include <src/base/mesh.h>
#include <src/nlp/instances/gdop/gdop.h>

namespace StrategiesGDOP {

// -- Base Strategy interfaces --

struct Initialization {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop) = 0; 
};

struct Simulation {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop, const ControlTrajectory& u) = 0;
};

struct MeshRefinement {
    virtual std::unique_ptr<Trajectory> operator()(GDOP& gdop) = 0;
};

// -- Default Strategy implementations --

struct DefaultConstantInitialization : public Initialization {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop) override;
};

struct DefaultNoSimulation : public Simulation {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop, const ControlTrajectory& u) override;
};

struct DefaultNoMeshRefinement : public MeshRefinement {
    std::unique_ptr<Trajectory> operator()(GDOP& gdop) override;
};

// -- Strategy container --

struct Container {
    std::unique_ptr<Initialization> initialization  = std::make_unique<DefaultConstantInitialization>();
    std::unique_ptr<Simulation>     simulation      = std::make_unique<DefaultNoSimulation>();
    std::unique_ptr<MeshRefinement> mesh_refinement = std::make_unique<DefaultNoMeshRefinement>();

    auto initialize(GDOP& gdop) {
        return (*initialization)(gdop);
    }

    auto simulate(GDOP& gdop, const ControlTrajectory& u) {
        return (*simulation)(gdop, u);
    }

    auto refine(GDOP& gdop) {
        return (*mesh_refinement)(gdop);
    }
};

} // namespace StrategiesGDOP

#endif // OPT_GDOP_STRATEGIES_H
