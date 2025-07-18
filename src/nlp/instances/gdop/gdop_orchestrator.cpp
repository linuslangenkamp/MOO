#include "gdop_orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    // reset strategies
    strategies->reinit();

    // initialize GDOP (creates sparsity, bounds, ...)
    gdop.init();

    // create inital guess
    auto inital_guess = strategies->initialize(gdop);
    gdop.set_initial_guess(std::move(inital_guess));
    gdop.init_starting_point();

optimizer:
    // set scaling
    auto scaling = strategies->create_scaling(gdop);
    gdop.set_scaling(scaling);

    // optimizer loop
    solver.optimize();

    // mesh refinement
    if (true) {
        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop);

        if (!mesh_update) { goto finalize; }

        // 2. interpolate with old mesh and MeshUpdate to new Mesh -> new initial guess | what about lambda interpolation?
        // 3. update mesh and update the buffers and stuff in GDOP
        gdop.mesh.update(std::move(mesh_update), gdop.collocation); // TODO: make thin wrapper in GDOP
        gdop.init();
        gdop.problem.resize_buffers();
        // 4. set new initial guess

        // 5. goto solver.optimize()
        goto optimizer;
    }

finalize:
    // verify and emit optimal solution
    strategies->verify(gdop, *gdop.optimal_solution);
    strategies->emit(gdop, *gdop.optimal_solution);
}

} // namespace GDOP
