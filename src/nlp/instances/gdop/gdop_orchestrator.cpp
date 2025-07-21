#include "gdop_orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    // initialize GDOP (creates sparsity, bounds, ...)
    gdop.init();

    // reset strategies
    strategies->reinit(gdop);

    // create initial guess
    auto initial_guess = strategies->initialize(gdop);

    // TODOS: - set solver flags after refinement, fix poly interpolation, fix verify / update, dont use interpolation - make own strategy for reinit

    for(;;) {
        // create inital guess
        gdop.set_initial_guess(std::move(initial_guess));
        gdop.init_starting_point();

        // set scaling
        auto scaling = strategies->create_scaling(gdop);
        gdop.set_scaling(scaling);

        // optimizer loop
        solver.optimize();

        // optional
        strategies->verify(gdop, *gdop.optimal_solution);

        // mesh refinement

        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop);

        if (!mesh_update) { break; }

        // 2. create refined Mesh
        auto refined_mesh = Mesh(std::move(mesh_update), gdop.collocation);

        // 3. interpolate to new mesh -> new initial guess | what about lambda interpolation?
        initial_guess = strategies->interpolate(gdop, refined_mesh, *gdop.optimal_solution);

        // 4. update mesh
        gdop.mesh.move_from(std::move(refined_mesh));

        // 5. update the buffers and init the GDOP
        gdop.init();
        gdop.problem.resize_buffers();
    }

    // emit optimal solution, maybe set verify only here
    strategies->emit(gdop, *gdop.optimal_solution);
}

} // namespace GDOP
