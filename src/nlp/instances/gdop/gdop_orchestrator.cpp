#include "gdop_orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    // initialize GDOP (creates sparsity, bounds, ...)
    gdop.init();

    // reset strategies
    strategies->reinit(gdop);

    // create initial guess
    auto initial_guess = strategies->initialize(gdop);

    // TODO: set solver flags after refinement
    // TODO: fix simulation-baed verify / update (ask Bernhard)
    // TODO: dont use interpolation - make own strategy for reinit (including lambda interpolation)

    for(;;) {
        // create inital guess
        gdop.set_initial_guess(std::move(initial_guess));
        gdop.init_starting_point();

        // set scaling
        auto scaling = strategies->create_scaling(gdop);
        gdop.set_scaling(scaling);

        // optimizer loop
        solver.optimize();

        // mesh refinement

        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop.mesh, gdop.collocation, *gdop.optimal_solution);

        if (!mesh_update) { break; }

        // 2. create refined Mesh
        auto refined_mesh = Mesh(std::move(mesh_update), gdop.collocation);

        // 3. interpolate to new mesh -> new initial guess | what about lambda interpolation?
        // reinit: initial_guess = strategies->interpolate(gdop.mesh, refined_mesh, gdop.collocation, *gdop.optimal_solution);
        initial_guess = strategies->initialize(gdop);

        // 4. reinit gdop with new mesh
        gdop.reinit(std::move(refined_mesh));
    }

    // verify by simulation
    strategies->verify(gdop, *gdop.optimal_solution);

    // emit optimal solution, maybe set verify only here
    strategies->emit(*gdop.optimal_solution->primals);

    gdop.optimal_solution->costates->to_csv("costates.csv");
    gdop.optimal_solution->lower_costates->to_csv("lower_costates.csv");
    gdop.optimal_solution->upper_costates->to_csv("upper_costates.csv");
}

} // namespace GDOP
