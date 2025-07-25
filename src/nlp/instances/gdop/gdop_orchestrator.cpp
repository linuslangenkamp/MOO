#include "gdop_orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    // initialize GDOP (creates sparsity, bounds, ...)
    gdop.init();

    // reset strategies
    strategies->reset(gdop);

    // create initial guess
    auto initial_guess = strategies->get_initial_guess(gdop);

    // TODO: fix simulation-based verify / update, probably interpolation error of controls

    for(;;) {
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
        //gdop.optimal_solution->costates->to_csv("costates.csv");

        // 3. interpolate (x*, lambda*, z*) to new mesh -> new initial guess
        initial_guess = strategies->get_refined_initial_guess(gdop.mesh, refined_mesh, gdop.collocation, *gdop.optimal_solution);
        solver.solver_settings.set(NLP::Option::WarmStart, true);

        initial_guess->costates->to_csv("costates_interp.csv");
        initial_guess->lower_costates->to_csv("lower_costates_interp.csv");
        initial_guess->upper_costates->to_csv("upper_costates_interp.csv");

        // 4. update gdop with new mesh
        gdop.update(std::move(refined_mesh));
    }

    // verify by simulation
    strategies->verify(gdop, *gdop.optimal_solution);

    // emit optimal solution, maybe set verify only here
    strategies->emit(*gdop.optimal_solution->primals);

    gdop.optimal_solution->costates->to_csv("costates_final.csv");
    gdop.optimal_solution->lower_costates->to_csv("lower_costates_final.csv");
    gdop.optimal_solution->upper_costates->to_csv("upper_costates_final.csv");

}

} // namespace GDOP
