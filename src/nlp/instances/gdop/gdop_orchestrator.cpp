#include "gdop_orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    // reset strategies
    strategies->reset(gdop);

    // create initial guess
    auto initial_guess = strategies->get_initial_guess(gdop);

    // set scaling
    gdop.set_scaling_factory(strategies->scaling_factory);

    // TODO: fix simulation-based verify / update, probably interpolation error of controls

    for(;;) {
        // set initial guess initially or after refinement
        gdop.set_initial_guess(std::move(initial_guess));

        // optimizer loop
        solver.optimize();

        // === mesh refinement ===

        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop.get_mesh(), *gdop.get_optimal_solution());

        if (!mesh_update) { break; }

        // 2. create refined Mesh
        auto refined_mesh = Mesh(std::move(mesh_update));

        // 3. interpolate (x*, lambda*, z*) to new mesh -> new initial guess
        initial_guess = strategies->get_refined_initial_guess(gdop.get_mesh(), refined_mesh, *gdop.get_optimal_solution());
        solver.solver_settings.set(NLP::Option::WarmStart, true);

        initial_guess->costates->to_csv("costates_interp.csv");
        initial_guess->lower_costates->to_csv("lower_costates_interp.csv");
        initial_guess->upper_costates->to_csv("upper_costates_interp.csv");

        // 4. update gdop with new mesh
        gdop.update(std::move(refined_mesh));
    }

    // verify by simulation
    strategies->verify(gdop, *gdop.get_optimal_solution());

    // emit optimal solution, maybe set verify only here
    strategies->emit(*gdop.get_optimal_solution()->primals);

    gdop.get_optimal_solution()->costates->to_csv("costates_final.csv");
    gdop.get_optimal_solution()->lower_costates->to_csv("lower_costates_final.csv");
    gdop.get_optimal_solution()->upper_costates->to_csv("upper_costates_final.csv");
}

} // namespace GDOP
