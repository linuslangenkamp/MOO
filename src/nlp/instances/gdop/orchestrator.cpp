#include "orchestrator.h"

namespace GDOP {

void MeshRefinementOrchestrator::optimize() {
    strategies->reset(gdop);

    auto initial_guess = strategies->get_initial_guess(gdop);

    gdop.set_scaling_factory(strategies->scaling_factory);

    for(;;) {
        gdop.set_initial_guess(std::move(initial_guess));

        solver.optimize();

        // === mesh refinement ===

        // 1. detect intervals and degrees (new vectors)
        auto mesh_update = strategies->detect(gdop.get_mesh(), *gdop.get_optimal_solution());

        if (!mesh_update) { break; }

        // 2. create refined Mesh
        auto refined_mesh = Mesh::create_from_mesh_update(std::move(mesh_update));

        // 3. interpolate (x*, lambda*, z*) to new mesh -> new initial guess
        initial_guess = strategies->get_refined_initial_guess(gdop.get_mesh(), *refined_mesh, *gdop.get_optimal_solution());
        solver.solver_settings.set(NLP::Option::WarmStart, true);

        initial_guess->costates->to_csv("costates_interp.csv");
        initial_guess->lower_costates->to_csv("lower_costates_interp.csv");
        initial_guess->upper_costates->to_csv("upper_costates_interp.csv");

        // 4. update gdop with new mesh
        gdop.update(refined_mesh);
    }

    strategies->verify(gdop, *gdop.get_optimal_solution());
    strategies->emit(*gdop.get_optimal_solution()->primals);

    gdop.get_optimal_solution()->costates->to_csv("costates_final.csv");
    gdop.get_optimal_solution()->lower_costates->to_csv("lower_costates_final.csv");
    gdop.get_optimal_solution()->upper_costates->to_csv("upper_costates_final.csv");
}

} // namespace GDOP
