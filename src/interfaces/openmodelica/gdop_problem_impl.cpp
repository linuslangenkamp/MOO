#include "gdop_problem_impl.h"

std::shared_ptr<Problem> create_gdop_om(DATA* data, std::shared_ptr<Mesh> mesh) {
    int number_states = data->modelData->nStates;
    int number_controls = data->modelData->nInputVars;
    int number_parameters = 0; // TODO: add this feature
    int number_path_constraints = data->modelData->nOptimizeConstraints;
    int number_boundary_constrs = data->modelData->nOptimizeFinalConstraints; // TODO: add *generic boundary* constraints later also at t=t0

    FixedVector<FunctionMR> mr(1);
    FixedVector<FunctionLFG> lfg(2);


    FixedVector<Bounds> g_bounds(0);
    FixedVector<Bounds> r_bounds(1);
    r_bounds[0].lb = 0;
    r_bounds[0].ub = 0;

    std::unique_ptr<FullSweep> fs(new FullSweepOM(std::move(lfg), mesh, g_bounds));
    std::unique_ptr<BoundarySweep> bs(new BoundarySweepOM(std::move(mr), mesh, r_bounds));


    FixedVector<Bounds> x_bounds(1);
    FixedVector<Bounds> u_bounds(0);
    FixedVector<Bounds> p_bounds(1);
    FixedVector<std::optional<double>> x0_fixed(1);
    FixedVector<std::optional<double>> xf_fixed(1);

    return std::make_shared<Problem>(std::move(fs), std::move(bs), 
                                     std::move(x_bounds), std::move(u_bounds), std::move(p_bounds),
                                     std::move(x0_fixed), std::move(xf_fixed));
}


