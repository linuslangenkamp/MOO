
#include <iostream>

#include "base/collocation.h"
#include "impl/gdop/gdop.h"
#include "impl/gdop/problem.h"
#include "interfaces/ipopt_solver.h"

#include "impl/gdop/test_problem_impl.h"

int main() {
    Collocation coll = Collocation();
    auto mesh = std::make_shared<Mesh>(Mesh::createEquidistantMeshFixedDegree(250, 100, 9, coll));

    std::unique_ptr<Collocation> radau = std::make_unique<Collocation>(coll);

    // r(x) = x0 - p
    FixedVector<FunctionMR> mr(1);
    mr[0].jac.dx0.push_back(JacobianSparsity({0, nullptr}));
    mr[0].jac.dp.push_back(JacobianSparsity({0, nullptr}));

    // integral (x - p)²
    FixedVector<FunctionLFG> lfg(2);
    lfg[0].jac.dx.push_back(JacobianSparsity({0, nullptr}));
    lfg[0].jac.dp.push_back(JacobianSparsity({0, nullptr}));
    lfg[0].hes.dx_dx.push_back(HessianSparsity({0, 0, nullptr}));
    lfg[0].hes.dp_dx.push_back(HessianSparsity({0, 0, nullptr}));
    lfg[0].hes.dp_dp.push_back(HessianSparsity({0, 0, nullptr}));
    
    // T' = -0.4 * T + 120.2 + 0.2 * sin(4*pi*time)
    lfg[1].jac.dx.push_back(JacobianSparsity({0, nullptr}));

    FixedVector<Bounds> g_bounds(0);
    FixedVector<Bounds> r_bounds(1);
    r_bounds[0].lb = 0;
    r_bounds[0].ub = 0;

    std::unique_ptr<FullSweep> fs(new FullSweepTestImpl(std::move(lfg), mesh, g_bounds));
    std::unique_ptr<BoundarySweep> bs(new BoundarySweepTestImpl(std::move(mr), mesh, r_bounds));

    FixedVector<Bounds> x_bounds(1);
    FixedVector<Bounds> u_bounds(0);
    FixedVector<Bounds> p_bounds(1);
    FixedVector<std::optional<double>> x0_fixed(1);
    FixedVector<std::optional<double>> xf_fixed(1);

    auto problem = std::make_shared<Problem>(std::move(fs), std::move(bs), std::move(x_bounds), std::move(u_bounds), std::move(p_bounds), std::move(x0_fixed), std::move(xf_fixed));

    // 0 guess
    std::shared_ptr<Trajectory> initial_guess(new Trajectory{{0, 100}, {{300, 301}}, {}, {300}, InterpolationMethod::LINEAR});

    GDOP gdop(problem, std::move(radau), mesh, initial_guess);
    
    IpoptSolver ipopt_solver(std::make_shared<GDOP>(std::move(gdop)), NULL);
    ipopt_solver.optimize();

    return 0;
}
