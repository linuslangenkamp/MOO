
#include <iostream>

#include "base/fLGR.h"
#include "impl/gdop/gdop.h"
#include "impl/gdop/problem.h"
#include "interfaces/ipopt_solver.h"

#include "impl/gdop/test_problem_impl.h"

int main() {
    auto radau = fLGR();
    auto mesh = std::make_shared<Mesh>(Mesh::create_equidistant_fixed_stages(1, 15000, 5, radau));

    std::unique_ptr<fLGR> p_radau = std::make_unique<fLGR>(radau);

    // M(x) = x2
    FixedVector<FunctionMR> mr(1);
    mr[0].jac.dxf.push_back(JacobianSparsity({1, nullptr}));

    // x1' = -(u + u² / 2) * x1
    FixedVector<FunctionLFG> lfg(2);
    lfg[0].jac.dx.push_back(JacobianSparsity({0, nullptr}));
    lfg[0].jac.du.push_back(JacobianSparsity({0, nullptr}));
    lfg[0].hes.du_dx.push_back(HessianSparsity({0, 0, nullptr}));
    lfg[0].hes.du_du.push_back(HessianSparsity({0, 0, nullptr}));

    // x2' = u * x1
    lfg[1].jac.dx.push_back(JacobianSparsity({0, nullptr}));
    lfg[1].jac.du.push_back(JacobianSparsity({0, nullptr}));
    lfg[1].hes.du_dx.push_back(HessianSparsity({0, 0, nullptr}));

    FixedVector<Bounds> g_bounds(0);
    FixedVector<Bounds> r_bounds(0);

    std::unique_ptr<FullSweep> fs(new FullSweepTestImpl(std::move(lfg), mesh, g_bounds));
    std::unique_ptr<BoundarySweep> bs(new BoundarySweepTestImpl(std::move(mr), mesh, r_bounds));

    FixedVector<Bounds> x_bounds(2);
    FixedVector<Bounds> u_bounds(1);
    FixedVector<Bounds> p_bounds(0);
    FixedVector<std::optional<f64>> x0_fixed(2);
    FixedVector<std::optional<f64>> xf_fixed(2);
    x0_fixed[0] = 1;
    x0_fixed[1] = 0;

    u_bounds[0].lb = 0;
    u_bounds[0].ub = 5;

    auto problem = std::make_shared<Problem>(std::move(fs), std::move(bs), std::move(x_bounds), std::move(u_bounds), std::move(p_bounds), std::move(x0_fixed), std::move(xf_fixed));

    // 0 guess
    std::shared_ptr<Trajectory> initial_guess(new Trajectory{{0, 1}, {{1, 1}, {0, 0}}, {{0.5, 0.5}}, {}, InterpolationMethod::LINEAR});

    GDOP gdop(problem, std::move(radau), mesh, initial_guess);
    
    IpoptSolver ipopt_solver(std::make_shared<GDOP>(std::move(gdop)), NULL);
    ipopt_solver.optimize();

    return 0;
}
