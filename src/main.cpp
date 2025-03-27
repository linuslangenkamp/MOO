#include <iostream>

#include "base/collocation.h"
#include "impl/gdop/gdop.h"
#include "impl/gdop/problem.h"
#include "interfaces/ipopt_solver.h"

#include "impl/gdop/test_problem_impl.h"

int main() {
    Mesh mesh = Mesh::createEquidistantMeshFixedDegree(1, 1, 100);
    std::unique_ptr<Collocation> radau = std::make_unique<Collocation>();

    // M(x) = x^2
    FixedVector<FunctionMR> mr(1);
    mr[0].jac.dxf.push_back(JacobianSparsity({0, nullptr}));
    mr[0].hes.dxf_dxf.push_back(HessianSparsity({0, 0, nullptr}));

    // f(x, u) = cos(x) + u
    FixedVector<FunctionLFG> lfg(1);
    lfg[0].jac.dx.push_back(JacobianSparsity({0, nullptr}));
    lfg[0].jac.du.push_back(JacobianSparsity({0, nullptr}));
    lfg[0].hes.dx_dx.push_back(HessianSparsity({0, 0, nullptr}));

    FixedVector<Bounds> g_bounds(0);
    FixedVector<Bounds> r_bounds(0);

    FullSweepTestImpl fs = FullSweepTestImpl(std::move(lfg), std::make_shared<Mesh>(mesh), g_bounds);
    BoundarySweepTestImpl bs = BoundarySweepTestImpl(std::move(mr), std::make_shared<Mesh>(mesh), r_bounds);

    FixedVector<Bounds> x_bounds(1);
    FixedVector<Bounds> u_bounds(1);
    FixedVector<Bounds> p_bounds(0);
    FixedVector<std::optional<double>> x0_fixed(1);
    x0_fixed[0] = 0.5;
    FixedVector<std::optional<double>> xf_fixed(1);
    x_bounds[0].lb = -1;
    x_bounds[0].ub = 1;
    u_bounds[0].lb = -1;
    u_bounds[0].ub = 1;

    Problem problem(std::make_unique<FullSweepTestImpl>(std::move(fs)), std::make_unique<BoundarySweepTestImpl>(std::move(bs)), 
        x_bounds, u_bounds, p_bounds, x0_fixed, xf_fixed);
    std::cout << sizeof(FunctionLFG) << std::endl;
    // 0 guess
    Trajectory initial_guess = {{0, 1}, {{0.5, 0.5}}, {{0, 1}}, {}, InterpolationMethod::LINEAR};

    GDOP gdop(std::make_shared<Problem>(std::move(problem)), std::move(radau), mesh, initial_guess);

    IpoptSolver ipopt_solver(std::make_shared<GDOP>(std::move(gdop)), NULL);
    ipopt_solver.optimize();

    return 0;
}