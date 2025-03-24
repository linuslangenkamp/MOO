#include <iostream>

#include "base/collocation.h"
#include "impl/gdop/gdop.h"
#include "impl/gdop/problem.h"
#include "interfaces/ipopt_solver.h"

#include "impl/gdop/test_problem_impl.h"

int main() {
    Mesh mesh = Mesh::createEquidistantMeshFixedDegree(10, 1, 3);
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

    FullSweepTestImpl fs = FullSweepTestImpl(lfg, std::make_shared<Mesh>(mesh), g_bounds);
    BoundarySweepTestImpl bs = BoundarySweepTestImpl(mr, std::make_shared<Mesh>(mesh), r_bounds);

    FixedVector<Bounds> x_bounds(1);
    FixedVector<Bounds> u_bounds(1);
    FixedVector<Bounds> p_bounds(0);
    FixedVector<std::optional<double>> x0_fixed(1);
    FixedVector<std::optional<double>> xf_fixed(1);
    x_bounds[0].lb = -1;
    x_bounds[0].ub = 1;
    u_bounds[0].lb = -1;
    u_bounds[0].ub = 1;

    Problem problem(std::make_unique<FullSweepTestImpl>(std::move(fs)), std::make_unique<BoundarySweepTestImpl>(std::move(bs)), 
        x_bounds, u_bounds, p_bounds, x0_fixed, xf_fixed);

    // 0 guess
    Trajectory initial_guess = {{0, 1}, {{0, 0.5}}, {{0, 1}}, {}, InterpolationMethod::LINEAR};

    GDOP gdop(std::make_shared<Problem>(std::move(problem)), std::move(radau), mesh, initial_guess);

    IpoptSolver ipopt_solver(std::make_shared<GDOP>(std::move(gdop)), NULL);
    ipopt_solver.optimize();

    /*FixedField<double, 3> x(1, 2, 5);
    FixedField<double, 1> q(2);
    std::vector<double> vec = {7, 8, 9, 10, 11};
    FixedVector<double> z(vec.begin(), vec.end());
    q.assign(&vec[2], 2);
    std::cout << q[1] << std::endl;
    q.print();
    Collocation radau = Collocation();
    BlockSparsity::createLowerTriangular(5, BlockType::Offset);
    Problem problem = Problem();
    problem.x_size = 1;
    problem.x_bounds = {{2, 4}};
    problem.x0_fixed = {std::nullopt};
    problem.xf_fixed = {std::nullopt};
    problem.u_size = 2;
    problem.u_bounds = {{2, 4}, {-1, 2}};
    problem.p_size = 0;
    problem.full->f_size = 1;
    problem.full->g_size = 1;
    problem.full->fg_size = 2;
    problem.full->g_bounds = {{1, 2}};
    problem.boundary->r_size = 1;
    problem.boundary->r_bounds = {{125, 2555}};
    Trajectory guess = {{0, 1}, {{1, 5}}, {{0, 2}, {3, 1}}, {1, 2}, InterpolationMethod::LINEAR};

    //NLP nlp = NLP(problem, radau, mesh, guess);
    //std::cout << nlp.number_vars << std::endl;
    //std::cout << nlp.number_constraints << std::endl;
    std::shared_ptr<GDOP> a;
    std::shared_ptr<std::unordered_map<std::string, std::string>> b;
    IpoptSolver ipopt = IpoptSolver(a, b);
    ipopt.optimize();*/
    return 0;
}