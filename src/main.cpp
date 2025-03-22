#include <iostream>

#include "collocation.h"
#include "gdop.h"
#include "problem.h"
#include "interfaces/ipopt_solver.h"

//  cmake --build build && cd build && ./gdopt_experimental && cd ..

int main() {
    FixedField<double, 3> x(1, 2, 5);
    FixedField<double, 1> q(2);
    std::vector<double> vec = {7, 8, 9, 10, 11};
    FixedVector<double> z(vec.begin(), vec.end());
    q.assign(&vec[2], 2);
    std::cout << q[1] << std::endl;
    q.print();
    Collocation radau = Collocation();
    Mesh mesh = Mesh::createEquidistantMeshFixedDegree(10, 1, 3);
    BlockSparsity::createLowerTriangular(5, BlockType::Offset);
    Problem problem = Problem();
    problem.x_size = 1;
    problem.x_bounds = {{2, 4}};
    problem.x0_fixed = {std::nullopt};
    problem.xf_fixed = {std::nullopt};
    problem.u_size = 2;
    problem.u_bounds = {{2, 4}, {-1, 2}};
    problem.p_size = 0;
    problem.full.f_size = 1;
    problem.full.g_size = 1;
    problem.full.fg_size = 2;
    problem.full.g_bounds = {{1, 2}};
    problem.boundary.r_size = 1;
    problem.boundary.r_bounds = {{125, 2555}};
    Trajectory guess = {{0, 1}, {{1, 5}}, {{0, 2}, {3, 1}}, {1, 2}, InterpolationMethod::LINEAR};

    //NLP nlp = NLP(problem, radau, mesh, guess);
    //std::cout << nlp.number_vars << std::endl;
    //std::cout << nlp.number_constraints << std::endl;
    std::shared_ptr<GDOP> a;
    std::shared_ptr<std::unordered_map<std::string, std::string>> b;
    IpoptSolver ipopt = IpoptSolver(a, b);
    ipopt.optimize();
    return 0;
}