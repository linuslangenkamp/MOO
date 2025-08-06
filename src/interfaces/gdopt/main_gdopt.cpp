#include <src/nlp/instances/gdop/orchestrator.h>
#include <src/nlp/solvers/ipopt/solver.h>
#include <src/nlp/instances/gdop/strategies.h>
#include <src/interfaces/c/structures.h>
#include <src/interfaces/c/problem.h>
#include <src/base/log.h>

extern f64* _rp;


// create config for the algorithm (for now basic)
class Config {

};

Config read_yaml() {
    return Config();
}

void set_global_configuration(Config& config) {

}

void set_global_runtime_parameters(Config& config) {
    /* set _rp */
}

int main_gdopt(int argc, char** argv) {
    LOG_PREFIX('*', "Entry point [OPT] - _main_OptimizationRuntime\n");

    auto config = read_yaml();
    set_global_configuration(config);

    set_global_runtime_parameters(config);

    c_problem_t* c_problem = get_update_c_problem();

    auto nlp_solver_settings = NLP::NLPSolverSettings(argc, argv);
    //nlp_solver_settings.set(NLP::Option::IpoptDerivativeTest, true);
    nlp_solver_settings.print();

    auto mesh = Mesh::create_equidistant_fixed_stages(1 /* tf */, 25 /* intervals */, 3 /* stages */);
    auto problem = C::create_gdop(c_problem, *mesh);

    // auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    auto strategies = std::make_unique<GDOP::Strategies>(GDOP::Strategies::default_strategies());
    auto gdop = GDOP::GDOP(problem);

    IpoptSolver::IpoptSolver ipopt_solver(gdop, nlp_solver_settings);

    auto orchestrator = GDOP::MeshRefinementOrchestrator(gdop, std::move(strategies), ipopt_solver);

    orchestrator.optimize();

    return 0;
}
