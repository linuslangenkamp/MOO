#include "ipopt_solver.h"


IpoptSolver::IpoptSolver(std::shared_ptr<NLP> nlp, std::shared_ptr<std::unordered_map<std::string, std::string>> solver_settings)
    : NLPSolver(nlp, solver_settings),
      app(IpoptApplicationFactory()),
      adapter(new IpoptAdapter(nlp)) {
    initIpoptApplication();
}

// simple wrapper to adapter
void IpoptSolver::optimize() {
    Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(adapter);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << "[Ipopt Interface] Optimization with succeeded!" << std::endl;
    } else {
        std::cout << "[Ipopt Interface] Optimization failed with status: " << status << std::endl;
    }
}

void IpoptSolver::initIpoptApplication() {
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << "[Ipopt Interface] Error during application initialization!" << std::endl;
        abort();
    }

    // set all the settings here
    app->Options()->SetNumericValue("tol", 1e-12);
    app->Options()->SetNumericValue("bound_push", 1e-2);
    app->Options()->SetNumericValue("bound_frac", 1e-2);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("fixed_variable_treatment", "make_parameter");
    app->Options()->SetIntegerValue("max_iter", 15);
    app->Options()->SetStringValue("linear_solver", "MUMPS");
    app->Options()->SetStringValue("timing_statistics", "yes");
}