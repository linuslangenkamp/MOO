#include "ipopt_solver.h"


IpoptSolver::IpoptSolver(std::shared_ptr<NLP> nlp, std::shared_ptr<std::unordered_map<std::string, std::string>> solver_settings)
    : NLPSolver(nlp, solver_settings),
      adapter(new IpoptAdapter(nlp)),
      app(IpoptApplicationFactory()) {
    init_IpoptApplication();
}

// simple wrapper to adapter
void IpoptSolver::optimize() {
    Ipopt::ApplicationReturnStatus status = app->OptimizeTNLP(adapter);

    if (status == Ipopt::Solve_Succeeded) {
        std::cout << "\n[Ipopt Interface] Optimization succeeded!" << std::endl;
    } else {
        std::cout << "[Ipopt Interface] Optimization failed with status: " << status << std::endl;
    }
}

void IpoptSolver::init_IpoptApplication() {
    Ipopt::ApplicationReturnStatus status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        std::cout << "[Ipopt Interface] Error during application initialization!" << std::endl;
        abort();
    }

    // set all the settings here

    // termination fallback
    app->Options()->SetIntegerValue("max_iter", 250);
    app->Options()->SetNumericValue("max_cpu_time", 3600);

    // numeric values
    app->Options()->SetNumericValue("tol", 1e-10);
    app->Options()->SetNumericValue("acceptable_tol", 1e-9);
    app->Options()->SetNumericValue("bound_push", 1e-2);
    app->Options()->SetNumericValue("bound_frac", 1e-2);
    app->Options()->SetNumericValue("alpha_red_factor", 0.6);
    
    // strategies
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("adaptive_mu_globalization", "kkt-error");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    app->Options()->SetStringValue("fixed_variable_treatment", "make_parameter");
    
    // subproblem
    app->Options()->SetStringValue("linear_solver", "MUMPS");

    // constant derivatives
    app->Options()->SetStringValue("grad_f_constant", "no");
    app->Options()->SetStringValue("jac_c_constant", "no");
    app->Options()->SetStringValue("jac_d_constant", "no");
    app->Options()->SetStringValue("hessian_constant", "no");

    // info
    app->Options()->SetStringValue("timing_statistics", "yes");
    app->Options()->SetIntegerValue("print_level", 5);

    // testings
    // app->Options()->SetStringValue("derivative_test", "second-order");
    // app->Options()->SetStringValue("hessian_approximation", "limited-memory");
}