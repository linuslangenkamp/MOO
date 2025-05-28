#include "ipopt_solver.h"


IpoptSolver::IpoptSolver(NLP& nlp, NLPSolverFlags& solver_flags)
    : NLPSolver(nlp, solver_flags),
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
    app->Options()->SetIntegerValue("max_iter", solver_flags.get_flag_f64_fallback("Iterations", 5000));
    app->Options()->SetNumericValue("max_cpu_time", solver_flags.get_flag_f64_fallback("CPUTime", 3600));

    // numeric values
    app->Options()->SetNumericValue("tol", solver_flags.get_flag_f64_fallback("Tolerance", 1e-10));
    app->Options()->SetNumericValue("acceptable_tol", solver_flags.get_flag_f64_fallback("Tolerance", 1e-10) * 1e3);
    app->Options()->SetNumericValue("bound_push", 1e-2);
    app->Options()->SetNumericValue("bound_frac", 1e-2);
    app->Options()->SetNumericValue("alpha_red_factor", 0.5);

    // strategies
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("adaptive_mu_globalization", "kkt-error");
    app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    app->Options()->SetStringValue("fixed_variable_treatment", "make_parameter");
    app->Options()->SetStringValue("bound_mult_init_method","constant");
    app->Options()->SetStringValue("dependency_detection_with_rhs", "yes");
    app->Options()->SetNumericValue("nu_init", 1e-9);
    app->Options()->SetNumericValue("eta_phi", 1e-10);

    // Hessian approximation
    if (solver_flags.check_flag("Hessian", "LBFGS")) {
        app->Options()->SetStringValue("hessian_approximation", "limited-memory");
    }

    // subproblem
    app->Options()->SetStringValue("linear_solver", solver_flags.get_flag_string_fallback("LinearSolver", "MUMPS"));

    // constant derivatives
    app->Options()->SetStringValue("grad_f_constant", "no");
    app->Options()->SetStringValue("jac_c_constant", "no");
    app->Options()->SetStringValue("jac_d_constant", "no");
    app->Options()->SetStringValue("hessian_constant", "no");

    // info
    app->Options()->SetStringValue("timing_statistics", "yes");
    app->Options()->SetIntegerValue("print_level", 5);

    // testing + validation
    if (solver_flags.check_flag("Ipopt_DerivativeTest", "true")) {
        app->Options()->SetStringValue("derivative_test", "second-order");
        app->Options()->SetNumericValue("derivative_test_tol", 1e-2);
        app->Options()->SetNumericValue("point_perturbation_radius", 0);
    }
}
