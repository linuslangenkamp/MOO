
#include "adapter.h"

#include "solver.h"

namespace IpoptSolver {

struct IpoptSolverData {
    Ipopt::SmartPtr<IpoptAdapter> adapter;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app;

    IpoptSolverData(NLP::NLP& nlp) : adapter(new IpoptAdapter(nlp)), app(IpoptApplicationFactory()) {}
};

IpoptSolver::IpoptSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings)
    : NLPSolver(nlp, solver_settings), data{new IpoptSolverData{nlp}}
{
    init_application();
}

IpoptSolver::~IpoptSolver() {
    delete data;
}

// simple wrapper to adapter
void IpoptSolver::optimize() {
    set_settings();

    Ipopt::ApplicationReturnStatus status = data->app->OptimizeTNLP(data->adapter);

    switch (status) {
        case Ipopt::Solve_Succeeded:
            LOG_SUCCESS("[Ipopt Interface] Optimization succeeded!");
            break;

        case Ipopt::Solved_To_Acceptable_Level:
            LOG_SUCCESS("[Ipopt Interface] Optimization succeeded (acceptable)!");
            break;

        case Ipopt::Infeasible_Problem_Detected:
            LOG_ERROR("[Ipopt Interface] Infeasible problem detected.");
            break;

        case Ipopt::Search_Direction_Becomes_Too_Small:
            LOG_WARNING("[Ipopt Interface] Search direction became too small.");
            break;

        case Ipopt::Diverging_Iterates:
            LOG_ERROR("[Ipopt Interface] Diverging iterates.");
            break;

        case Ipopt::User_Requested_Stop:
            LOG_WARNING("[Ipopt Interface] Optimization stopped by user request.");
            break;

        case Ipopt::Feasible_Point_Found:
            LOG_ERROR("[Ipopt Interface] Feasible point found.");
            break;

        case Ipopt::Maximum_Iterations_Exceeded:
            LOG_WARNING("[Ipopt Interface] Maximum iterations exceeded.");
            break;

        case Ipopt::Restoration_Failed:
            LOG_ERROR("[Ipopt Interface] Restoration failed.");
            break;

        case Ipopt::Error_In_Step_Computation:
            LOG_ERROR("[Ipopt Interface] Error in step computation.");
            break;

        case Ipopt::Maximum_CpuTime_Exceeded:
            LOG_WARNING("[Ipopt Interface] Maximum CPU time exceeded.");
            break;

        case Ipopt::Maximum_WallTime_Exceeded:
            LOG_WARNING("[Ipopt Interface] Maximum wall time exceeded.");
            break;

        case Ipopt::Not_Enough_Degrees_Of_Freedom:
            LOG_ERROR("[Ipopt Interface] Not enough degrees of freedom.");
            break;

        case Ipopt::Invalid_Problem_Definition:
            LOG_ERROR("[Ipopt Interface] Invalid problem definition.");
            break;

        case Ipopt::Invalid_Option:
            LOG_ERROR("[Ipopt Interface] Invalid option.");
            break;

        case Ipopt::Invalid_Number_Detected:
            LOG_ERROR("[Ipopt Interface] Invalid number detected.");
            break;

        case Ipopt::Unrecoverable_Exception:
            LOG_ERROR("[Ipopt Interface] Unrecoverable exception occurred.");
            break;

        case Ipopt::NonIpopt_Exception_Thrown:
            LOG_ERROR("[Ipopt Interface] Non-Ipopt exception thrown.");
            break;

        case Ipopt::Insufficient_Memory:
            LOG_ERROR("[Ipopt Interface] Insufficient memory.");
            break;

        case Ipopt::Internal_Error:
            LOG_ERROR("[Ipopt Interface] Internal error.");
            break;

        default:
            LOG_ERROR("[Ipopt Interface] Unknown return status: {}", static_cast<int>(status));
            break;
    }
}

void IpoptSolver::init_application() {
    Ipopt::ApplicationReturnStatus status = data->app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        LOG_ERROR("[Ipopt Interface] Error during application initialization!");
        abort();
    }
}

void IpoptSolver::set_settings() {
    // --- termination ---
    data->app->Options()->SetIntegerValue("max_iter", solver_settings.get_or_default<int>(NLP::Option::Iterations));
    data->app->Options()->SetNumericValue("max_cpu_time", solver_settings.get_or_default<f64>(NLP::Option::CPUTime));

    // --- tolerances, step sizes ---
    f64 tol = solver_settings.get_or_default<f64>(NLP::Option::Tolerance);
    data->app->Options()->SetNumericValue("tol", tol);
    data->app->Options()->SetNumericValue("acceptable_tol", tol * 1e3);
    data->app->Options()->SetNumericValue("bound_push", 1e-2);
    data->app->Options()->SetNumericValue("bound_frac", 1e-2);
    data->app->Options()->SetNumericValue("alpha_red_factor", 0.5);

    // --- strategy settings ---
    data->app->Options()->SetStringValue("mu_strategy", "adaptive");
    data->app->Options()->SetStringValue("adaptive_mu_globalization", "kkt-error");
    data->app->Options()->SetStringValue("nlp_scaling_method", "gradient-based");
    data->app->Options()->SetStringValue("fixed_variable_treatment", "make_parameter");

    // --- hessian options ---
    NLP::HessianOption hess_opt = solver_settings.get_or_default<NLP::HessianOption>(NLP::Option::Hessian);
    switch (hess_opt) {
        case NLP::HessianOption::LBFGS:
            data->app->Options()->SetStringValue("hessian_approximation", "limited-memory");
            break;
        case NLP::HessianOption::CONST:
            data->app->Options()->SetStringValue("hessian_constant", "yes");
            break;
        case NLP::HessianOption::Exact:
            data->app->Options()->SetStringValue("hessian_approximation", "exact");
            break;
    }

    // --- warm start ---
    if (solver_settings.option_is_true(NLP::Option::WarmStart)) {
        data->app->Options()->SetStringValue("warm_start_init_point", "yes");
        data->app->Options()->SetStringValue("mu_strategy", "monotone");
        data->app->Options()->SetNumericValue("mu_init", 1e-14);
        data->app->Options()->SetNumericValue("warm_start_bound_push", 1e-8);
        data->app->Options()->SetNumericValue("warm_start_bound_frac", 1e-8);
        data->app->Options()->SetNumericValue("warm_start_slack_bound_push", 1e-8);
        data->app->Options()->SetNumericValue("warm_start_slack_bound_frac", 1e-8);
        data->app->Options()->SetNumericValue("warm_start_mult_bound_push", 1e-8);
    }

    // --- linear solver ---
    NLP::LinearSolverOption linear_solver = solver_settings.get_or_default<NLP::LinearSolverOption>(NLP::Option::LinearSolver);
    switch (linear_solver) {
        case NLP::LinearSolverOption::MUMPS: data->app->Options()->SetStringValue("linear_solver", "mumps"); break;
        case NLP::LinearSolverOption::MA27:  data->app->Options()->SetStringValue("linear_solver", "ma27");  break;
        case NLP::LinearSolverOption::MA57:  data->app->Options()->SetStringValue("linear_solver", "ma57");  break;
        case NLP::LinearSolverOption::MA77:  data->app->Options()->SetStringValue("linear_solver", "ma77");  break;
        case NLP::LinearSolverOption::MA86:  data->app->Options()->SetStringValue("linear_solver", "ma86");  break;
        case NLP::LinearSolverOption::MA97:  data->app->Options()->SetStringValue("linear_solver", "ma97");  break;
    }

    // --- constant derivatives (assumed false for now) ---
    data->app->Options()->SetStringValue("grad_f_constant", "no");
    data->app->Options()->SetStringValue("jac_c_constant", "no");
    data->app->Options()->SetStringValue("jac_d_constant", "no");

    // --- info ---
    data->app->Options()->SetStringValue("timing_statistics", "yes");
    data->app->Options()->SetIntegerValue("print_level", 5);

    // --- derivative test (optional) ---
    if (solver_settings.option_is_true(NLP::Option::IpoptDerivativeTest)) {
        data->app->Options()->SetStringValue("derivative_test", "second-order");
        data->app->Options()->SetNumericValue("derivative_test_tol", 1e-2);
        data->app->Options()->SetNumericValue("point_perturbation_radius", 0);
    }
}

} // namespace IpoptSolver
