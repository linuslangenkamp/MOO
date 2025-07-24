#ifndef OPT_NLP_H
#define OPT_NLP_H

#include <base/fixed_vector.h>

#include "nlp_scaling.h"

namespace NLP {

// TODO: maybe design this interface, such that we provide the pointers like in Ipopt, so its clear what needs to be done in which function
//       will make the generic interface way more modular and sane :: make stuff private

/* generic NLP base class - can be used with the generic NLP_Solver interface
 * allows for a nice abstraction / to use any NLP solver when implementing this NLP
 * 
 *    min f(x)
 * s.t.
 *    g^{LB} <= g(x) <= g^{UB}
 *    x^{LB} <=  x   <= x^{UB}
 * 
 * Callbacks:
 *    eval_f()      => f(x)
 *    eval_g()      => g(x)
 *    eval_grad_f() => df/dx
 *    eval_jac_g()  => dg/dx
 *    eval_hes()    => sigma_f * d²f/dx² + sum_{i=1}^{m} lambda_i * d²(g_i)/dx²
 */
class NLP {
public:
    NLP() = default;

    // NLP stuff itself
    int number_vars = 0;        // total number of variables in the NLP
    int number_constraints = 0; // total number of constraints in the NLP
    int nnz_jac = 0;            // nnz Jacobian in the NLP
    int nnz_hes = 0;            // nnz Hessian in the NLP

    // current iterates
    FixedVector<f64> curr_x;       // current NLP primal variables
    FixedVector<f64> curr_lambda;  // current NLP dual variables
    f64              curr_sigma_f; // current objective weight in hessian

    // variable bounds
    FixedVector<f64> x_lb; // lower bounds on NLP variables
    FixedVector<f64> x_ub; // upper bounds on NLP variables

    // optimal dual multipliers (set after finalize)
    FixedVector<f64> z_lb; // dual multipliers of lower bounds on NLP variables
    FixedVector<f64> z_ub; // dual multipliers of upper bounds on NLP variables

    // nlp function data
    f64 curr_obj;               // current NLP objective value
    FixedVector<f64> curr_grad; // current NLP gradient of the objective function
    FixedVector<f64> curr_g;    // current NLP constraint function evaluation
    FixedVector<f64> curr_jac;  // current NLP jacobian of the constraints
    FixedVector<f64> curr_hes;  // current NLP hessian of the lagrangian

    // scaled constraint bounds
    FixedVector<f64> g_lb; // lower bounds on NLP constraints
    FixedVector<f64> g_ub; // upper bounds on NLP constraints

    // COO sparsity patterns
    FixedVector<int> i_row_jac; // row COO of the Jacobian
    FixedVector<int> j_col_jac; // column COO of the Jacobian
    FixedVector<int> i_row_hes; // row COO of the Hessian
    FixedVector<int> j_col_hes; // column COO of the Hessian

    // TODO: add a generic block BFGS routine, which can calculate blocks of the Lagrangian Hessian $\nabla_{xx} \mathcal{L}_{AA -> BB}$

    std::shared_ptr<Scaling> scaling = std::make_shared<NoScaling>();                 // generic scaling routine
    void set_scaling(std::shared_ptr<Scaling> new_scaling) { scaling = new_scaling; } // set the generic scaling routine

    virtual void eval_f(bool new_x) = 0;                    // fill curr_obj
    virtual void eval_g(bool new_x) = 0;                    // fill curr_g
    virtual void eval_grad_f(bool new_x) = 0;               // fill curr_grad
    virtual void eval_jac_g(bool new_x) = 0;                // fill curr_jac
    virtual void eval_hes(bool new_x, bool new_lambda) = 0; // fill curr_hes
    virtual void finalize_solution() = 0;                   // finalize the solution

    inline void update_unscale_dual_bounds(const f64* z_L, const f64* z_U) {
        scaling->unscale_x(z_L, z_lb.raw(), number_vars);
        scaling->unscale_x(z_U, z_ub.raw(), number_vars);
    }

    inline void update_unscale_curr_x(bool new_x, const f64* x) {
        if (new_x) scaling->unscale_x(x, curr_x.raw(), number_vars);
    }

    inline void update_unscale_curr_lambda(bool new_lambda, const f64* lambda) {
        if (new_lambda) scaling->scale_g(lambda, curr_lambda.raw(), number_constraints);
    }

    inline void update_unscale_curr_sigma_f(const f64* sigma_f) {
        scaling->scale_f(sigma_f, &curr_sigma_f);
    }

    inline void unscale_objective(const f64* obj) {
        scaling->unscale_f(obj, &curr_obj);
    }
};

} // namespace NLP

#endif  // OPT_NLP_H
