#ifndef OPT_NLP_H
#define OPT_NLP_H

#include <base/fixed_vector.h>

// TODO: maybe add the scaler to the NLP and only do the setup in the generic NLP implementation, then the scaler is always called on the iterates

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
    FixedVector<F64> init_x;       // initial NLP primal variables
    FixedVector<F64> curr_x;       // current NLP primal variables
    FixedVector<F64> curr_lambda;  // current NLP dual variables
    F64              curr_sigma_f; // current objective weight in hessian

    // variable bounds
    FixedVector<F64> x_lb; // lower bounds on NLP variables
    FixedVector<F64> x_ub; // upper bounds on NLP variables

    // nlp function data
    F64 curr_obj;                   // current NLP objective value
    FixedVector<F64> curr_grad;     // current NLP gradient of the objective function
    FixedVector<F64> curr_g;        // current NLP constraint function evaluation
    FixedVector<F64> curr_jac;      // current NLP jacobian of the constraints
    FixedVector<F64> curr_hes;      // current NLP hessian of the lagrangian

    // scaled constraint bounds
    FixedVector<F64> g_lb; // lower bounds on NLP constraints
    FixedVector<F64> g_ub; // upper bounds on NLP constraints

    // COO sparsity patterns
    FixedVector<int> i_row_jac; // row COO of the Jacobian
    FixedVector<int> j_col_jac; // column COO of the Jacobian
    FixedVector<int> i_row_hes; // row COO of the Hessian
    FixedVector<int> j_col_hes; // column COO of the Hessian
    
    // TODO: add a generic block BFGS routine, which can calculate blocks of the Lagrangian Hessian $\nabla_{xx} \mathcal{L}_{AA -> BB}$
    
    virtual void eval_f(const F64* nlp_solver_x, bool new_x) = 0;      // fill curr_obj
    virtual void eval_g(const F64* nlp_solver_x, bool new_x) = 0;      // fill curr_g
    virtual void eval_grad_f(const F64* nlp_solver_x, bool new_x) = 0; // fill curr_grad
    virtual void eval_jac_g(const F64* nlp_solver_x, bool new_x) = 0;  // fill curr_jac
    virtual void eval_hes(const F64* nlp_solver_x, const F64* nlp_solver_lambda, F64 sigma, bool new_x, bool new_lambda) = 0; // fill curr_hes
};

#endif  // OPT_NLP_H
