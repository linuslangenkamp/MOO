#ifndef OPT_NLP_H
#define OPT_NLP_H

#include "../base/fixed_vector.h"


// generic NLP base class - can be used with the generic NLP_Solver interface
// allows for a nice abstraction / to use any NLP solver when implementing this NLP
// also all NLP members must be implemented for the solvers anyway :)
// thus, GDOP is just an implementation of this NLP
class NLP {
public:
    NLP() = default;

    // NLP stuff itself
    int number_vars;         // total number of variables in the NLP
    int number_constraints;  // total number of constraints in the NLP
    int nnz_jac = 0;         // nnz Jacobian in the NLP
    int nnz_hes = 0;         // nnz Hessian in the NLP

    // current iterates
    FixedVector<double> init_x;       // initial NLP primal variables
    FixedVector<double> curr_x;       // current NLP primal variables
    FixedVector<double> curr_lambda;  // current NLP dual variables
    double              curr_sigma_f; // current objective weight in hessian

    // variable bounds
    FixedVector<double> x_lb; // lower bounds on NLP variables
    FixedVector<double> x_ub; // upper bounds on NLP variables

    // nlp function data
    double curr_obj;                // current NLP objective value
    FixedVector<double> curr_grad;  // current NLP gradient of the objective function
    FixedVector<double> curr_g;     // current NLP constraint function evaluation
    FixedVector<double> curr_jac;   // current NLP jacobian of the constraints
    FixedVector<double> der_jac;    // constant NLP derivative matrix part of the jacobian
    FixedVector<double> curr_hes;   // current NLP hessian of the lagrangian

    // scaled constraint bounds
    FixedVector<double> g_lb; // lower bounds on NLP constraints
    FixedVector<double> g_ub; // upper bounds on NLP constraints

    // COO sparsity patterns
    FixedVector<int> i_row_jac; // row COO of the Jacobian
    FixedVector<int> j_col_jac; // column COO of the Jacobian
    FixedVector<int> i_row_hes; // row COO of the Hessian
    FixedVector<int> j_col_hes; // column COO of the Hessian

    virtual void eval_f_safe(const double* nlp_solver_x, bool new_x) = 0;
    virtual void eval_g_safe(const double* nlp_solver_x, bool new_x) = 0;
    virtual void eval_grad_f_safe(const double* nlp_solver_x, bool new_x) = 0;
    virtual void eval_jac_g_safe(const double* nlp_solver_x, bool new_x) = 0;
    virtual void eval_hes_safe(const double* nlp_solver_x, const double* nlp_solver_lambda, double sigma, bool new_x, bool new_lambda) = 0;
};

#endif  // OPT_NLP_H
