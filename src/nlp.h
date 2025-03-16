#ifndef OPT_NLP_H
#define OPT_NLP_H

#include "linalg.h"
#include "problem.h"
#include "mesh.h"
#include "collocation.h"
#include "util.h"
#include "nlp_state.h"

struct NLP {
    NLP(Problem& problem, Collocation& collocation, Mesh& mesh, Trajectory& guess)
        : problem(std::make_shared<Problem>(problem)),
          collocation(std::make_shared<Collocation>(collocation)),
          mesh(std::make_shared<Mesh>(mesh)),
          guess(std::make_shared<Trajectory>(guess)) {
        init();
    }

    // NLP stuff itself
    int number_vars;         // total number of variables in the NLP
    int number_constraints;  // total number of constraints in the NLP
    int nnz_jac = 0;         // nnz Jacobian in the NLP
    int nnz_hes = 0;         // nnz Hessian in the NLP

    // current iterates
    std::unique_ptr<double[]> curr_x;

    // scaled variable bounds
    std::unique_ptr<double[]> x_lb;
    std::unique_ptr<double[]> x_ub;

    // nlp function data
    double curr_obj;                      // current objective value
    std::unique_ptr<double[]> curr_grad;  // current gradient of the objective function
    std::unique_ptr<double[]> curr_g;     // current constraint function evaluation
    std::unique_ptr<double[]> curr_jac;   // current jacobian of the constraints
    std::unique_ptr<double[]> der_jac;    // constant derivative matrix part of the jacobian
    std::unique_ptr<double[]> curr_hes;   // current hessian of the lagrangian

    // scaled constraint bounds
    std::unique_ptr<double[]> g_lb;
    std::unique_ptr<double[]> g_ub;

    // COO sparsity patterns
    std::unique_ptr<int[]> i_row_jac;
    std::unique_ptr<int[]> j_col_jac;
    std::unique_ptr<int[]> i_row_hes;
    std::unique_ptr<int[]> j_col_hes;

    // structures
    std::shared_ptr<Mesh> mesh;                // grid / mesh
    std::shared_ptr<Problem> problem;          // continuous GDOP
    std::shared_ptr<Collocation> collocation;  // collocation data
    std::shared_ptr<Trajectory> guess;         // initial guess / trajectory, will be interpolated accordingly
    NLP_State evaluation_state;                // simple state to check which actions are / have to be performed for an iteration

    // offsets
    int off_x;         // offset #xVars
    int off_u;         // offset #uVars
    int off_p;         // offset #pVars
    int off_xu;        // number of vars for one collocation grid point
    int off_last_xu;   // last collocation grid point *x_nm, u_nm
    int off_xu_total;  // first parameter variable index
    int off_fg_total;  // first boundary constraint index

    // note off_acc_xu[0][0] = off_x for time t = first collocation node, since there are no controls at time t=0
    std::vector<std::vector<int>> off_acc_xu;  // offset to NLP_X first index of (x, u)(t_ij), i.e. NLP_X[off_acc_xu[i][j]] = x[i][j], u[i][j]
    std::vector<std::vector<int>> off_acc_fg;  // offset to NLP_G first index of (f, g)(t_ij), i.e. NLP_G[off_acc_fg[i][j]] = f[i][j], g[i][j]
    std::vector<int> off_acc_jac;              // offset to NLP_JAC_G first index of nabla (f, g)(t_ij)

    /* TODO: maybe set these in fullsweep or boundarysweep, they dont need to be here, alloc them there w.r.t. macros in nlp.cpp!
    // evaluation data (problem double* point to these arrays)
    std::unique_ptr<double[]> eval_data_LFG;
    std::unique_ptr<double[]> jac_data_LFG;
    std::unique_ptr<double[]> hes_data_LFG;

    std::unique_ptr<double[]> eval_data_MR;
    std::unique_ptr<double[]>  jac_data_MR;
    std::unique_ptr<double[]>  hes_data_MR;
    */

    void init();
    void initSizesOffsets();
    void initBounds();
    void initStartingPoint();
    void initJacobian();
    void calculateJacobianNonzeros();
    void initJacobianSparsityPattern();
    void initHessianSparsity();

    // get callback data
    void callback_evaluation();
    void callback_jacobian();
    void callback_hessian();

    // nlp solver calls
    void check_new_x(const double* nlp_solver_x, bool new_x);
    void eval_f();
    void eval_f_safe(const double* nlp_solver_x, bool new_x);
    void eval_g();
    void eval_g_safe(const double* nlp_solver_x, bool new_x);
    void eval_grad_f();
    void eval_grad_f_safe(const double* nlp_solver_x, bool new_x);

    /* TODO: add external scaler class which can perform, no, nominal, adaptive scaling
    // TODO: use these later, fill one time and then scale at the end of calculations
    // these have the same sizes as the curr_'s, just divide element wise
    //std::vector<double> curr_x_unscaled;
    double curr_obj_nominal = 1;
    std::vector<double> curr_grad_nominal;
    std::vector<double> curr_g_nominal;
    std::vector<double> curr_jac_values_nominal;
    std::vector<double> curr_hes_values_nominal;
    */
};

#endif  // OPT_NLP_H
