#ifndef OPT_NLP_H
#define OPT_NLP_H

#include "linalg.h"
#include "problem.h"
#include "mesh.h"
#include "collocation.h"
#include "util.h"

struct NLP {
    NLP(Problem problem, Collocation collocation, Mesh mesh)
        : problem(std::make_shared<Problem>(problem)),
          collocation(std::make_shared<Collocation>(collocation)),
          mesh(std::make_shared<Mesh>(mesh)) {
        init();
    }

    // NLP stuff itself
    int off_x;  // offset #xVars
    int off_u;  // offset #uVars
    int off_p;  // offset #pVars
    int off_xu; // number of vars for one collocation grid point

    // note off_acc_xu[0][0] = off_x for time t = first collocation node, since there are no controls at time t=0
    std::vector<std::vector<int>> off_acc_xu;  // offset to *(x, u)(t_ij): NLP_X[off_acc_xu[i][j]],
    int off_xu_total;                          // first parameter variable
    int number_vars;                           // total number of variables in the NLP
    int number_constraints = 0;
    int nnz_jac = 0;
    int nnz_hes = 0;
     
    // current iterates
    std::vector<double> curr_x_scaled;
    //std::vector<double> curr_x_unscaled;

    // scaled variable bounds
    std::vector<double> x_bounds_u;
    std::vector<double> x_bounds_l;

    // iterate nominals
    //std::vector<double> x_nominal; // shift scaling vector

    // nlp function data
    double curr_obj;
    std::vector<double> curr_grad;
    std::vector<double> curr_g;
    std::vector<double> curr_jac_values;
    std::vector<double> curr_hes_values;

    // scaled constraint bounds
    std::vector<double> g_bounds_u;
    std::vector<double> g_bounds_l;

    // COO sparsity patterns
    std::vector<int> i_row_jac;
    std::vector<int> j_col_jac;
    std::vector<int> i_row_hes;
    std::vector<int> j_col_hes;

    // structures
    std::shared_ptr<Mesh> mesh;
    std::shared_ptr<Problem> problem;
    std::shared_ptr<Collocation> collocation;

    // evaluation data (problem double* point to these arrays)
    std::vector<double> eval_data_LFG;
    std::vector<double> jac_data_LFG_D; // contains the derivative matrix blocks, which must be calculated only once -> use memcpy _D -> stdjac
    std::vector<double> jac_data_LFG;
    std::vector<double> hes_data_LFG;

    // offset for iterating over the data array
    int nnz_jac_LFG = 0;
    int nnz_hes_LFG = 0;

    std::vector<double> eval_data_MR;
    std::vector<double> jac_data_MR;
    std::vector<double> hes_data_MR;

    void init();
    void initSizesOffsets();

    /* TODO: add external scaler class which can perform, no, nominal, adaptive scaling
    // TODO: use these later, fill one time and then scale at the end of calculations
    // these have the same sizes as the curr_'s, just divide element wise
    double curr_obj_nominal = 1;
    std::vector<double> curr_grad_nominal;
    std::vector<double> curr_g_nominal;
    std::vector<double> curr_jac_values_nominal;
    std::vector<double> curr_hes_values_nominal;
    */
};

#endif  // OPT_NLP_H
