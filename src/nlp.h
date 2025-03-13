#ifndef OPT_NLP_H
#define OPT_NLP_H

#include "linalg.h"
#include "problem.h"
#include "collocation.h"
#include "util.h"

struct NLP {
    // NLP stuff itself
    uint32_t number_vars = 0;
    uint32_t number_constraints = 0;
    uint32_t nnz_jac = 0;
    uint32_t nnz_hes = 0;

    // current iterates
    std::unique_ptr<double> curr_x_scaled;
    std::unique_ptr<double> curr_x_unscaled;

    // scaling arras: iterate nominals
    std::unique_ptr<double> x_nominal; // shift scaling vector

    // nlp function data
    double curr_obj;
    std::unique_ptr<double> curr_grad;
    std::unique_ptr<double> curr_g;
    std::unique_ptr<double> curr_jac_values;
    std::unique_ptr<double> curr_hes_values;

    // use these later, fill one time and then scale at the end of calculations
    double curr_obj_nominal = 1;
    std::unique_ptr<double> curr_grad_nominal;
    std::unique_ptr<double> curr_g_nominal;
    std::unique_ptr<double> curr_jac_values_nominal;
    std::unique_ptr<double> curr_hes_values_nominal;

    // structures
    Problem problem;
    Collocation collocation;

    // evaluation data (problem double* point to these arrays)
    std::unique_ptr<double> eval_data_LFG;
    std::unique_ptr<double> jac_data_LFG_D; // contains the derivative matrix blocks, which must be calculated only once -> use memcpy _D -> stdjac
    std::unique_ptr<double> jac_data_LFG;
    std::unique_ptr<double> hes_data_LFG;

    uint32_t block_size_jac_LFG = 0;
    uint32_t block_size_hes_LFG = 0;

    std::unique_ptr<double> eval_data_MR;
    std::unique_ptr<double> jac_data_MR;
    std::unique_ptr<double> hes_data_MR;
        
    // scaling methods
    void scaleCurrentIterate();   // x_unscaled -> x_scaled
    void unscaleCurrentIterate(); // x_scaled   -> x_unscaled
    void scaleConstraintEvaluation();
    void scaleJacobianEvaluation();
    void scaleHessianEvaluation();
};

#endif  // OPT_NLP_H
