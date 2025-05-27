#include <iomanip> 

#include "ipopt_adapter.h"


bool IpoptAdapter::get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) {
    n = nlp.number_vars;
    m = nlp.number_constraints;
    nnz_jac_g = nlp.nnz_jac;
    nnz_h_lag = nlp.nnz_hes;
    index_style = IndexStyleEnum::C_STYLE;
    return true;
};

bool IpoptAdapter::get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) {
    nlp.x_lb.write_to(x_l);
    nlp.x_ub.write_to(x_u);
    nlp.g_lb.write_to(g_l);
    nlp.g_ub.write_to(g_u);
    return true;
};

bool IpoptAdapter::get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L,
                                      Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) {
    assert(init_x      == true);
    assert(init_lambda == false);
    assert(init_z      == false);
    nlp.init_x.write_to(x);
    return true;
};

bool IpoptAdapter::get_scaling_parameters(Ipopt::Number& obj_scaling, bool& use_x_scaling, Ipopt::Index n, Ipopt::Number* x_scaling,
                                          bool& use_g_scaling, Ipopt::Index m, Ipopt::Number* g_scaling) {
    use_x_scaling = false;
    use_g_scaling = false;
    return true;
};

bool IpoptAdapter::eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) {
    nlp.eval_f(x, new_x);
    obj_value = nlp.curr_obj;
    return true;
};

bool IpoptAdapter::eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) {
    nlp.eval_grad_f(x, new_x);
    nlp.curr_grad.write_to(grad_f);
    return true;
};

bool IpoptAdapter::eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) {
    nlp.eval_g(x, new_x);
    nlp.curr_g.write_to(g);
    return true;
};

bool IpoptAdapter::eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac,
                              Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) {
    if (!values) {
        nlp.i_row_jac.write_to(iRow);
        nlp.j_col_jac.write_to(jCol);
    }
    else {
        nlp.eval_jac_g(x, new_x);
        nlp.curr_jac.write_to(values);
    }
    return true;
};

bool IpoptAdapter::eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
                          bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) {
    if (!values) {
        nlp.i_row_hes.write_to(iRow);
        nlp.j_col_hes.write_to(jCol);
    }
    else {
        nlp.eval_hes(x, lambda, obj_factor, new_x, new_lambda);
        nlp.curr_hes.write_to(values);
    }
    return true;
};

void IpoptAdapter::finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L,
                                     const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                                     Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) {
    // TODO: to implement
};

bool IpoptAdapter::intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size,
    Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) {
    // TODO: to implement
    return true;
};

