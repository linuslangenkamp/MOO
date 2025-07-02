#ifndef OPT_IPOPT_ADAPTER_H
#define OPT_IPOPT_ADAPTER_H

#include <IpTNLP.hpp>
#include <IpIpoptData.hpp>
#include <IpDenseVector.hpp>
#include <IpSmartPtr.hpp>
#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"

#include <nlp/nlp.h>


class IpoptAdapter : public Ipopt::TNLP {
public:
    explicit IpoptAdapter(NLP& nlp) : nlp(nlp) {}

    NLP& nlp;

    // Attention: this part is super confusing, since the lambda and sigma_f *unscale* is actually a scale with g and f!!
    // this is the case, because the augmented Lagrangian does not allow for to scale f and g independently after evaluation
    // so the only possibility is to apply this scaling a priori by updating the duals and sigma_f
    inline void update_unscale_curr_x(bool new_x, const f64* x) {
        if (new_x) nlp.scaling->unscale_x(x, nlp.curr_x.raw(), nlp.number_vars);
    }

    inline void update_unscale_curr_lambda(bool new_lambda, const f64* lambda) {
        if (new_lambda) nlp.scaling->scale_g(lambda, nlp.curr_lambda.raw(), nlp.number_constraints);
    }

    inline void update_unscale_curr_sigma_f(const f64* sigma_f) {
        nlp.scaling->scale_f(sigma_f, &nlp.curr_sigma_f);
    }

    bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style) override;

    bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u) override;

    bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda) override;

    bool get_scaling_parameters(Ipopt::Number& obj_scaling, bool& use_x_scaling, Ipopt::Index n, Ipopt::Number* x_scaling, bool& use_g_scaling, Ipopt::Index m, Ipopt::Number* g_scaling) override;

    bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value) override;

    bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f) override;

    bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g) override;

    bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values) override;

    bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
                Ipopt::Index* jCol, Ipopt::Number* values) override;

    void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                           Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq) override;

    bool intermediate_callback(Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value, Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu, Ipopt::Number d_norm, Ipopt::Number regularization_size,
     Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq);
};

#endif // OPT_IPOPT_ADAPTER_H
