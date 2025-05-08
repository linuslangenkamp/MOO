#include "test_problem_impl.h"
#include <math.h>
// a first cpp implementation of the GDOP Problem


FullSweepTestImpl::FullSweepTestImpl(FixedVector<FunctionLFG>&& lfg_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds)
    : FullSweep(std::move(lfg_in), mesh, g_bounds, true, 1, 0, 1, 0, 1) {
        // integral (x - p)²
        lfg[0].eval = &eval_buffer[0];
        lfg[0].jac.dx[0].value = &jac_buffer[0];
        lfg[0].jac.dp[0].value = &jac_buffer[1];
        lfg[0].hes.dx_dx[0].value = &hes_buffer[0];
        lfg[0].hes.dp_dx[0].value = &hes_buffer[1];
        lfg[0].hes.dp_dp[0].value = &hes_buffer[2];

        // T' = -0.4 * T + 120.2 + 0.2 * sin(4*pi*time)
        lfg[1].eval = &eval_buffer[1];
        lfg[1].jac.dx[0].value = &jac_buffer[2];
    }

void FullSweepTestImpl::callback_eval(const double* xu_nlp, const double* p) {
    int sum_ij = 0;
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            const double* xu_ij = xu_nlp + (x_size + u_size) * sum_ij;
            double T = xu_ij[0];
            eval_buffer[2 * sum_ij] = (T - p[0] - 0.015) * (T - p[0] - 0.015);
            eval_buffer[2 * sum_ij + 1] = -0.4 * T + 120.2 + 0.2 * sin(4*M_PI*mesh->t[i][j]);
            sum_ij++;
        }
    }
}

void FullSweepTestImpl::callback_jac(const double* xu_nlp, const double* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const double* xu_ij = xu_nlp + (x_size + u_size) * i;
        double T = xu_ij[0];
        jac_buffer[3 * i] = 2 * (T - p[0] - 0.015);
        jac_buffer[3 * i + 1] = -2 * (T - p[0] - 0.015);
        jac_buffer[3 * i + 2] = -0.4;
    }
}

void FullSweepTestImpl::callback_hes(const double* xu_nlp, const double* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        hes_buffer[3 * i] = 2;
        hes_buffer[3 * i + 1] = -2;
        hes_buffer[3 * i + 2] = 2;
    }
}

BoundarySweepTestImpl::BoundarySweepTestImpl(FixedVector<FunctionMR>&& mr_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds)
    : BoundarySweep(std::move(mr_in), mesh, r_bounds, false, 1, 1, 1) {
        mr[0].eval = &eval_buffer[0];
        mr[0].jac.dx0[0].value = &jac_buffer[0];
        mr[0].jac.dp[0].value = &jac_buffer[1];
    }

void BoundarySweepTestImpl::callback_eval(const double* x0_nlp, const double* xf_nlp, const double* p) {
    eval_buffer[0] = x0_nlp[0] - p[0];
}

void BoundarySweepTestImpl::callback_jac(const double* x0_nlp, const double* xf_nlp, const double* p) {
    jac_buffer[0] = 1;
    jac_buffer[1] = -1;
}

void BoundarySweepTestImpl::callback_hes(const double* x0_nlp, const double* xf_nlp, const double* p) {

};
