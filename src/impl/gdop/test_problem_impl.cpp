#include "test_problem_impl.h"

// a first cpp implementation of the GDOP Problem


FullSweepTestImpl::FullSweepTestImpl(FixedVector<FunctionLFG>& lfg_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds)
    : FullSweep(lfg_in, mesh, g_bounds, false, 1, 0, 1, 1, 0) {
        lfg[0].eval = &eval_buffer[0];
        lfg[0].jac.dx[0].value = &jac_buffer[0];
        lfg[0].jac.du[0].value = &jac_buffer[1];
        lfg[0].hes.dx_dx[0].value = &hes_buffer[0];
    }

void FullSweepTestImpl::callbackEval(const double* xu_nlp, const double* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const double* xu_ij = xu_nlp + (x_size + u_size) * i;
        eval_buffer[i] = cos(xu_ij[0]) + xu_ij[1]; 
    }
}

void FullSweepTestImpl::callbackJac(const double* xu_nlp, const double* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const double* xu_ij = xu_nlp + (x_size + u_size) * i;
        jac_buffer[jac_size * i] = -sin(xu_ij[0]);
        jac_buffer[jac_size * i + 1] = 1;
    }
}

void FullSweepTestImpl::callbackHes(const double* xu_nlp, const double* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const double* xu_ij = xu_nlp + (x_size + u_size) * i;
        hes_buffer[hes_size * i] = -cos(xu_ij[0]);
    }
}

BoundarySweepTestImpl::BoundarySweepTestImpl(FixedVector<FunctionMR>& mr_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds)
    : BoundarySweep(mr_in, mesh, r_bounds, true, 0, 1, 0) {
        mr[0].eval = &eval_buffer[0];
        mr[0].jac.dxf[0].value = &jac_buffer[0];
        mr[0].hes.dxf_dxf[0].value = &hes_buffer[0];
    }

void BoundarySweepTestImpl::callbackEval(const double* x0_nlp, const double* xf_nlp, const double* p) {
    eval_buffer[0] = xf_nlp[0] * xf_nlp[0];
}

void BoundarySweepTestImpl::callbackJac(const double* x0_nlp, const double* xf_nlp, const double* p) {
    jac_buffer[0] = 2 * xf_nlp[0];
}

void BoundarySweepTestImpl::callbackHes(const double* x0_nlp, const double* xf_nlp, const double* p) {
    hes_buffer[0] = 2;
};

