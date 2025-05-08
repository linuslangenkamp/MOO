#include "test_problem_impl.h"

// a first cpp implementation of the GDOP Problem


FullSweepTestImpl::FullSweepTestImpl(FixedVector<FunctionLFG>&& lfg_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds)
    : FullSweep(std::move(lfg_in), mesh, g_bounds, false, 2, 0, 2, 1, 0) {
        // x1' = -(u + u² / 2) * x1
        lfg[0].eval = &eval_buffer[0];
        lfg[0].jac.dx[0].value = &jac_buffer[0];
        lfg[0].jac.du[0].value = &jac_buffer[1];
        lfg[0].hes.du_dx[0].value = &hes_buffer[0];
        lfg[0].hes.du_du[0].value = &hes_buffer[1];

        // x2' = u * x1
        lfg[1].eval = &eval_buffer[1];
        lfg[1].jac.dx[0].value = &jac_buffer[2];
        lfg[1].jac.du[0].value = &jac_buffer[3];
        lfg[1].hes.du_dx[0].value = &hes_buffer[2];
    }

void FullSweepTestImpl::callback_eval(const f64* xu_nlp, const f64* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const f64* xu_ij = xu_nlp + (x_size + u_size) * i;
        f64 x1 = xu_ij[0];
        f64 u = xu_ij[2];
        eval_buffer[2 * i] = -(u + u*u / 2) * x1;
        eval_buffer[2 * i + 1] = u * x1;
    }
}

void FullSweepTestImpl::callbackJac(const f64* xu_nlp, const f64* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const f64* xu_ij = xu_nlp + (x_size + u_size) * i;
        f64 x1 = xu_ij[0];
        f64 u = xu_ij[2];
        jac_buffer[4 * i] = -(u + u*u / 2);
        jac_buffer[4 * i + 1] = -(1 + u) * x1;
        jac_buffer[4 * i + 2] = u;
        jac_buffer[4 * i + 3] = x1;
    }
}

void FullSweepTestImpl::callback_hes(const f64* xu_nlp, const f64* p) {
    for (int i = 0; i < mesh->node_count; i++) {
        const f64* xu_ij = xu_nlp + (x_size + u_size) * i;
        f64 x1 = xu_ij[0];
        f64 u = xu_ij[2];
        hes_buffer[3 * i] = -(1 + u);
        hes_buffer[3 * i + 1] = -x1;
        hes_buffer[3 * i + 2] = 1;
    }
}

BoundarySweepTestImpl::BoundarySweepTestImpl(FixedVector<FunctionMR>&& mr_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds)
    : BoundarySweep(std::move(mr_in), mesh, r_bounds, true, 0, 2, 0) {
        mr[0].eval = &eval_buffer[0];
        mr[0].jac.dxf[0].value = &jac_buffer[0];
    }

void BoundarySweepTestImpl::callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    eval_buffer[0] = -xf_nlp[1];
}

void BoundarySweepTestImpl::callbackJac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    jac_buffer[0] = -1;
}

void BoundarySweepTestImpl::callback_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {

};
