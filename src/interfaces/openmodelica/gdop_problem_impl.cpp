#include "gdop_problem_impl.h"

// Dummy implementation of FullSweep_OM
FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds, 
                           bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size)
    : FullSweep(std::move(lfg), mesh, g_bounds, has_lagrange, f_size, g_size, x_size, u_size, p_size) {
    // Dummy constructor (does nothing)
}
void FullSweep_OM::callback_eval(const double* xu_nlp, const double* p) {
    // Dummy evaluation (does nothing)
}

void FullSweep_OM::callback_jac(const double* xu_nlp, const double* p) {
    // Dummy Jacobian (does nothing)
}

void FullSweep_OM::callback_hes(const double* xu_nlp, const double* p) {
    // Dummy Hessian (does nothing)
}

// Dummy implementation of BoundarySweep_OM
BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds,
                                   bool has_mayer, int r_size, int x_size, int p_size)
    : BoundarySweep(std::move(mr), mesh, r_bounds, has_mayer, r_size, x_size, p_size) {
    // Dummy constructor (does nothing)
}

void BoundarySweep_OM::callback_eval(const double* x0_nlp, const double* xf_nlp, const double* p) {
    // Dummy evaluation (does nothing)
}

void BoundarySweep_OM::callback_jac(const double* x0_nlp, const double* xf_nlp, const double* p) {
    // Dummy Jacobian (does nothing)
}

void BoundarySweep_OM::callback_hes(const double* x0_nlp, const double* xf_nlp, const double* p) {
    // Dummy Hessian (does nothing)
}

std::shared_ptr<Problem> create_gdop_om(DATA* data, std::shared_ptr<Mesh> mesh) {
    // sizes
    int size_x = data->modelData->nStates;
    int size_u = data->modelData->nInputVars;
    int size_p = 0; // TODO: add this feature

    // TODO: figure this out we get some derivative? ptrs i guess?
    short der_index_m = -1;
    short der_indices_l[2] = {-1, -1};
    double* address_m;
    double* address_l;
    // this is really ugly IMO, fix this when ready for master!
    bool mayer = (data->callback->mayer(data, &address_m, &der_index_m) >= 0);
    bool lagrange = (data->callback->lagrange(data, &address_l, &der_indices_l[0], &der_indices_l[1]) >= 0);

    int size_g = data->modelData->nOptimizeConstraints;
    int size_r = data->modelData->nOptimizeFinalConstraints; // TODO: add *generic boundary* constraints later also at t=t0

    // create functions and bounds
    FixedVector<FunctionMR> mr((int)mayer + size_r);
    FixedVector<FunctionLFG> lfg((int)lagrange + size_x + size_g); // size_f == size_x

    /* create sparsity patterns */

    FixedVector<Bounds> x_bounds(size_x);
    FixedVector<Bounds> u_bounds(size_u);
    FixedVector<Bounds> p_bounds(size_p);

    FixedVector<Bounds> g_bounds(size_g);
    FixedVector<Bounds> r_bounds(size_r);

    FixedVector<std::optional<double>> x0_fixed(size_x);
    FixedVector<std::optional<double>> xf_fixed(size_x);
    /* set bounds and initial values */

    /*
    std::unique_ptr<FullSweep> fs(new FullSweep_OM(std::move(lfg), mesh, g_bounds));
    std::unique_ptr<BoundarySweep> bs(new BoundarySweep_OM(std::move(mr), mesh, r_bounds));
    std::make_shared<Problem>(std::move(fs), std::move(bs), 
                                     std::move(x_bounds), std::move(u_bounds), std::move(p_bounds),
                                     std::move(x0_fixed), std::move(xf_fixed))
    */

    return NULL;
}


