#include "gdop_problem.h"

// Dummy implementation of FullSweep_OM
FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds, 
                           bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size)
    : FullSweep(std::move(lfg), mesh, std::move(g_bounds), has_lagrange, f_size, g_size, x_size, u_size, p_size) {
    // Dummy constructor (does nothing)
}
void FullSweep_OM::callback_eval(const F64* xu_nlp, const F64* p) {
    // Dummy evaluation (does nothing)
}

void FullSweep_OM::callback_jac(const F64* xu_nlp, const F64* p) {
    // Dummy Jacobian (does nothing)
}

void FullSweep_OM::callback_hes(const F64* xu_nlp, const F64* p) {
    // Dummy Hessian (does nothing)
}

// Dummy implementation of BoundarySweep_OM
BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                                   bool has_mayer, int r_size, int x_size, int p_size)
    : BoundarySweep(std::move(mr), mesh, std::move(r_bounds), has_mayer, r_size, x_size, p_size) {
    // Dummy constructor (does nothing)
}

void BoundarySweep_OM::callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    // Dummy evaluation (does nothing)
}

void BoundarySweep_OM::callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    // Dummy Jacobian (does nothing)
}

void BoundarySweep_OM::callback_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    // Dummy Hessian (does nothing)
}

Problem* create_gdop(DATA* data, threadData_t* threadData, Mesh& mesh) {
    /* create small info struct <-> same purpose as DATA* in OpenModeica */
    InfoGDOP info;

    /* variable sizes */
    info.x_size = data->modelData->nStates;
    info.u_size = data->modelData->nInputVars;
    info.xu_size = info.x_size + info.u_size;
    info.p_size = 0; // TODO: add this feature

    /* variable bounds */
    FixedVector<Bounds> x_bounds(info.x_size);
    FixedVector<Bounds> u_bounds(info.u_size);
    FixedVector<Bounds> p_bounds(info.p_size);

    for (int x = 0; x < info.x_size; x++) {
        x_bounds[x].lb = data->modelData->realVarsData[x].attribute.min;
        x_bounds[x].ub = data->modelData->realVarsData[x].attribute.max;
    }

    /* new generated function getInputVarIndices, just fills the index list of all optimizable inputs */
    auto u_indices_real_vars = std::make_unique<int>(info.u_size);
    data->callback->getInputVarIndices(data, u_indices_real_vars.get());
    for (int u = 0; u < info.u_size; u++) {
        int u_index = *(u_indices_real_vars.get() + u);
        u_bounds[u].lb = data->modelData->realVarsData[u_index].attribute.min;
        u_bounds[u].ub = data->modelData->realVarsData[u_index].attribute.max;
    }

    /* constraint sizes */
    info.f_size = info.x_size;
    info.g_size = data->modelData->nOptimizeConstraints;
    info.r_size = data->modelData->nOptimizeFinalConstraints; // TODO: add *generic boundary* constraints later also at t=t0

    /* constraint bounds */
    FixedVector<Bounds> g_bounds(info.g_size);
    FixedVector<Bounds> r_bounds(info.r_size);

    info.real_vars_g_start_index = data->modelData->nVariablesReal - (info.g_size + info.r_size);
    info.real_vars_r_start_index = data->modelData->nVariablesReal - info.r_size;
    for (int g = 0; g < info.g_size; g++) {
        g_bounds[g].lb = data->modelData->realVarsData[info.real_vars_g_start_index + g].attribute.min;
        g_bounds[g].ub = data->modelData->realVarsData[info.real_vars_g_start_index + g].attribute.max;
    }

    for (int r = 0; r < info.r_size; r++) {
        r_bounds[r].lb = data->modelData->realVarsData[info.real_vars_r_start_index + r].attribute.min;
        r_bounds[r].ub = data->modelData->realVarsData[info.real_vars_r_start_index + r].attribute.max;
    }

    /* for now we ignore xf fixed (need some steps in Backend to detect)
     * and also ignore x0 non fixed, since too complicated
     * => assume x(t_0) = x0 fixed, x(t_f) free to r constraint / maybe the old BE can do that already?!
     * option: generate fixed final states individually
     */
    FixedVector<std::optional<F64>> x0_fixed(info.x_size);
    FixedVector<std::optional<F64>> xf_fixed(info.x_size);

    /* set *fixed* initial, final states */
    for (int x = 0; x < info.x_size; x++) {
        x0_fixed[x] = data->modelData->realVarsData[x].attribute.start;
    }

    // TODO: figure this out; we get some derivative ptrs i guess?
    short der_index_mayer_realVars = -1;
    short der_indices_lagrange_realVars[2] = {-1, -1};

    // this is really ugly IMO, fix this when ready for master!
    info.mayer_exists = (data->callback->mayer(data, &info.address_mayer_real_vars, &der_index_mayer_realVars) >= 0);
    info.lagrange_exists = (data->callback->lagrange(data, &info.address_lagrange_real_vars, &der_indices_lagrange_realVars[0], &der_indices_lagrange_realVars[1]) >= 0);

    /* create functions and bounds */
    FixedVector<FunctionMR> mr((int)info.mayer_exists + info.r_size);
    FixedVector<FunctionLFG> lfg((int)info.lagrange_exists + info.f_size + info.g_size);

    /* create CSC <-> COO exchange, init jacobians */
    ExchangeJacobians exc_jac(data, threadData, info);
    
    /* init (OPT) */
    exc_jac.init_jac(data, threadData, info, lfg, mr);

    // set lfg jac ptrs

    print_jacobian_sparsity(exc_jac.C, true, "C");

    //exc_jac.C_mayer_coo->col.print();
    //exc_jac.C_mayer_coo->csc_to_coo.print();
    /*
    print_jacobian_sparsity(exc_jac.A, true, "A");
    print_jacobian_sparsity(exc_jac.D, true, "D");
    */

    return NULL;
}
