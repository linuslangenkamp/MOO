#include "gdop_problem.h"

// Dummy implementation of FullSweep_OM
FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds, 
                           DATA* data, threadData_t* threadData, InfoGDOP& info)
    : FullSweep(std::move(lfg), mesh, std::move(g_bounds), info.lagrange_exists, info.f_size,
                info.g_size, info.x_size, info.u_size, info.p_size),
      data(data), threadData(threadData), info(info) {
}
void FullSweep_OM::callback_eval(const F64* xu_nlp, const F64* p) {
    set_parameters(data, threadData, info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(data, threadData, info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(data, threadData, info, mesh.t[i][j]);
            eval_current_point(data, threadData, info);
            eval_lfg_write_to_buffer(data, threadData, info, eval_buffer);
        }
    }
}

void FullSweep_OM::callback_jac(const F64* xu_nlp, const F64* p) {
    set_parameters(data, threadData, info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(data, threadData, info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(data, threadData, info, mesh.t[i][j]);
            // eval
            // jacobian
            // write jacobian
        }
    }
}

void FullSweep_OM::callback_hes(const F64* xu_nlp, const F64* p) {
}

// Dummy implementation of BoundarySweep_OM
BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                                   DATA* data, threadData_t* threadData, InfoGDOP& info)
    : BoundarySweep(std::move(mr), mesh, std::move(r_bounds), info.mayer_exists,
                    info.r_size, info.x_size, info.p_size),
      data(data), threadData(threadData), info(info) {
}

void BoundarySweep_OM::callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
    set_parameters(data, threadData, info, p);
    set_states(data, threadData, info, xf_nlp);
    set_time(data, threadData, info, mesh.tf);
    eval_current_point(data, threadData, info);
    eval_mr_write_to_buffer(data, threadData, info, eval_buffer);
}

void BoundarySweep_OM::callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
}

void BoundarySweep_OM::callback_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p) {
}

Problem create_gdop(DATA* data, threadData_t* threadData, InfoGDOP& info, Mesh& mesh) {
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
    info.u_indices_real_vars = FixedVector<int>(info.u_size);
    data->callback->getInputVarIndices(data, info.u_indices_real_vars.raw());
    for (int u = 0; u < info.u_size; u++) {
        int u_index = info.u_indices_real_vars[u];
        u_bounds[u].lb = data->modelData->realVarsData[u_index].attribute.min;
        u_bounds[u].ub = data->modelData->realVarsData[u_index].attribute.max;
    }

    /* constraint sizes */
    info.f_size = info.x_size;
    info.index_der_x_real_vars = info.x_size;
    info.g_size = data->modelData->nOptimizeConstraints;
    info.r_size = data->modelData->nOptimizeFinalConstraints; // TODO: add *generic boundary* constraints later also at t=t0

    // TODO: figure this out; we get some derivative ptrs i guess?
    short der_index_mayer_realVars = -1;
    short der_indices_lagrange_realVars[2] = {-1, -1};

    // this is really ugly IMO, fix this when ready for master!
    info.mayer_exists = (data->callback->mayer(data, &info.__address_mayer_real_vars, &der_index_mayer_realVars) >= 0);
    if (info.mayer_exists) {
        info.index_mayer_real_vars = (int)(info.__address_mayer_real_vars - data->localData[0]->realVars);
    }

    info.lagrange_exists = (data->callback->lagrange(data, &info.__address_lagrange_real_vars, &der_indices_lagrange_realVars[0], &der_indices_lagrange_realVars[1]) >= 0);
    if (info.lagrange_exists) {
        info.index_lagrange_real_vars = (int)(info.__address_lagrange_real_vars - data->localData[0]->realVars);
    }

    /* constraint bounds */
    FixedVector<Bounds> g_bounds(info.g_size);
    FixedVector<Bounds> r_bounds(info.r_size);

    info.index_g_real_vars = data->modelData->nVariablesReal - (info.g_size + info.r_size);
    info.index_r_real_vars = data->modelData->nVariablesReal - info.r_size;
    for (int g = 0; g < info.g_size; g++) {
        g_bounds[g].lb = data->modelData->realVarsData[info.index_g_real_vars + g].attribute.min;
        g_bounds[g].ub = data->modelData->realVarsData[info.index_g_real_vars + g].attribute.max;
    }

    for (int r = 0; r < info.r_size; r++) {
        r_bounds[r].lb = data->modelData->realVarsData[info.index_r_real_vars + r].attribute.min;
        r_bounds[r].ub = data->modelData->realVarsData[info.index_r_real_vars + r].attribute.max;
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

    /* create functions and bounds */
    FixedVector<FunctionMR> mr((int)info.mayer_exists + info.r_size);
    FixedVector<FunctionLFG> lfg((int)info.lagrange_exists + info.f_size + info.g_size);

    /* create CSC <-> COO exchange, init jacobians */
    ExchangeJacobians exc_jac(data, threadData, info);
    
    /* init (OPT) */
    init_eval(data, threadData, info, lfg, mr);
    init_jac(data, threadData, info, exc_jac, lfg, mr);

    FullSweep_OM full(std::move(lfg), mesh, std::move(g_bounds), data, threadData, info);
    BoundarySweep_OM boundary(std::move(mr), mesh, std::move(r_bounds), data, threadData, info);
    
    Problem problem = Problem{
        full,
        boundary,
        mesh,
        std::move(x_bounds),
        std::move(u_bounds),
        std::move(p_bounds),
        std::move(x0_fixed),
        std::move(xf_fixed)
    };

    return problem; 
}
