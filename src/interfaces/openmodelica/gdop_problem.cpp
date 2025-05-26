#include "gdop_problem.h"

// Dummy implementation of FullSweep_OM
FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, std::unique_ptr<AugmentedHessianLFG> aug_hes, std::unique_ptr<AugmentedParameterHessian> aug_pp_hes,
        Collocation& collocation, Mesh& mesh, FixedVector<Bounds>&& g_bounds, DATA* data, threadData_t* threadData, InfoGDOP& info)
    : FullSweep(std::move(lfg), std::move(aug_hes), std::move(aug_pp_hes), collocation, mesh, std::move(g_bounds), info.lagrange_exists, info.f_size,
                info.g_size, info.x_size, info.u_size, info.p_size),
      data(data), threadData(threadData), info(info) {
}
void FullSweep_OM::callback_eval(const f64* xu_nlp, const f64* p) {
    set_parameters(data, threadData, info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(data, threadData, info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(data, threadData, info, mesh.t[i][j]);
            eval_current_point(data, threadData, info);
            eval_lfg_write(data, threadData, info, &eval_buffer[eval_size * mesh.acc_nodes[i][j]]);
        }
    }
}

void FullSweep_OM::callback_jac(const f64* xu_nlp, const f64* p) {
    set_parameters(data, threadData, info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(data, threadData, info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(data, threadData, info, mesh.t[i][j]);
            eval_current_point(data, threadData, info);
            /* TODO: check if B matrix does hold additional ders */
            jac_eval_write_as_csc(data, threadData, info, info.exc_jac->B, &jac_buffer[jac_size * mesh.acc_nodes[i][j]]);
        }
    }
}

void FullSweep_OM::callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, f64* lambda) {
    set_parameters(data, threadData, info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(data, threadData, info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(data, threadData, info, mesh.t[i][j]);
            eval_current_point(data, threadData, info);
            /* TODO: check if B matrix does hold additional ders */
            if (has_lagrange) {
                /* OpenModelica sorts the Functions as fLg, we have to swap the order for lambda
                 * thus, use the workspace buffer from info.exc_hes */
                info.exc_hes->B_lambda[x_size] = lagrange_factors[i][j];
                for (int f = 0; f < f_size; f++) {
                    info.exc_hes->B_lambda[f] = lambda[fg_size * mesh.acc_nodes[i][j] + f];
                }
                for (int g = 0; g < g_size; g++) {
                    /* Lagrange offset */
                    info.exc_hes->B_lambda[f_size + 1 + g] = lambda[fg_size * mesh.acc_nodes[i][j] + f_size + g];
                }
                /* set wrapper lambda */
                info.exc_hes->B_args.lambda = info.exc_hes->B_lambda.raw();
            }
            else {
                /* set wrapper lambda */
                info.exc_hes->B_args.lambda = &lambda[fg_size * mesh.acc_nodes[i][j]];
            }
            // FIXME: TODO: Make this use the previous Jacobian
            __richardsonExtrapolation(info.exc_hes->B_extr, __forwardDiffHessianWrapper, &info.exc_hes->B_args,
                                      1e-8, 1, 2, 1, &aug_hes_buffer[aug_hes_size * mesh.acc_nodes[i][j]]);
        }
    }
}

BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, std::unique_ptr<AugmentedHessianMR> aug_hes, Mesh& mesh,
                                   FixedVector<Bounds>&& r_bounds, DATA* data, threadData_t* threadData, InfoGDOP& info)
    : BoundarySweep(std::move(mr), std::move(aug_hes), mesh, std::move(r_bounds), info.mayer_exists,
                    info.r_size, info.x_size, info.p_size),
      data(data), threadData(threadData), info(info) {
}

void BoundarySweep_OM::callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    set_parameters(data, threadData, info, p);
    set_states(data, threadData, info, xf_nlp);
    set_time(data, threadData, info, mesh.tf);
    eval_current_point(data, threadData, info);
    eval_mr_write(data, threadData, info, eval_buffer.raw());
}

void BoundarySweep_OM::callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    set_parameters(data, threadData, info, p);
    set_states(data, threadData, info, xf_nlp);
    set_time(data, threadData, info, mesh.tf);
    eval_current_point(data, threadData, info);
    /* TODO: check if C matrix does hold additional ders */
    /* derivative of mayer to jacbuffer[0] ... jac_buffer[exc_jac.D_coo.nnz_offset - 1] */
    if (has_mayer) {
        jac_eval_write_first_row_as_csc(data, threadData, info, info.exc_jac->C, info.exc_jac->C_buffer.raw(),
                                        jac_buffer.raw(), info.exc_jac->C_coo);
    }

    if (info.exc_jac->D_exists) {
        /* TODO: check if D matrix does hold additional ders */
        jac_eval_write_as_csc(data, threadData, info, info.exc_jac->D, &jac_buffer[info.exc_jac->D_coo.nnz_offset]);
    }
}

void BoundarySweep_OM::callback_aug_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p, const f64 mayer_factor, f64* lambda) {
    set_parameters(data, threadData, info, p);
    set_states(data, threadData, info, xf_nlp);
    set_time(data, threadData, info, mesh.tf);
    eval_current_point(data, threadData, info);
    aug_hes_buffer.fill_zero();

    if (has_mayer) {
        int index_mayer = info.x_size + (int)(info.lagrange_exists);
        info.exc_hes->C_lambda[index_mayer] = mayer_factor;
        __richardsonExtrapolation(info.exc_hes->C_extr, __forwardDiffHessianWrapper, &info.exc_hes->C_args,
                                  1e-8, 1, 2, 1, info.exc_hes->C_buffer.raw());
        for (auto& [index_C, index_buffer] : info.exc_hes->C_to_Mr_buffer) {
            aug_hes_buffer[index_buffer] += info.exc_hes->C_buffer[index_C];
        }
    }

    if (r_size != 0) {
        info.exc_hes->D_args.lambda = lambda;
        __richardsonExtrapolation(info.exc_hes->D_extr, __forwardDiffHessianWrapper, &info.exc_hes->D_args,
                                  1e-8, 1, 2, 1, info.exc_hes->D_buffer.raw());
        for (auto& [index_D, index_buffer] : info.exc_hes->D_to_Mr_buffer) {
            aug_hes_buffer[index_buffer] += info.exc_hes->D_buffer[index_D];
        }
    }
}

Problem create_gdop(DATA* data, threadData_t* threadData, InfoGDOP& info, Mesh& mesh, Collocation& collocation) {
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
    data->callback->getInputVarIndicesInOptimization(data, info.u_indices_real_vars.raw()); // TODO: use these everywhere, maybe we need to add a buffer in OM C-SimRuntime for it
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

    /* this is really ugly IMO, fix this when ready for master! */
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
     * option: generate fixed final states individually */
    FixedVector<std::optional<f64>> x0_fixed(info.x_size);
    FixedVector<std::optional<f64>> xf_fixed(info.x_size);

    /* set *fixed* initial, final states */
    for (int x = 0; x < info.x_size; x++) {
        x0_fixed[x] = data->modelData->realVarsData[x].attribute.start;
    }

    /* create functions and bounds */
    FixedVector<FunctionMR> mr((int)info.mayer_exists + info.r_size);
    FixedVector<FunctionLFG> lfg((int)info.lagrange_exists + info.f_size + info.g_size);

    /* create CSC <-> COO exchange, init jacobians */
    info.exc_jac = std::make_unique<ExchangeJacobians>(data, threadData, info);

    /* create HESSIAN_PATTERNs and allocate buffers for extrapolation / evaluation */
    info.exc_hes = std::make_unique<ExchangeHessians>(data, threadData, info);
    auto aug_hes_lfg = std::make_unique<AugmentedHessianLFG>();
    auto aug_hes_lfg_pp = std::make_unique<AugmentedParameterHessian>();
    auto aug_hes_mr = std::make_unique<AugmentedHessianMR>();

    /* init (OPT) */
    init_eval(data, threadData, info, lfg, mr);
    init_jac(data, threadData, info, lfg, mr);
    init_hes(data, threadData, info, *aug_hes_lfg, *aug_hes_lfg_pp, *aug_hes_mr, mr);
    
    return Problem(
        std::make_unique<FullSweep_OM>(std::move(lfg), std::move(aug_hes_lfg), std::move(aug_hes_lfg_pp), collocation, mesh, std::move(g_bounds), data, threadData, info),
        std::make_unique<BoundarySweep_OM>(std::move(mr), std::move(aug_hes_mr), mesh, std::move(r_bounds), data, threadData, info),
        mesh,
        std::move(x_bounds),
        std::move(u_bounds),
        std::move(p_bounds),
        std::move(x0_fixed),
        std::move(xf_fixed)
    );
}

// TODO: make a NLP initializer class with different init methods: simulation, bionic, constant, ...
Trajectory create_constant_guess(DATA* data, threadData_t* threadData, InfoGDOP& info) {
    std::vector<f64> t = {0, info.tf};
    std::vector<std::vector<f64>> x_guess;
    std::vector<std::vector<f64>> u_guess;
    std::vector<f64> p;
    InterpolationMethod interpolation = InterpolationMethod::LINEAR;{};

    for (int x = 0; x < info.x_size; x++) {
        x_guess.push_back({data->modelData->realVarsData[x].attribute.start, data->modelData->realVarsData[x].attribute.start});
    }

    for (int u : info.u_indices_real_vars) {
        u_guess.push_back({data->modelData->realVarsData[u].attribute.start, data->modelData->realVarsData[u].attribute.start});
    }
    // TODO: add p

    return Trajectory{t, x_guess, u_guess, p, interpolation};
}
