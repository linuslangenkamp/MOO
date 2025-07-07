#include "gdop_problem.h"

// use this field when needing some data object inside OpenModelica
// callbacks / interfaces but no nice void* field exists
static void *_global_reference_data_field = nullptr;

// set the global pointer
void set_global_reference_data(void *reference_data) {
    assert(!_global_reference_data_field);         // assert nullptr
    _global_reference_data_field = reference_data;
}

// get the global pointer
void* get_global_reference_data() {
    assert(_global_reference_data_field);          // assert non nullptr
    return _global_reference_data_field;
}

// clear the global pointer
void clear_global_reference_data() {
    _global_reference_data_field = nullptr;
}

FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, std::unique_ptr<AugmentedHessianLFG> aug_hes, std::unique_ptr<AugmentedParameterHessian> aug_pp_hes,
        Collocation& collocation, Mesh& mesh, FixedVector<Bounds>&& g_bounds, InfoGDOP& info)
    : FullSweep(std::move(lfg), std::move(aug_hes), std::move(aug_pp_hes), collocation, mesh, std::move(g_bounds), info.lagrange_exists, info.f_size,
                info.g_size, info.x_size, info.u_size, info.p_size), info(info) {
}
void FullSweep_OM::callback_eval(const f64* xu_nlp, const f64* p) {
    set_parameters(info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(info, mesh.t[i][j]);
            eval_current_point(info);
            eval_lfg_write(info, &eval_buffer[eval_size * mesh.acc_nodes[i][j]]);
        }
    }
}

void FullSweep_OM::callback_jac(const f64* xu_nlp, const f64* p) {
    set_parameters(info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(info, mesh.t[i][j]);
            eval_current_point(info);
            /* TODO: check if B matrix does hold additional ders */
            jac_eval_write_as_csc(info, info.exc_jac->B, &jac_buffer[jac_size * mesh.acc_nodes[i][j]]);
        }
    }
}

void FullSweep_OM::callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, f64* lambda) {
    set_parameters(info, p);
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            set_states_inputs(info, &xu_nlp[info.xu_size * mesh.acc_nodes[i][j]]);
            set_time(info, mesh.t[i][j]);
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
            /* set previous Jacobian *CSC* OpenModelica buffer */
            info.exc_hes->B_args.jac_csc = &jac_buffer[jac_size * mesh.acc_nodes[i][j]];

            /* call Hessian */
            __richardsonExtrapolation(info.exc_hes->B_extr, __forwardDiffHessianWrapper, &info.exc_hes->B_args,
                                      1e-6, 1, 2, 1, &aug_hes_buffer[aug_hes_size * mesh.acc_nodes[i][j]]);
        }
    }
}

BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, std::unique_ptr<AugmentedHessianMR> aug_hes, Mesh& mesh,
                                   FixedVector<Bounds>&& r_bounds, InfoGDOP& info)
    : BoundarySweep(std::move(mr), std::move(aug_hes), mesh, std::move(r_bounds), info.mayer_exists,
                    info.r_size, info.x_size, info.p_size), info(info) {}

void BoundarySweep_OM::callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    set_parameters(info, p);
    set_states(info, xf_nlp);
    set_time(info, mesh.tf);
    eval_current_point(info);
    eval_mr_write(info, eval_buffer.raw());
}

void BoundarySweep_OM::callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    set_parameters(info, p);
    set_states(info, xf_nlp);
    set_time(info, mesh.tf);
    eval_current_point(info);
    /* TODO: check if C matrix does hold additional ders */
    /* derivative of mayer to jacbuffer[0] ... jac_buffer[exc_jac.D_coo.nnz_offset - 1] */
    if (has_mayer) {
        jac_eval_write_first_row_as_csc(info, info.exc_jac->C, info.exc_jac->C_buffer.raw(),
                                        jac_buffer.raw(), info.exc_jac->C_coo);
    }

    if (info.exc_jac->D_exists) {
        /* TODO: check if D matrix does hold additional ders */
        jac_eval_write_as_csc(info, info.exc_jac->D, &jac_buffer[info.exc_jac->D_coo.nnz_offset]);
    }
}

void BoundarySweep_OM::callback_aug_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p, const f64 mayer_factor, f64* lambda) {
    set_parameters(info, p);
    set_states(info, xf_nlp);
    set_time(info, mesh.tf);
    aug_hes_buffer.fill_zero();

    if (has_mayer) {
        int index_mayer = info.x_size + (int)(info.lagrange_exists);
        info.exc_hes->C_lambda[index_mayer] = mayer_factor;
        __richardsonExtrapolation(info.exc_hes->C_extr, __forwardDiffHessianWrapper, &info.exc_hes->C_args,
                                  1e-6, 1, 2, 1, info.exc_hes->C_buffer.raw());
        for (auto& [index_C, index_buffer] : info.exc_hes->C_to_Mr_buffer) {
            aug_hes_buffer[index_buffer] += info.exc_hes->C_buffer[index_C];
        }
    }

    if (r_size != 0) {
        /* set duals and precomputed Jacobian D */
        info.exc_hes->D_args.lambda = lambda;
        info.exc_hes->D_args.jac_csc = &jac_buffer[info.exc_jac->D_coo.nnz_offset];

        __richardsonExtrapolation(info.exc_hes->D_extr, __forwardDiffHessianWrapper, &info.exc_hes->D_args,
                                  1e-6, 1, 2, 1, info.exc_hes->D_buffer.raw());
        for (auto& [index_D, index_buffer] : info.exc_hes->D_to_Mr_buffer) {
            aug_hes_buffer[index_buffer] += info.exc_hes->D_buffer[index_D];
        }
    }
}

Problem create_gdop(InfoGDOP& info, Mesh& mesh, Collocation& collocation) {
    DATA* data = info.data;
    // at first call init for all start values
    initialize_model(info);
    
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
    info.mayer_exists = (data->callback->mayer(data, &info.address_mayer_real_vars, &der_index_mayer_realVars) >= 0);
    if (info.mayer_exists) {
        info.index_mayer_real_vars = (int)(info.address_mayer_real_vars - data->localData[0]->realVars);
    }

    info.lagrange_exists = (data->callback->lagrange(data, &info.address_lagrange_real_vars, &der_indices_lagrange_realVars[0], &der_indices_lagrange_realVars[1]) >= 0);
    if (info.lagrange_exists) {
        info.index_lagrange_real_vars = (int)(info.address_lagrange_real_vars - data->localData[0]->realVars);
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
    info.exc_jac = std::make_unique<ExchangeJacobians>(info);

    /* create HESSIAN_PATTERNs and allocate buffers for extrapolation / evaluation */
    info.exc_hes = std::make_unique<ExchangeHessians>(info);
    auto aug_hes_lfg = std::make_unique<AugmentedHessianLFG>();
    auto aug_hes_lfg_pp = std::make_unique<AugmentedParameterHessian>();
    auto aug_hes_mr = std::make_unique<AugmentedHessianMR>();

    /* init (OPT) */
    init_eval(info, lfg, mr);
    init_jac(info, lfg, mr);
    init_hes(info, *aug_hes_lfg, *aug_hes_lfg_pp, *aug_hes_mr, mr);

    return Problem(
        std::make_unique<FullSweep_OM>(std::move(lfg), std::move(aug_hes_lfg), std::move(aug_hes_lfg_pp), collocation, mesh, std::move(g_bounds), info),
        std::make_unique<BoundarySweep_OM>(std::move(mr), std::move(aug_hes_mr), mesh, std::move(r_bounds), info),
        mesh,
        std::move(x_bounds),
        std::move(u_bounds),
        std::move(p_bounds),
        std::move(x0_fixed),
        std::move(xf_fixed)
    );
}

// TODO: make a NLP initializer class with different init methods: simulation, bionic, constant, ...
//       always returns a unique_ptr<Trajectory>
std::unique_ptr<Trajectory> create_constant_guess(InfoGDOP& info) {
    DATA* data = info.data;

    std::vector<f64> t = {0, info.tf};
    std::vector<std::vector<f64>> x_guess;
    std::vector<std::vector<f64>> u_guess;
    std::vector<f64> p;
    InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    for (int x = 0; x < info.x_size; x++) {
        x_guess.push_back({data->modelData->realVarsData[x].attribute.start, data->modelData->realVarsData[x].attribute.start});
    }

    for (int u : info.u_indices_real_vars) {
        u_guess.push_back({data->modelData->realVarsData[u].attribute.start, data->modelData->realVarsData[u].attribute.start});
    }
    // TODO: add p

    return std::make_unique<Trajectory>(Trajectory{t, x_guess, u_guess, p, interpolation});
}

static int control_trajectory_input_function(DATA* data, threadData_t* threadData) {
    AuxiliaryControls* aux_controls = (AuxiliaryControls*)get_global_reference_data();

    ControlTrajectory& controls = aux_controls->controls;
    InfoGDOP& info = aux_controls->info;
    f64* u_interpolation = aux_controls->u_interpolation.raw();

    // important use data here not info.data
    // for some reason a new data object is created for each solve
    // but not important for now
    f64 time = data->localData[0]->timeValue;

    // transform back to GDOP time (always [0, tf])
    time -= info.start_time;

    controls.interpolate(time, u_interpolation);
    set_inputs(info, u_interpolation);

    return 0;
}

static void trajectory_xut_emit(simulation_result* sim_result, DATA* data, threadData_t *threadData)
{
    AuxiliaryTrajectory* aux = (AuxiliaryTrajectory*)sim_result->storage; // exploit void* field
    InfoGDOP& info = aux->info;
    Trajectory& trajectory = aux->trajectory;
    SOLVER_INFO* solver_info = aux->solver_info;

    for (int x_idx = 0; x_idx < info.x_size; x_idx++) {
        trajectory.x[x_idx].push_back(data->localData[0]->realVars[x_idx]);
    }

    for (int u_idx = 0; u_idx < info.u_size; u_idx++) {
        int u = info.u_indices_real_vars[u_idx];
        trajectory.u[u_idx].push_back(data->localData[0]->realVars[u]);
    }

    trajectory.t.push_back(solver_info->currentTime - info.start_time);
}

static void trajectory_p_emit(simulation_result* sim_result, DATA* data, threadData_t *threadData)
{
    // TODO: parameters
    AuxiliaryTrajectory* aux = (AuxiliaryTrajectory*)sim_result->storage;
    InfoGDOP& info = aux->info;
    Trajectory& trajectory = aux->trajectory;

    for (int p_idx = 0; p_idx < info.p_size; p_idx++) {
        trajectory.p.push_back(data->localData[0]->realVars[0 /* parameter index */]);
    }
}

// sets state and control initial values
// from e.g. initial equations / parameters
void initialize_model(InfoGDOP& info) {
    externalInputallocate(info.data);
    initializeModel(info.data, info.threadData, "", "", info.start_time);
}

std::unique_ptr<Trajectory> simulate(InfoGDOP& info, SOLVER_METHOD solver,
                                     int num_steps, ControlTrajectory& controls) {
    DATA* data = info.data;
    threadData_t* threadData = info.threadData;
    SOLVER_INFO solver_info;
    SIMULATION_INFO *simInfo = data->simulationInfo;

    simInfo->numSteps = num_steps;
    simInfo->stepSize = info.tf / (f64)num_steps;
    simInfo->useStopTime = 1;
    solver_info.solverMethod = solver;

    // allocate and reserve trajectory vectors
    std::vector<f64> t;
    t.reserve(num_steps + 1);

    std::vector<std::vector<f64>> x_sim(info.x_size);
    for (auto& v : x_sim) v.reserve(num_steps + 1);

    std::vector<std::vector<f64>> u_sim(info.u_size);
    for (auto& v : u_sim) v.reserve(num_steps + 1);

    std::vector<f64> p_sim(info.p_size);

    // create Trajectory object
    auto trajectory = std::make_unique<Trajectory>(Trajectory{t, x_sim, u_sim, p_sim, InterpolationMethod::LINEAR});

    // auxiliary data (passed as void* in storage member of sim_result)
    auto aux_trajectory = std::make_unique<AuxiliaryTrajectory>(AuxiliaryTrajectory{*trajectory, info, &solver_info});

    // define global sim_result
    sim_result.filename = NULL;
    sim_result.numpoints = 0;
    sim_result.cpuTime = 0;
    sim_result.storage = aux_trajectory.get();
    sim_result.emit = trajectory_xut_emit;
    sim_result.init = nullptr;
    sim_result.writeParameterData = trajectory_p_emit;
    sim_result.free = nullptr;

    // init simulation
    initializeSolverData(data, threadData, &solver_info);
    setZCtol(fmin(simInfo->stepSize, simInfo->tolerance));
    initialize_model(info);
    data->real_time_sync.enabled = FALSE;

    // create an auxiliary object (stored in global void*)
    // since the input_function interface offers no additional argument
    FixedVector<f64> u_interpolation_buffer = FixedVector<f64>(info.u_size);
    auto aux_controls = AuxiliaryControls{controls, info, u_interpolation_buffer};
    set_global_reference_data(&aux_controls);

    // set the new input function
    auto generated_input_function = data->callback->input_function;
    data->callback->input_function = control_trajectory_input_function;

    // set controls for start time
    // emit for time = t_0
    controls.interpolate(0.0, u_interpolation_buffer.raw());
    set_inputs(info, u_interpolation_buffer.raw());
    trajectory_xut_emit(&sim_result, data, threadData);

    // simulation with custom emit
    data->callback->performSimulation(data, threadData, &solver_info);

    // reset to previous input function (from generated code)
    data->callback->input_function = generated_input_function;

    // set global aux data to nullptr
    clear_global_reference_data();

    return trajectory;
}

void emit_trajectory_om(Trajectory& trajectory, InfoGDOP& info) {
    DATA* data = info.data;
    threadData_t* threadData = info.threadData;

    // setting default emitter
    if (!data->modelData->resultFileName) {
        std::string result_file = std::string(data->modelData->modelFilePrefix) + "_res." + data->simulationInfo->outputFormat;
        data->modelData->resultFileName = GC_strdup(result_file.c_str());
    }
    data->simulationInfo->numSteps = trajectory.t.size();
    initializeResultData(data, threadData, 0);
    sim_result.writeParameterData(&sim_result, data, threadData);

    // allocate contiguous array for xu
    FixedVector<f64> xu(trajectory.x.size() + trajectory.u.size());
    for (size_t i = 0; i < trajectory.t.size(); i++) {
        // move trajectory data in contiguous array
        for (size_t x_index = 0; x_index < trajectory.x.size(); x_index++) {
            xu[x_index] = trajectory.x[x_index][i];
        }
        for (size_t u_index = 0; u_index < trajectory.u.size(); u_index++) {
            xu[trajectory.x.size() + u_index] = trajectory.u[u_index][i];
        }

        // evaluate all algebraic variables
        set_time(info, info.start_time + trajectory.t[i]);
        set_states_inputs(info, xu.raw());
        eval_current_point(info);

        // emit point
        sim_result.emit(&sim_result, data, threadData);
    }

    sim_result.free(&sim_result, data, threadData);
}

NominalScaling create_gdop_om_nominal_scaling(GDOP& gdop, InfoGDOP& info) {
    // x, g, f of the NLP { min f(x) s.t. g_l <= g(x) <= g_l }
    auto x_nominal = FixedVector<f64>(gdop.number_vars);
    auto g_nominal = FixedVector<f64>(gdop.number_constraints);
    f64  f_nominal = 1;

    // get problem sizes
    auto x_size  = info.x_size;
    auto u_size  = info.u_size;
    auto xu_size = info.xu_size;
    auto f_size = info.f_size;
    auto g_size = info.g_size;
    auto r_size = info.r_size;
    auto fg_size = f_size + g_size;

    auto real_vars_data = info.data->modelData->realVarsData;

    auto has_mayer = gdop.problem.boundary->has_mayer;
    auto has_lagrange = gdop.problem.full->has_lagrange;

    if (has_mayer && has_lagrange) {
        f_nominal = (real_vars_data[info.index_mayer_real_vars].attribute.nominal +
                     real_vars_data[info.index_lagrange_real_vars].attribute.nominal) / 2;
    }
    else if (has_lagrange) {
        f_nominal = real_vars_data[info.index_lagrange_real_vars].attribute.nominal;
    }
    else if (has_mayer) {
        f_nominal = real_vars_data[info.index_mayer_real_vars].attribute.nominal;
    }
    else {
        // use default 1 - no objective set!
    }

    // x(t_0)
    for (int x = 0; x < info.x_size; x++) {
        x_nominal[x] = real_vars_data[x].attribute.nominal;
    }

    // (x, u)_(t_node)
    for (int node = 0; node < gdop.mesh.node_count; node++) {
        for (int x = 0; x < x_size; x++) {
            x_nominal[x_size + node * xu_size + x] = real_vars_data[x].attribute.nominal;
        }

        for (int u = 0; u < u_size; u++) {
            int u_real_vars = info.u_indices_real_vars[u];
            x_nominal[2 * x_size + node * xu_size + u] = real_vars_data[u_real_vars].attribute.nominal;
        }
    }

    for (int node = 0; node < gdop.mesh.node_count; node++) {
        for (int f = 0; f < f_size; f++) {
            g_nominal[node * fg_size + f] = x_nominal[f]; // reuse x nominal for dynamic for now!
        }

        for (int g = 0; g < g_size; g++) {
            g_nominal[f_size + node * fg_size + g] = real_vars_data[info.index_g_real_vars + g].attribute.nominal;
        }
    }

    for (int r = 0; r < r_size; r++) {
        g_nominal[gdop.off_fg_total + r] = real_vars_data[info.index_r_real_vars + r].attribute.nominal;
    }

    return NominalScaling(std::move(x_nominal), std::move(g_nominal), f_nominal);
}
