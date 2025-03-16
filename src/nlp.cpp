#include "nlp.h"

void NLP::init() {
    // init all the structures for optimization
    // n, m, nnz_jac, nnz_hess
    // (create scaling), create NLP bounds
    // sparsity pattern
    initSizesOffsets();
    initBounds();
    initStartingPoint();
}

void NLP::initSizesOffsets() {
    off_x = problem->x_size;
    off_u = problem->u_size;
    off_p = problem->p_size;
    off_xu = off_x + off_u;
    off_acc_xu = mesh->createAccOffsetXU(off_x, off_xu);          // variables  x_ij offset
    off_last_xu = off_acc_xu.back().back();                       // variables final grid point x_ij
    off_xu_total = off_last_xu + off_xu;                          // first parameter
    number_vars = off_xu_total + problem->p_size;
    off_acc_fg = mesh->createAccOffsetFG(problem->full.fg_size);  // constraint f_ij offset
    off_fg_total = mesh->node_count * problem->full.fg_size;      // constraint r_0 offset
    number_constraints = problem->boundary.r_size + off_fg_total;
}

void NLP::initBounds() {
    // TODO: perform scaling before the bounds are set!

    curr_x    = std::make_unique<double[]>(number_vars);
    curr_grad = std::make_unique<double[]>(number_vars);
    x_lb      = std::make_unique<double[]>(number_vars);
    x_ub      = std::make_unique<double[]>(number_vars);
    curr_g    = std::make_unique<double[]>(number_constraints);
    g_lb      = std::make_unique<double[]>(number_constraints);
    g_ub      = std::make_unique<double[]>(number_constraints);

    // standard bounds, but checking for x0_fixed or xf_fixed
    for (int x_index = 0; x_index < off_x; x_index++) {
        x_lb[x_index] = problem->x0_fixed[x_index] ? *problem->x0_fixed[x_index] : problem->x_bounds[x_index].lb;
        x_ub[x_index] = problem->x0_fixed[x_index] ? *problem->x0_fixed[x_index] : problem->x_bounds[x_index].ub;
    }

    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (i == mesh->intervals - 1 && j == mesh->nodes[i] - 1) {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_lb[off_acc_xu[i][j] + x_index] = problem->xf_fixed[x_index] ? *problem->xf_fixed[x_index] : problem->x_bounds[x_index].lb;
                    x_ub[off_acc_xu[i][j] + x_index] = problem->xf_fixed[x_index] ? *problem->xf_fixed[x_index] : problem->x_bounds[x_index].ub;
                }
            } 
            else {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_lb[off_acc_xu[i][j] + x_index] = problem->x_bounds[x_index].lb;
                    x_ub[off_acc_xu[i][j] + x_index] = problem->x_bounds[x_index].ub;
                }
            }
            for (int u_index = 0; u_index < off_u; u_index++) {
                x_lb[off_acc_xu[i][j] + off_x + u_index] = problem->u_bounds[u_index].lb;
                x_ub[off_acc_xu[i][j] + off_x + u_index] = problem->u_bounds[u_index].ub;
            }
        }
    }

    for (int p_index = 0; p_index < off_p; p_index++) {
        x_lb[off_xu_total + p_index] = problem->p_bounds[p_index].lb;
        x_ub[off_xu_total + p_index] = problem->p_bounds[p_index].ub;
    }

    // standard constraint bounds
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full.f_size; f_index++) {
                g_lb[off_acc_fg[i][j] + f_index] = 0;
                g_ub[off_acc_fg[i][j] + f_index] = 0;
            }
            for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
                g_lb[off_acc_fg[i][j] + problem->full.f_size + g_index] = problem->full.g_bounds[g_index].lb;
                g_ub[off_acc_fg[i][j] + problem->full.f_size + g_index] = problem->full.g_bounds[g_index].ub;
            }
        }
    }

    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {
        g_lb[off_fg_total + r_index] = problem->boundary.r_bounds[r_index].lb;
        g_ub[off_fg_total + r_index] = problem->boundary.r_bounds[r_index].ub;
    }
}

void NLP::initStartingPoint() {
    Trajectory new_guess = guess->interpolate(*mesh, *collocation);

    for (int x_index = 0; x_index < off_x; x_index++) {
        curr_x[x_index] = new_guess.x[x_index][0];
    }

    int index = 1;
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int x_index = 0; x_index < off_x; x_index++) {
                curr_x[off_acc_xu[i][j] + x_index] = new_guess.x[x_index][index];
            }
            for (int u_index = 0; u_index < off_u; u_index++) {
                curr_x[off_acc_xu[i][j] + off_x + u_index] = new_guess.u[u_index][index];
            }
            index++;
        }
    }
    for (int p_index = 0; p_index < off_p; p_index++) {
        curr_x[off_xu_total + p_index] = new_guess.p[p_index];
    }
}

/* nlp function evaluations happen in two stages:
 * 1. fill the input buffers -> evaluate all continuous callback functions and fill the output buffers
 * 2. evaluate the nlp function by accessing the structures and double* defined in problem, simply use the filled buffers
 * 
 * since all functions are evaluated in step 2, this can always be executed in parallel even if the callbacks were not parallel!
 * because of this structure, before every nlp evaluation check_new_x has to be performed and step 1 has to be executed in case
 */

void NLP::check_new_x(const double* nlp_solver_x, bool new_x) {
    evaluation_state.check_reset(new_x);
    if (!evaluation_state.x_unscaled) {
        // Scaler.scale(nlp_solver_x, curr_x), perform scaling here, memcpy nlp_solver_x -> unscaled -> scale
        // rn we just memset the data to curr_x
        std::memcpy(curr_x.get(), nlp_solver_x, number_vars * sizeof(double));
        evaluation_state.x_unscaled = true;
    }
}

void NLP::eval_f_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (evaluation_state.eval_f) {
        return;
    }
    else {
        callback_evaluation();
        eval_f();
    }
}

void NLP::eval_g_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (evaluation_state.eval_g) {
        return;
    }
    else {
        callback_evaluation();
        eval_g();
    }
}


void NLP::eval_grad_f_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (evaluation_state.grad_f) {
        return;
    }
    else {
        callback_jacobian();
        eval_grad_f();
    }
};

void NLP::eval_f() {
    double mayer = 0;
    if (problem->boundary.has_mayer) {
        mayer = problem->boundary.getEvalM();
    };

    double lagrange = 0;
    if (problem->full.has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                lagrange += mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * problem->full.getEvalL(mesh->acc_nodes[i][j]);
            }
        }
    }

    curr_obj = mayer + lagrange;
}


void NLP::eval_grad_f() {
    std::memset(curr_grad.get(), 0, number_vars * sizeof(double));

    if (problem->full.has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                for (auto& dL_dx : problem->full.lfg[0].jac.dx) {
                    curr_grad[off_acc_xu[i][j] + dL_dx.index] = mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*dL_dx.value);
                }
                for (auto& dL_du : problem->full.lfg[0].jac.du) {
                    curr_grad[off_acc_xu[i][j] + off_x + dL_du.index] = mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*dL_du.value);
                }
                for (auto& dL_dp : problem->full.lfg[0].jac.dp) {
                    curr_grad[off_xu_total + dL_dp.index] += mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*dL_dp.value);
                }
            }
        }
    }
    if (problem->boundary.has_mayer) {
        for (auto& dM_dx0 : problem->boundary.mr[0].jac.dx0) {
            curr_grad[dM_dx0.index] = (*dM_dx0.value);
        }
        for (auto& dM_dxf : problem->boundary.mr[0].jac.dxf) {
            curr_grad[off_last_xu + dM_dxf.index] += (*dM_dxf.value);
        }
        for (auto& dM_dp : problem->boundary.mr[0].jac.dp) {
            curr_grad[off_xu_total + dM_dp.index] += (*dM_dp.value);
        }
    }
};

void NLP::eval_g() {
    std::memset(curr_g.get(), 0, number_constraints * sizeof(double));
    for (int i = 0; i < mesh->intervals; i++) {
        for (int f = problem->full.f_index_start; f < problem->full.f_index_end; f++) {
            collocation->diff_matrix_multiply(mesh->nodes[i], off_x, off_xu, problem->full.fg_size,
                                              &curr_x[off_acc_xu[i][0]],   // x_{i-1, m_{i-1}} == x_{i, 0}, base point states
                                              &curr_x[off_acc_xu[i][1]],   // x_{i, 1}, collocation point states
                                              &curr_g[off_acc_fg[i][0]]);  // constraint start index 
        }
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full.f_size; f_index++) {
                curr_g[off_acc_fg[i][j] + f_index] -= mesh->delta_t[i] * problem->full.getEvalF(f_index, mesh->acc_nodes[i][j]);
            }
            for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
                curr_g[off_acc_fg[i][j] + problem->full.f_size + g_index] += problem->full.getEvalG(g_index, mesh->acc_nodes[i][j]);
            }
        }
    }
    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {
        curr_g[off_fg_total + r_index] = problem->boundary.getEvalR(r_index);
    }
}

// TODO: how to perform the data filling and evaluations
void NLP::callback_evaluation() {
    // TODO: how to interface here
    // problem->full.fillInputData();
    evaluation_state.eval_f = true;
    evaluation_state.eval_g = true;
}

void NLP::callback_jacobian() {
    // TODO: how to interface here
    // problem->full.fillInputData();
    evaluation_state.grad_f = true;
    evaluation_state.jac_g = true;
}

void NLP::callback_hessian() {
    // TODO: how to interface here
    // problem->full.fillInputData();
    evaluation_state.hes_lag = true;
}
