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
    off_acc_xu = mesh->createAccOffsetXU(off_x, off_xu);
    off_xu_total = off_acc_xu.back().back() + off_xu;
    number_vars = off_xu_total + problem->p_size;
    number_constraints = problem->boundary.r_size + (problem->full.f_size + problem->full.g_size) * mesh->node_count;
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
                g_lb[mesh->acc_nodes[i][j] + f_index] = 0;
                g_ub[mesh->acc_nodes[i][j] + f_index] = 0;
            }
            for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
                g_lb[mesh->acc_nodes[i][j] + problem->full.f_size + g_index] = problem->full.g_bounds[g_index].lb;
                g_ub[mesh->acc_nodes[i][j] + problem->full.f_size + g_index] = problem->full.g_bounds[g_index].ub;
            }
        }
    }

    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {
        g_lb[mesh->node_count + r_index] = problem->boundary.r_bounds[r_index].lb;
        g_ub[mesh->node_count + r_index] = problem->boundary.r_bounds[r_index].ub;
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

void NLP::eval_f_safe(const double* nlp_solver_x, bool new_x) {
    evaluation_state.check_reset(new_x);
    if (!evaluation_state.x_unscaled) {
        // Scaler.scale(nlp_solver_x, curr_x), perform scaling here, memcpy nlp_solver_x -> unscaled -> scale
        // rn we just memset the data to curr_x
        std::memcpy(curr_x.get(), nlp_solver_x, number_vars * sizeof(double));
        evaluation_state.x_unscaled = true;
    }
    if (evaluation_state.eval_f) {
        return;
    }
    else {
        callback_evaluation();
        eval_f();
    }
}

void NLP::eval_f() {
    double mayer = 0;
    if (problem->boundary.has_mayer) {
        mayer = (*problem->boundary.mr[0].eval);
    };

    double lagrange = 0;
    if (problem->full.has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                lagrange += mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(problem->full.lfg[0].eval + mesh->acc_nodes[i][j]));
            }
        }
    }

    curr_obj = mayer + lagrange;
}

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
