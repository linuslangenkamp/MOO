#include "nlp.h"

void NLP::init() {
    // init all the structures for optimization
    // n, m, nnz_jac, nnz_hess
    // (create scaling), create NLP bounds
    // sparsity pattern
    initSizesOffsets();
    initBounds();
    initInitialPoint();
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

    x_bounds_l = std::make_unique<double[]>(number_vars);
    x_bounds_u = std::make_unique<double[]>(number_vars);
    g_bounds_l = std::make_unique<double[]>(number_constraints);
    g_bounds_u = std::make_unique<double[]>(number_constraints);

    // just standard bounds, but checking for x0_fixed or xf_fixed
    for (int x_index = 0; x_index < off_x; x_index++) {
        x_bounds_l[x_index] = problem->x0_fixed[x_index] ? *problem->x0_fixed[x_index] : problem->x_bounds[x_index].LB;
        x_bounds_u[x_index] = problem->x0_fixed[x_index] ? *problem->x0_fixed[x_index] : problem->x_bounds[x_index].UB;
    }

    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (i == mesh->intervals - 1 && j == mesh->nodes[i] - 1) {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_bounds_l[off_acc_xu[i][j] + x_index] = problem->xf_fixed[x_index] ? *problem->xf_fixed[x_index] : problem->x_bounds[x_index].LB;
                    x_bounds_u[off_acc_xu[i][j] + x_index] = problem->xf_fixed[x_index] ? *problem->xf_fixed[x_index] : problem->x_bounds[x_index].UB;
                }
            } 
            else {
                for (int x_index = 0; x_index < off_x; x_index++) {
                    x_bounds_l[off_acc_xu[i][j] + x_index] = problem->x_bounds[x_index].LB;
                    x_bounds_u[off_acc_xu[i][j] + x_index] = problem->x_bounds[x_index].UB;
                }
            }
            for (int u_index = 0; u_index < off_u; u_index++) {
                x_bounds_l[off_acc_xu[i][j] + off_x + u_index] = problem->u_bounds[u_index].LB;
                x_bounds_u[off_acc_xu[i][j] + off_x + u_index] = problem->u_bounds[u_index].UB;
            }
        }
    }

    for (int p_index = 0; p_index < off_p; p_index++) {
        x_bounds_l[off_xu_total + p_index] = problem->p_bounds[p_index].LB;
        x_bounds_u[off_xu_total + p_index] = problem->p_bounds[p_index].UB;
    }

    // standard constraint bounds
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full.f_size; f_index++) {
                g_bounds_l[mesh->acc_nodes[i][j] + f_index] = 0;
                g_bounds_u[mesh->acc_nodes[i][j] + f_index] = 0;
            }
            for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
                g_bounds_l[mesh->acc_nodes[i][j] + problem->full.f_size + g_index] = problem->full.g_bounds[g_index].LB;
                g_bounds_u[mesh->acc_nodes[i][j] + problem->full.f_size + g_index] = problem->full.g_bounds[g_index].UB;
            }
        }
    }

    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {
        g_bounds_l[mesh->node_count + r_index] = problem->boundary.r_bounds[r_index].LB;
        g_bounds_u[mesh->node_count + r_index] = problem->boundary.r_bounds[r_index].UB;
    }
}

void NLP::initInitialPoint() {

}