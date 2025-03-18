#include "nlp.h"

// some convenience macros, define the offsets / number of elements in the callback function data arrays
#define EVAL_OFFSET_IJ problem->full.eval_size * mesh->acc_nodes[i][j]
#define JAC_OFFSET_IJ problem->full.jac_size * mesh->acc_nodes[i][j]
#define HES_OFFSET_IJ problem->full.hes_size * mesh->acc_nodes[i][j]

void NLP::init() {
    initSizesOffsets();
    initBounds();
    initStartingPoint();
    initJacobian();
    initHessian();
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

void NLP::initJacobian() {
    calculateJacobianNonzeros();
    initJacobianSparsityPattern();
}

void NLP::calculateJacobianNonzeros() {
    // stage 1: calculate nnz of blocks and number of collisions, where df_k / dx_k != 0. these are contained by default because of the D-Matrix
    int nnz_r = 0;
    int nnz_f_g = 0;
    int diagonal_collisions = 0;
    for (int f_index = 0; f_index < problem->full.f_size; f_index++) {
        for (const auto& df_k_dx : problem->full.lfg[problem->full.f_index_start + f_index].jac.dx) {
            if (df_k_dx.index == f_index) {
                diagonal_collisions++;
            }
        }
        nnz_f_g += problem->full.lfg[problem->full.f_index_start + f_index].jac.dx.size() +
                   problem->full.lfg[problem->full.f_index_start + f_index].jac.du.size() +
                   problem->full.lfg[problem->full.f_index_start + f_index].jac.dp.size();
    }
    for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
            nnz_f_g += problem->full.lfg[problem->full.g_index_start + g_index].jac.dx.size() +
                       problem->full.lfg[problem->full.g_index_start + g_index].jac.du.size() +
                       problem->full.lfg[problem->full.g_index_start + g_index].jac.dp.size();
    }

    // nnz of block i can be calculated as m_i * ((m_i + 2) * #f + #g - coll(df_i, dx_i)), where m_i is the number of nodes on that interval
    off_acc_jac.reserve(mesh->intervals + 1);
    off_acc_jac.emplace_back(0);
    for (int i = 0; i < mesh->intervals; i++) {
        off_acc_jac.emplace_back(off_acc_jac[i - 1] + mesh->nodes[i] * ((mesh->nodes[i] + 2) * problem->full.f_size + problem->full.g_size - diagonal_collisions));
    }

    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {
            nnz_r += problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dx0.size() +
                     problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dxf.size() +
                     problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dp.size();
    }

    nnz_jac = off_acc_jac.back() + nnz_r;

    // allocate memory
    curr_jac  = std::make_unique<double[]>(nnz_jac);
    der_jac   = std::make_unique<double[]>(nnz_jac);
    i_row_jac = std::make_unique<int[]>(nnz_jac);
    j_col_jac = std::make_unique<int[]>(nnz_jac);
    std::memset(der_jac.get(), 0, nnz_jac * sizeof(double));
}

void NLP::initJacobianSparsityPattern() {
    // stage 2: calculate the sparsity pattern i_row_jac, j_col_jac and the constant differentiation matrix part der_jac
    for (int i = 0; i < mesh->intervals; i++) {
        int nnz_index = off_acc_jac[i]; // make local var: possible block parallelization
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full.f_size; f_index++) {
                int eqn_index = off_acc_fg[i][j] + f_index;

                // dColl / dx for x_{i-1, m_{i-1}} base point states
                i_row_jac[nnz_index] = eqn_index;
                j_col_jac[nnz_index] = (i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]) + f_index;
                der_jac[nnz_index]   = collocation->D[mesh->nodes[i]][j + 1][0];
                nnz_index++;

                // dColl / dx for x_{i, j} collocation point states
                for (int k = 0; k < mesh->nodes[i]; k++) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][k] + f_index;
                    der_jac[nnz_index]   = collocation->D[mesh->nodes[i]][j + 1][k + 1];
                    nnz_index++;
                }

                // df / dx
                for (auto& df_dx : problem->full.lfg[problem->full.f_index_start + f_index].jac.dx) {
                    if (df_dx.index != f_index) {
                        i_row_jac[nnz_index] = eqn_index;
                        j_col_jac[nnz_index] = off_acc_xu[i][j] + df_dx.index;
                        nnz_index++;
                    }
                }

                // df / du
                for (auto& df_du : problem->full.lfg[problem->full.f_index_start + f_index].jac.du) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + off_x + df_du.index;
                    nnz_index++;
                }

                // df / dp
                for (auto& df_dp : problem->full.lfg[problem->full.f_index_start + f_index].jac.dp) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xu_total + df_dp.index;
                    nnz_index++;
                }
            }

            for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
                int eqn_index = off_acc_fg[i][j] + problem->full.f_size + g_index;

                // dg / dx
                for (auto& dg_dx : problem->full.lfg[problem->full.g_index_start + g_index].jac.dx) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + dg_dx.index;
                    nnz_index++;
                }

                // dg / du
                for (auto& dg_du : problem->full.lfg[problem->full.g_index_start + g_index].jac.du) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + off_x + dg_du.index;
                    nnz_index++;
                }

                // dg / dp
                for (auto& dg_dp : problem->full.lfg[problem->full.g_index_start + g_index].jac.dp) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xu_total+ dg_dp.index;
                    nnz_index++;
                }
            }
        }
        assert(nnz_index == off_acc_jac[i + 1]);
    }

    int nnz_index = off_acc_jac.back();
    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {
        int eqn_index = off_fg_total + r_index;

        // dr / dx0
        for (auto& dr_dx0 : problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dx0) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = dr_dx0.index;
            nnz_index++;
        }

        // dg / dxf
        for (auto& dr_dxf : problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dxf) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_last_xu + dr_dxf.index;
            nnz_index++;
        }

        // dg / dp
        for (auto& dr_dp: problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dp) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_xu_total + dr_dp.index;
            nnz_index++;
        }
    }
    assert(nnz_index == nnz_jac);
}

void NLP::initHessian() {
    initSparsityHessian();
}

void NLP::initSparsityHessian() {
    // takes O(nnz(A) + nnz(B) + ...+ nnz(H))

    // stage 1: calculate IndexSet and nnz
    OrderedIndexSet A, B, C, D, E, F, G, H;
    for (auto& mr : problem->boundary.mr) {
        A.insertSparsity(mr.hes.dx0_dx0, 0, 0);
        C.insertSparsity(mr.hes.dxf_dx0, 0, 0);
        D.insertSparsity(mr.hes.dxf_dxf, 0, 0);
        E.insertSparsity(mr.hes.dp_dx0,  0, 0);
        G.insertSparsity(mr.hes.dp_dxf,  0, 0);
        H.insertSparsity(mr.hes.dp_dp,   0, 0);
    }
    for (auto& lfg : problem->full.lfg) {
        B.insertSparsity(lfg.hes.dx_dx, 0, 0);
        B.insertSparsity(lfg.hes.du_dx, off_x, 0);
        B.insertSparsity(lfg.hes.du_du, off_x, off_x);
        D.insertSparsity(lfg.hes.dx_dx, 0, 0);
        D.insertSparsity(lfg.hes.du_dx, off_x, 0);
        D.insertSparsity(lfg.hes.du_du, off_x, off_x);
        F.insertSparsity(lfg.hes.dp_dx, 0, 0);
        F.insertSparsity(lfg.hes.dp_du, 0, off_x);
        G.insertSparsity(lfg.hes.dp_dx, 0, 0);
        G.insertSparsity(lfg.hes.dp_du, 0, off_x);
        H.insertSparsity(lfg.hes.dp_dp, 0, 0);
    }

    nnz_hes = (B.size() + F.size()) * (mesh->node_count - 1) + A.size() + C.size() + D.size() + E.size() + G.size() + H.size();

    i_row_hes = std::make_unique<int[]>(nnz_hes);
    j_col_hes = std::make_unique<int[]>(nnz_hes);
    curr_hes  = std::make_unique<double[]>(nnz_hes);

    // stage 2: build int** (row, col) -> int index for all block structures
    // can be exact (no offset needed) or non exact (offset for full block or even rowwise needed)
    int hes_nnz_counter = 0;

    // A: exact
    for (auto& [row, col] : A.set) {
        hes_a.insert(row, col, hes_nnz_counter++);
    }

    // B: non exact, thus local counter
    int block_b_nnz = 0;
    for (auto& [row, col] : B.set) {
        hes_b.insert(row, col, block_b_nnz++);
    }
    hes_b.off_prev = hes_a.nnz;  // set size of A block as offset
    hes_nnz_counter += block_b_nnz * (mesh->node_count - 1);

    // C, D: exact with row dependence
    int c_index = 0;
    int d_index = 0;
    std::vector<std::pair<int, int>> C_flat(C.set.begin(), C.set.end());
    std::vector<std::pair<int, int>> D_flat(D.set.begin(), D.set.end());
    for (int x_index = 0; x_index < problem->x_size; x_index++) {
        while (C_flat[c_index].first == x_index) {
            hes_c.insert(C_flat[c_index].first, C_flat[c_index].second, hes_nnz_counter++);
            c_index++;
        }
        while (D_flat[d_index].first == x_index) {
            hes_d.insert(D_flat[d_index].first, D_flat[d_index].second, hes_nnz_counter++);
            d_index++;
        }
    }
    for (;d_index < D_flat.size(); d_index++) {
        hes_d.insert(D_flat[d_index].first, D_flat[d_index].second, hes_nnz_counter++);
    }

    // E, F, G, H: partially exact (all except F) with row dependence
    int e_index = 0;
    int f_index = 0;
    int g_index = 0;
    int h_index = 0;
    std::vector<std::pair<int, int>> E_flat(E.set.begin(), E.set.end());
    std::vector<std::pair<int, int>> F_flat(F.set.begin(), F.set.end());
    std::vector<std::pair<int, int>> G_flat(G.set.begin(), G.set.end());
    std::vector<std::pair<int, int>> H_flat(H.set.begin(), H.set.end());
    for (int p_index = 0; p_index < problem->p_size; p_index++) {
        while (E_flat[e_index].first == p_index) {
            hes_e.insert(E_flat[e_index].first, E_flat[e_index].second, hes_nnz_counter++);
            e_index++;
        }

        int row_f_nnz = 0;
        hes_f.row_offset_prev[p_index] = hes_nnz_counter; // E_{p_index, :} offset
        while (F_flat[f_index].first == p_index) {
            hes_f.insert(F_flat[f_index].first, F_flat[f_index].second, row_f_nnz++);
            f_index++;
        }
        hes_f.row_size[p_index] = row_f_nnz; // F_{p_index, :} size -> offset for next F blocks
        hes_nnz_counter += (mesh->node_count - 1) * row_f_nnz;

        while (G_flat[g_index].first == p_index) {
            hes_g.insert(G_flat[g_index].first, G_flat[g_index].second, hes_nnz_counter++);
            g_index++;
        }

        while (H_flat[h_index].first == p_index) {
            hes_h.insert(H_flat[h_index].first, H_flat[h_index].second, hes_nnz_counter++);
            h_index++;
        }
    }
    assert(hes_nnz_counter == nnz_hes);
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

 void NLP::eval_jac_g_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (evaluation_state.jac_g) {
        return;
    }
    else {
        callback_jacobian();
        eval_jac_g();
    }
 }

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
                    curr_grad[off_acc_xu[i][j] + dL_dx.index] = mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(dL_dx.value + JAC_OFFSET_IJ));
                }
                for (auto& dL_du : problem->full.lfg[0].jac.du) {
                    curr_grad[off_acc_xu[i][j] + off_x + dL_du.index] = mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(dL_du.value + JAC_OFFSET_IJ));
                }
                for (auto& dL_dp : problem->full.lfg[0].jac.dp) {
                    curr_grad[off_xu_total + dL_dp.index] += mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(dL_dp.value + JAC_OFFSET_IJ));
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
                                              &curr_x[i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]],  // x_{i-1, m_{i-1}} base point states
                                              &curr_x[off_acc_xu[i][0]],                                        // collocation point states
                                              &curr_g[off_acc_fg[i][0]]);                                       // constraint start index 
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

void NLP::eval_jac_g() {
    // copy constant part into curr_jac
    std::memcpy(curr_jac.get(), der_jac.get(), nnz_jac * sizeof(double));

    for (int i = 0; i < mesh->intervals; i++) {
        int nnz_index = off_acc_jac[i]; // make local var: possible block parallelization
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full.f_size; f_index++) {
                // index of possible collision in df_i / dx_i
                int collision_index = nnz_index + j + 1;

                // offset of constant differentiation matrix part
                nnz_index += mesh->nodes[i] + 1;

                // df / dx
                for (auto& df_dx : problem->full.lfg[problem->full.f_index_start + f_index].jac.dx) {
                    if (df_dx.index != f_index) {
                        curr_jac[nnz_index++] = -mesh->delta_t[i] * (*(df_dx.value + JAC_OFFSET_IJ));
                    } 
                    else {
                        curr_jac[collision_index] -= mesh->delta_t[i] * (*(df_dx.value + JAC_OFFSET_IJ));
                    }
                }

                // df / du
                for (auto& df_du : problem->full.lfg[problem->full.f_index_start + f_index].jac.du) {
                    curr_jac[nnz_index++] = -mesh->delta_t[i] * (*(df_du.value + JAC_OFFSET_IJ));
                }

                // df / dp
                for (auto& df_dp : problem->full.lfg[problem->full.f_index_start + f_index].jac.dp) {
                    curr_jac[nnz_index++] = -mesh->delta_t[i] * (*(df_dp.value + JAC_OFFSET_IJ));
                }
            }

            for (int g_index = 0; g_index < problem->full.g_size; g_index++) {
                // dg / dx
                for (auto& dg_dx : problem->full.lfg[problem->full.g_index_start + g_index].jac.dx) {
                    curr_jac[nnz_index++] = (*(dg_dx.value + JAC_OFFSET_IJ));
                }

                // dg / du
                for (auto& dg_du : problem->full.lfg[problem->full.g_index_start + g_index].jac.du) {
                    curr_jac[nnz_index++] = (*(dg_du.value + JAC_OFFSET_IJ));
                }

                // dg / dp
                for (auto& dg_dp : problem->full.lfg[problem->full.g_index_start + g_index].jac.dp) {
                    curr_jac[nnz_index++] = (*(dg_dp.value + JAC_OFFSET_IJ));
                }
            }
        }
        assert(nnz_index == off_acc_jac[i + 1]);
    }

    int nnz_index = off_acc_jac.back();
    for (int r_index = 0; r_index < problem->boundary.r_size; r_index++) {

        // dr / dx0
        for (auto& dr_dx0 : problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dx0) {
            curr_jac[nnz_index++] = (*(dr_dx0.value));
        }

        // dg / dxf
        for (auto& dr_dxf : problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dxf) {
            curr_jac[nnz_index++] = (*(dr_dxf.value));
        }

        // dg / dp
        for (auto& dr_dp: problem->boundary.mr[problem->boundary.r_index_start + r_index].jac.dp) {
            curr_jac[nnz_index++] = (*(dr_dp.value));
        }
    }
    assert(nnz_index == nnz_jac);
};

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
