#include "gdop.h"

// some convenience macros, define the offsets / number of elements in the callback function data arrays
#define EVAL_OFFSET_IJ problem->full->eval_size * mesh->acc_nodes[i][j]
#define JAC_OFFSET_IJ problem->full->jac_size * mesh->acc_nodes[i][j]
#define HES_OFFSET_IJ problem->full->hes_size * mesh->acc_nodes[i][j]

void GDOP::init() {
    initSizesOffsets();
    initBuffers();
    initBounds();
    initJacobian();
    initHessian();
    initStartingPoint();
}

void GDOP::initSizesOffsets() {
    off_x = problem->x_size;
    off_u = problem->u_size;
    off_p = problem->p_size;
    off_xu = off_x + off_u;
    off_acc_xu = mesh->createAccOffsetXU(off_x, off_xu);          // variables  x_ij offset
    off_last_xu = off_acc_xu.back().back();                       // variables final grid point x_ij
    off_xu_total = off_last_xu + off_xu;                          // first parameter
    number_vars = off_xu_total + problem->p_size;
    off_acc_fg = mesh->createAccOffsetFG(problem->full->fg_size);  // constraint f_ij offset
    off_fg_total = mesh->node_count * problem->full->fg_size;      // constraint r_0 offset
    number_constraints = problem->boundary->r_size + off_fg_total;
}

void GDOP::initBuffers() {
    curr_x      = FixedVector<double>(number_vars);
    init_x      = FixedVector<double>(number_vars);
    curr_grad   = FixedVector<double>(number_vars);
    x_lb        = FixedVector<double>(number_vars);
    x_ub        = FixedVector<double>(number_vars);
    curr_lambda = FixedVector<double>(number_constraints);
    curr_g      = FixedVector<double>(number_constraints);
    g_lb        = FixedVector<double>(number_constraints);
    g_ub        = FixedVector<double>(number_constraints);
}

void GDOP::initBounds() {
    // TODO: perform scaling before the bounds are set!

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
            for (int f_index = 0; f_index < problem->full->f_size; f_index++) {
                g_lb[off_acc_fg[i][j] + f_index] = 0;
                g_ub[off_acc_fg[i][j] + f_index] = 0;
            }
            for (int g_index = 0; g_index < problem->full->g_size; g_index++) {
                g_lb[off_acc_fg[i][j] + problem->full->f_size + g_index] = problem->full->g_bounds[g_index].lb;
                g_ub[off_acc_fg[i][j] + problem->full->f_size + g_index] = problem->full->g_bounds[g_index].ub;
            }
        }
    }

    for (int r_index = 0; r_index < problem->boundary->r_size; r_index++) {
        g_lb[off_fg_total + r_index] = problem->boundary->r_bounds[r_index].lb;
        g_ub[off_fg_total + r_index] = problem->boundary->r_bounds[r_index].ub;
    }
}

void GDOP::initStartingPoint() {
    Trajectory new_guess = guess->interpolate(*mesh, *collocation);

    for (int x_index = 0; x_index < off_x; x_index++) {
        init_x[x_index] = new_guess.x[x_index][0];
    }

    int index = 1;
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int x_index = 0; x_index < off_x; x_index++) {
                init_x[off_acc_xu[i][j] + x_index] = new_guess.x[x_index][index];
            }
            for (int u_index = 0; u_index < off_u; u_index++) {
                init_x[off_acc_xu[i][j] + off_x + u_index] = new_guess.u[u_index][index];
            }
            index++;
        }
    }
    for (int p_index = 0; p_index < off_p; p_index++) {
        init_x[off_xu_total + p_index] = new_guess.p[p_index];
    }
}

void GDOP::initJacobian() {
    initJacobianNonzeros();
    initJacobianSparsityPattern();
}

void GDOP::initJacobianNonzeros() {
    // stage 1: calculate nnz of blocks and number of collisions, where df_k / dx_k != 0. these are contained by default because of the D-Matrix
    int nnz_f = 0;
    int nnz_g = 0;
    int nnz_r = 0;
    int diagonal_collisions = 0;
    for (int f_index = 0; f_index < problem->full->f_size; f_index++) {
        for (const auto& df_k_dx : problem->full->lfg[problem->full->f_index_start + f_index].jac.dx) {
            if (df_k_dx.col == f_index) {
                diagonal_collisions++;
            }
        }
        nnz_f += problem->full->lfg[problem->full->f_index_start + f_index].jac.nnz();
    }
    for (int g_index = 0; g_index < problem->full->g_size; g_index++) {
            nnz_g += problem->full->lfg[problem->full->g_index_start + g_index].jac.nnz();
    }

    // nnz of block i can be calculated as m_i * ((m_i + 2) * #f + #g - coll(df_i, dx_i)), where m_i is the number of nodes on that interval
    off_acc_jac_fg = FixedVector<int>(mesh->intervals + 1);
    for (int i = 0; i < mesh->intervals; i++) {
        off_acc_jac_fg[i+1] = off_acc_jac_fg[i] + mesh->nodes[i] * ((mesh->nodes[i] + 1) * off_x + nnz_f + nnz_g - diagonal_collisions);
    }

    for (int r_index = 0; r_index < problem->boundary->r_size; r_index++) {
            nnz_r += problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dx0.size() +
                     problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dxf.size() +
                     problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dp.size();
    }

    nnz_jac = off_acc_jac_fg.back() + nnz_r;

    // allocate memory
    curr_jac  = FixedVector<double>(nnz_jac);
    der_jac   = FixedVector<double>(nnz_jac);
    i_row_jac = FixedVector<int>(nnz_jac);
    j_col_jac = FixedVector<int>(nnz_jac);
}

void GDOP::initJacobianSparsityPattern() {
    // stage 2: calculate the sparsity pattern i_row_jac, j_col_jac and the constant differentiation matrix part der_jac
    for (int i = 0; i < mesh->intervals; i++) {
        int nnz_index = off_acc_jac_fg[i]; // make local var: possible block parallelization
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full->f_size; f_index++) {
                int eqn_index = off_acc_fg[i][j] + f_index;

                // dColl / dx for x_{i-1, m_{i-1}} base point states (k = -1, prev state)
                i_row_jac[nnz_index] = eqn_index;
                j_col_jac[nnz_index] = (i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]) + f_index;
                der_jac[nnz_index]   = collocation->D[mesh->nodes[i]][j + 1][0];
                nnz_index++;

                // dColl / dx for x_{i, j} collocation point states up to the collision
                // this means the derivative matrix linear combinations x_ik have k < j
                for (int k = 0; k < j; k++) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][k] + f_index;
                    der_jac[nnz_index]   = collocation->D[mesh->nodes[i]][j + 1][k + 1];
                    nnz_index++;
                }

                // df / dx
                for (auto& df_dx : problem->full->lfg[problem->full->f_index_start + f_index].jac.dx) {
                    if (df_dx.col != f_index) {
                        i_row_jac[nnz_index] = eqn_index;
                        j_col_jac[nnz_index] = off_acc_xu[i][j] + df_dx.col;
                        nnz_index++;
                    }
                    else {
                        // handle the diagonal collision of the diagonal jacobian block
                        // this means the derivative matrix linear combinations x_ik have k == j and df_k / dx_k != 0
                        i_row_jac[nnz_index] = eqn_index;
                        j_col_jac[nnz_index] = off_acc_xu[i][j] + df_dx.col; // here df_dx and f_index collide!! (could also be written with dx_dx = f_index), which is nz element of the derivaitve matrix part
                        der_jac[nnz_index] = collocation->D[mesh->nodes[i]][j + 1][j + 1];
                        nnz_index++;
                    }
                }

                // df / du
                for (auto& df_du : problem->full->lfg[problem->full->f_index_start + f_index].jac.du) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + off_x + df_du.col;
                    nnz_index++;
                }

                // dColl / dx for x_{i, j} collocation point states after the collision / diagonal block in block jacobian
                // this means the derivative matrix linear combinations x_ik have k > j
                for (int k = j + 1; k < mesh->nodes[i]; k++) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][k] + f_index;
                    der_jac[nnz_index]   = collocation->D[mesh->nodes[i]][j + 1][k + 1];
                    nnz_index++;
                }

                // df / dp
                for (auto& df_dp : problem->full->lfg[problem->full->f_index_start + f_index].jac.dp) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xu_total + df_dp.col;
                    nnz_index++;
                }
            }

            for (int g_index = 0; g_index < problem->full->g_size; g_index++) {
                int eqn_index = off_acc_fg[i][j] + problem->full->f_size + g_index;

                // dg / dx
                for (auto& dg_dx : problem->full->lfg[problem->full->g_index_start + g_index].jac.dx) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + dg_dx.col;
                    nnz_index++;
                }

                // dg / du
                for (auto& dg_du : problem->full->lfg[problem->full->g_index_start + g_index].jac.du) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_acc_xu[i][j] + off_x + dg_du.col;
                    nnz_index++;
                }

                // dg / dp
                for (auto& dg_dp : problem->full->lfg[problem->full->g_index_start + g_index].jac.dp) {
                    i_row_jac[nnz_index] = eqn_index;
                    j_col_jac[nnz_index] = off_xu_total + dg_dp.col;
                    nnz_index++;
                }
            }
        }
        assert(nnz_index == off_acc_jac_fg[i + 1]);
    }

    int nnz_index = off_acc_jac_fg.back();
    for (int r_index = 0; r_index < problem->boundary->r_size; r_index++) {
        int eqn_index = off_fg_total + r_index;

        // dr / dx0
        for (auto& dr_dx0 : problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dx0) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = dr_dx0.col;
            nnz_index++;
        }

        // dg / dxf
        for (auto& dr_dxf : problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dxf) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_last_xu + dr_dxf.col;
            nnz_index++;
        }

        // dg / dp
        for (auto& dr_dp: problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dp) {
            i_row_jac[nnz_index] = eqn_index;
            j_col_jac[nnz_index] = off_xu_total + dr_dp.col;
            nnz_index++;
        }
    }
    assert(nnz_index == nnz_jac);
}

void GDOP::initHessian() {
    // takes O(nnz(A) + nnz(B) + ...+ nnz(H)) for creation of ** Maps and O(nnz(Hessian)) for creation of Hessian sparsity pattern

    // stage 1: calculate IndexSet and nnz
    OrderedIndexSet A, B, C, D, E, F, G, H;
    for (auto& mr : problem->boundary->mr) {
        A.insertSparsity(mr.hes.dx0_dx0, 0, 0);
        C.insertSparsity(mr.hes.dxf_dx0, 0, 0);
        D.insertSparsity(mr.hes.dxf_dxf, 0, 0);
        E.insertSparsity(mr.hes.dp_dx0,  0, 0);
        G.insertSparsity(mr.hes.dp_dxf,  0, 0);
        H.insertSparsity(mr.hes.dp_dp,   0, 0);
    }
    for (auto& lfg : problem->full->lfg) {
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

    i_row_hes = FixedVector<int>(nnz_hes);
    j_col_hes = FixedVector<int>(nnz_hes);
    curr_hes  = FixedVector<double>(nnz_hes);

    // stage 2: build int** (row, col) -> int index for all block structures
    // can be exact (no offset needed) or non exact (offset for full block or even rowwise needed)
    // also init sparsity pattern (i_row_hes, j_col_hes) for all exact blocks
    int hes_nnz_counter = 0;

    // A: exact
    for (auto& [row, col] : A.set) {
        i_row_hes[hes_nnz_counter] = row; // x_0
        j_col_hes[hes_nnz_counter] = col; // x_0
        hes_a.insert(row, col, hes_nnz_counter++);
    }

    // B: non exact, thus local counter
    int block_b_nnz = 0;
    for (auto& [row, col] : B.set) {
        hes_b.insert(row, col, block_b_nnz++);
    }
    hes_b.off_prev = hes_a.nnz;  // set size of A block as offset

    // init B hessian pattern O(node_count * nnz(L_{xu, xu} ∪ f_{xu, xu} ∪ g_{xu, xu})) - expensive, parallel execution should be possible
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (i != mesh->intervals - 1 && j != mesh->nodes[mesh->intervals - 1] - 1) {
                for (auto& [row, col] : B.set) {
                    int xu_hes_index = hes_b.access(row, col, mesh->acc_nodes[i][j]);
                    i_row_hes[xu_hes_index] = off_acc_xu[i][j] + row; // xu_{ij}
                    j_col_hes[xu_hes_index] = off_acc_xu[i][j] + col; // xu_{ij}
                }
            }
        }
    }

    hes_nnz_counter += block_b_nnz * (mesh->node_count - 1);

    // C, D: exact with row dependence
    int c_index = 0;
    int d_index = 0;
    FixedVector<std::pair<int, int>> C_flat(C.set.begin(), C.set.end());
    FixedVector<std::pair<int, int>> D_flat(D.set.begin(), D.set.end());
    for (int x_index = 0; x_index < problem->x_size; x_index++) {
        while (c_index < C_flat.size() && C_flat[c_index].first == x_index) {
            i_row_hes[hes_nnz_counter] = off_last_xu + C_flat[c_index].first; // x_{nm}
            j_col_hes[hes_nnz_counter] = C_flat[c_index].second;              // x_0
            hes_c.insert(C_flat[c_index].first, C_flat[c_index].second, hes_nnz_counter++);
            c_index++;
        }
        while (d_index < D_flat.size() && D_flat[d_index].first == x_index) {
            i_row_hes[hes_nnz_counter] = off_last_xu + D_flat[d_index].first;  // x_{nm}
            j_col_hes[hes_nnz_counter] = off_last_xu + D_flat[d_index].second; // x_{nm}
            hes_d.insert(D_flat[d_index].first, D_flat[d_index].second, hes_nnz_counter++);
            d_index++;
        }
    }
    for (;d_index < D_flat.size(); d_index++) {
        i_row_hes[hes_nnz_counter] = off_last_xu + D_flat[d_index].first;  // u_{nm}
        j_col_hes[hes_nnz_counter] = off_last_xu + D_flat[d_index].second; // xu_{nm}
        hes_d.insert(D_flat[d_index].first, D_flat[d_index].second, hes_nnz_counter++);
    }

    // E, F, G, H: partially exact (all except F) with row dependence
    int e_index = 0;
    int f_index = 0;
    int g_index = 0;
    int h_index = 0;
    FixedVector<std::pair<int, int>> E_flat(E.set.begin(), E.set.end());
    FixedVector<std::pair<int, int>> F_flat(F.set.begin(), F.set.end());
    FixedVector<std::pair<int, int>> G_flat(G.set.begin(), G.set.end());
    FixedVector<std::pair<int, int>> H_flat(H.set.begin(), H.set.end());
    for (int p_index = 0; p_index < problem->p_size; p_index++) {
        while (e_index < E_flat.size() && E_flat[e_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + E_flat[e_index].first; // p
            j_col_hes[hes_nnz_counter] = E_flat[e_index].second;               // x_0
            hes_e.insert(E_flat[e_index].first, E_flat[e_index].second, hes_nnz_counter++);
            e_index++;
        }
        

        int row_f_nnz = 0;
        hes_f.row_offset_prev[p_index] = hes_nnz_counter; // E_{p_index, :} offset
        while (f_index < F_flat.size() && F_flat[f_index].first == p_index) {
            hes_f.insert(F_flat[f_index].first, F_flat[f_index].second, row_f_nnz++);
            f_index++;
        }
        
        hes_f.row_size[p_index] = row_f_nnz; // F_{p_index, :} size -> offset for next F blocks
        hes_nnz_counter += (mesh->node_count - 1) * row_f_nnz;
        while (f_index < F_flat.size() && G_flat[g_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + G_flat[g_index].first; // p
            j_col_hes[hes_nnz_counter] = off_last_xu + G_flat[g_index].second; // xu_{nm}
            hes_g.insert(G_flat[g_index].first, G_flat[g_index].second, hes_nnz_counter++);
            g_index++;
        }
        
        while (g_index < H_flat.size() && H_flat[h_index].first == p_index) {
            i_row_hes[hes_nnz_counter] = off_xu_total + G_flat[g_index].first;  // p
            j_col_hes[hes_nnz_counter] = off_xu_total + G_flat[g_index].second; // p
            hes_h.insert(H_flat[h_index].first, H_flat[h_index].second, hes_nnz_counter++);
            h_index++;
        }
        
    }
    assert(hes_nnz_counter == nnz_hes);

    // init F hessian pattern O(node_count * nnz(L_{p, xu} ∪ f_{p, xu} ∪ g_{p, xu})) - expensive, parallel execution should be possible
    for (int i = 0; i < mesh->intervals; i++) {
        for (int j = 0; j < mesh->nodes[i]; j++) {
            if (i != mesh->intervals - 1 && j != mesh->nodes[mesh->intervals - 1] - 1) {
                for (auto& [row, col] : F.set) {
                    int xu_hes_index = hes_b.access(row, col, mesh->acc_nodes[i][j]);
                    i_row_hes[xu_hes_index] = off_xu_total + row;     // p
                    j_col_hes[xu_hes_index] = off_acc_xu[i][j] + col; // xu_{ij}
                }
            }
        }
    }
}

/* nlp function evaluations happen in two stages:
 * 1. fill the input buffers -> evaluate all continuous callback functions and fill the output buffers
 * 2. evaluate the nlp function by accessing the structures and double* defined in problem, simply use the filled buffers
 * 
 * since all functions are evaluated in step 2, this can always be executed in parallel even if the callbacks were not parallel!
 * because of this structure, before every nlp evaluation check_new_x has to be performed and step 1 has to be executed in case
 */

void GDOP::check_new_x(const double* nlp_solver_x, bool new_x) {
    evaluation_state.check_reset_x(new_x);
    if (!evaluation_state.x_set_unscaled) {
        // Scaler.scale(nlp_solver_x, curr_x), perform scaling here, memcpy nlp_solver_x -> unscaled -> scale
        // rn we just memset the data to curr_x
        curr_x.assign(nlp_solver_x, number_vars);
        evaluation_state.x_set_unscaled = true;
    }
}

void GDOP::check_new_lambda(const double* nlp_solver_lambda, const bool new_lambda) {
    evaluation_state.check_reset_lambda(new_lambda);
    if (!evaluation_state.lambda_set) {
        curr_lambda.assign(nlp_solver_lambda, number_constraints);
        evaluation_state.lambda_set = true;
    }
}

void GDOP::check_new_sigma(const double sigma_f) {
    if (sigma_f != curr_sigma_f) {
        curr_sigma_f = sigma_f;
        evaluation_state.hes_lag = false;
    }
}

void GDOP::eval_f_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (!evaluation_state.eval_f) {
        callback_evaluation();
    }
    eval_f();
}

void GDOP::eval_g_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (!evaluation_state.eval_g) {
        callback_evaluation();
    }
    eval_g();
}

void GDOP::eval_grad_f_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (!evaluation_state.grad_f) {
        callback_jacobian();
    }
    eval_grad_f();
};

 void GDOP::eval_jac_g_safe(const double* nlp_solver_x, bool new_x) {
    check_new_x(nlp_solver_x, new_x);
    if (!evaluation_state.jac_g) {
        callback_jacobian();
    }
    eval_jac_g();
 }

 void GDOP::eval_hes_safe(const double* nlp_solver_x, const double* nlp_solver_lambda, double sigma_f, bool new_x, bool new_lambda) {
    check_new_x(nlp_solver_x, new_x);
    check_new_lambda(nlp_solver_lambda, new_lambda);
    check_new_sigma(sigma_f);
    if (!evaluation_state.hes_lag) {
        callback_hessian();
    }
    eval_hes();
 }

void GDOP::eval_f() {
    double mayer = 0;
    if (problem->boundary->has_mayer) {
        mayer = problem->boundary->getEvalM();
    };

    double lagrange = 0;
    if (problem->full->has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                lagrange += mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * problem->full->getEvalL(mesh->acc_nodes[i][j]);
            }
        }
    }

    curr_obj = mayer + lagrange;
}


void GDOP::eval_grad_f() {
    curr_grad.fill_zero();
    if (problem->full->has_lagrange) {
        for (int i = 0; i < mesh->intervals; i++) {
            for (int j = 0; j < mesh->nodes[i]; j++) {
                for (auto& dL_dx : problem->full->lfg[0].jac.dx) {
                    curr_grad[off_acc_xu[i][j] + dL_dx.col] = mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(dL_dx.value + JAC_OFFSET_IJ));
                }
                for (auto& dL_du : problem->full->lfg[0].jac.du) {
                    curr_grad[off_acc_xu[i][j] + off_x + dL_du.col] = mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(dL_du.value + JAC_OFFSET_IJ));
                }
                for (auto& dL_dp : problem->full->lfg[0].jac.dp) {
                    curr_grad[off_xu_total + dL_dp.col] += mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j] * (*(dL_dp.value + JAC_OFFSET_IJ));
                }
            }
        }
    }
    if (problem->boundary->has_mayer) {
        for (auto& dM_dx0 : problem->boundary->mr[0].jac.dx0) {
            curr_grad[dM_dx0.col] = (*dM_dx0.value);
        }
        for (auto& dM_dxf : problem->boundary->mr[0].jac.dxf) {
            curr_grad[off_last_xu + dM_dxf.col] += (*dM_dxf.value);
        }
        for (auto& dM_dp : problem->boundary->mr[0].jac.dp) {
            curr_grad[off_xu_total + dM_dp.col] += (*dM_dp.value);
        }
    }
};

void GDOP::eval_g() {
    curr_g.fill_zero();
    for (int i = 0; i < mesh->intervals; i++) {
        for (int f = problem->full->f_index_start; f < problem->full->f_index_end; f++) {
            collocation->diff_matrix_multiply(mesh->nodes[i], off_x, off_xu, problem->full->fg_size,
                                              &curr_x[i == 0 ? 0 : off_acc_xu[i - 1][mesh->nodes[i - 1] - 1]],  // x_{i-1, m_{i-1}} base point states
                                              &curr_x[off_acc_xu[i][0]],                                        // collocation point states
                                              &curr_g[off_acc_fg[i][0]]);                                       // constraint start index 
        }
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full->f_size; f_index++) {
                curr_g[off_acc_fg[i][j] + f_index] -= mesh->delta_t[i] * problem->full->getEvalF(f_index, mesh->acc_nodes[i][j]);
            }
            for (int g_index = 0; g_index < problem->full->g_size; g_index++) {
                curr_g[off_acc_fg[i][j] + problem->full->f_size + g_index] += problem->full->getEvalG(g_index, mesh->acc_nodes[i][j]);
            }
        }
    }
    for (int r_index = 0; r_index < problem->boundary->r_size; r_index++) {
        curr_g[off_fg_total + r_index] = problem->boundary->getEvalR(r_index);
    }
}

void GDOP::eval_jac_g() {
    // copy constant part into curr_jac
    curr_jac = der_jac;

    for (int i = 0; i < mesh->intervals; i++) {
        int nnz_index = off_acc_jac_fg[i]; // make local var: possible block parallelization
        for (int j = 0; j < mesh->nodes[i]; j++) {
            for (int f_index = 0; f_index < problem->full->f_size; f_index++) {
                // offset for all leading diagonal matrix blocks: d_{jk} * I at t_ik with j < k
                nnz_index += j + 1;

                // df / dx, j == k collision is included implicitly (thus, use -= and not = - like for the others)
                for (auto& df_dx : problem->full->lfg[problem->full->f_index_start + f_index].jac.dx) {
                    curr_jac[nnz_index++] -= mesh->delta_t[i] * (*(df_dx.value + JAC_OFFSET_IJ));
                }

                // df / du
                for (auto& df_du : problem->full->lfg[problem->full->f_index_start + f_index].jac.du) {
                    curr_jac[nnz_index++] = -mesh->delta_t[i] * (*(df_du.value + JAC_OFFSET_IJ));
                }

                // offset for all remaining diagonal matrix blocks: d_{jk} * I at t_ik with k > j
                 nnz_index += mesh->nodes[i] - j - 1;

                // df / dp
                for (auto& df_dp : problem->full->lfg[problem->full->f_index_start + f_index].jac.dp) {
                    curr_jac[nnz_index++] = -mesh->delta_t[i] * (*(df_dp.value + JAC_OFFSET_IJ));
                }
            }

            for (int g_index = 0; g_index < problem->full->g_size; g_index++) {
                // dg / dx
                for (auto& dg_dx : problem->full->lfg[problem->full->g_index_start + g_index].jac.dx) {
                    curr_jac[nnz_index++] = (*(dg_dx.value + JAC_OFFSET_IJ));
                }

                // dg / du
                for (auto& dg_du : problem->full->lfg[problem->full->g_index_start + g_index].jac.du) {
                    curr_jac[nnz_index++] = (*(dg_du.value + JAC_OFFSET_IJ));
                }

                // dg / dp
                for (auto& dg_dp : problem->full->lfg[problem->full->g_index_start + g_index].jac.dp) {
                    curr_jac[nnz_index++] = (*(dg_dp.value + JAC_OFFSET_IJ));
                }
            }
        }
        assert(nnz_index == off_acc_jac_fg[i + 1]);
    }

    int nnz_index = off_acc_jac_fg.back();
    for (int r_index = 0; r_index < problem->boundary->r_size; r_index++) {

        // dr / dx0
        for (auto& dr_dx0 : problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dx0) {
            curr_jac[nnz_index++] = (*(dr_dx0.value));
        }

        // dg / dxf
        for (auto& dr_dxf : problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dxf) {
            curr_jac[nnz_index++] = (*(dr_dxf.value));
        }

        // dg / dp
        for (auto& dr_dp: problem->boundary->mr[problem->boundary->r_index_start + r_index].jac.dp) {
            curr_jac[nnz_index++] = (*(dr_dp.value));
        }
    }
    assert(nnz_index == nnz_jac);
};

void GDOP::eval_hes() {
    curr_hes.fill_zero();
    int lagrange_offset = (int) problem->full->has_lagrange;  // correction for array accesses
    int mayer_offset    = (int) problem->boundary->has_mayer; // correction for array accesses
    for (int i = 0; i < mesh->intervals; i++) {
        // insert buffer for parallel execution here / pass the buffer in updateHessian(double*) calls
        for (int j = 0; j < mesh->nodes[i]; j++) {
            // make sure to be in the right ptr_map region, B and F are only valid for i,j != n,m
            //                                              D and G are only valid for i,j == n,m
            const BlockSparsity* ptr_map_xu_xu;
            const BlockSparsity* ptr_map_p_xu;
            if (i != mesh->intervals - 1 && j != mesh->nodes[mesh->intervals - 1] - 1) {
                ptr_map_xu_xu = &hes_b;
                ptr_map_p_xu  = &hes_f;
            }
            else {
                ptr_map_xu_xu = &hes_d;
                ptr_map_p_xu  = &hes_g;
            }
            if (problem->full->has_lagrange) {
                updateHessianLFG(curr_hes, problem->full->lfg[0].hes, i, j, ptr_map_xu_xu, ptr_map_p_xu,
                                 curr_sigma_f * mesh->delta_t[i] * collocation->b[mesh->nodes[i]][j]);
            }
            for (int f_index = problem->full->f_index_start; f_index < problem->full->f_index_end; f_index++) {
                updateHessianLFG(curr_hes, problem->full->lfg[f_index].hes, i, j, ptr_map_xu_xu, ptr_map_p_xu,
                                 -curr_lambda[off_acc_fg[i][j] + f_index - lagrange_offset] * mesh->delta_t[i]);
            }
            for (int g_index = problem->full->g_index_start; g_index < problem->full->g_index_end; g_index++) {
                updateHessianLFG(curr_hes, problem->full->lfg[g_index].hes, i, j, ptr_map_xu_xu, ptr_map_p_xu,
                                 curr_lambda[off_acc_fg[i][j] + g_index - lagrange_offset]);
            }
        }
    }
    if (problem->boundary->has_mayer) {
        updateHessianMR(curr_hes, problem->boundary->mr[0].hes, curr_sigma_f);
    }

    for (int r_index = problem->boundary->r_index_start; r_index < problem->boundary->r_index_end; r_index++) {
        updateHessianMR(curr_hes, problem->boundary->mr[r_index].hes, curr_lambda[off_fg_total + r_index - mayer_offset]);
    }
}

void GDOP::updateHessianLFG(FixedVector<double>& values, const HessianLFG& hes, const int i, const int j, const BlockSparsity* ptr_map_xu_xu,
                           const BlockSparsity* ptr_map_p_xu, const double factor) {
    const int block_count = mesh->acc_nodes[i][j];
    for (const auto& dx_dx : hes.dx_dx) {
        values[ptr_map_xu_xu->access(dx_dx.row, dx_dx.col, block_count)] += factor * (*(dx_dx.value + HES_OFFSET_IJ));
    }
    for (const auto& du_dx : hes.du_dx) {
        values[ptr_map_xu_xu->access(off_x + du_dx.row, du_dx.col, block_count)] += factor * (*(du_dx.value + HES_OFFSET_IJ));
    }
    for (const auto& du_du : hes.du_du) {
        values[ptr_map_xu_xu->access(off_x + du_du.row, off_x + du_du.col, block_count)] += factor * (*(du_du.value + HES_OFFSET_IJ));
    }
    for (const auto& dp_dx : hes.dp_dx) {
        values[ptr_map_p_xu->access(dp_dx.row, dp_dx.col, block_count)] += factor * (*(dp_dx.value + HES_OFFSET_IJ));
    }
    for (const auto& dp_du : hes.dp_du) {
        values[ptr_map_p_xu->access(dp_du.row, off_x + dp_du.col, block_count)] += factor * (*(dp_du.value + HES_OFFSET_IJ));
    }
    for (const auto& dp_dp : hes.dp_dp) {
        values[hes_h.access(dp_dp.row, dp_dp.col)] += factor * (*(dp_dp.value + HES_OFFSET_IJ));
    } 
}

void GDOP::updateHessianMR(FixedVector<double>& values, const HessianMR& hes, const double factor) {
    for (const auto& dx0_dx0 : hes.dx0_dx0) {
        values[hes_a.access(dx0_dx0.row, dx0_dx0.col)] += factor * (*(dx0_dx0.value));
    }
    for (const auto& dxf_dx0 : hes.dxf_dx0) {
        values[hes_c.access(dxf_dx0.row, dxf_dx0.col)] += factor * (*(dxf_dx0.value));
    }
    for (const auto& dxf_dxf : hes.dxf_dxf) {
        values[hes_d.access(dxf_dxf.row, dxf_dxf.col)] += factor * (*(dxf_dxf.value));
    }
    for (const auto& dp_dx0 : hes.dp_dx0) {
        values[hes_e.access(dp_dx0.row,  dp_dx0.col)]  += factor * (*(dp_dx0.value));
    }
    for (const auto& dp_dxf : hes.dp_dxf) {
        values[hes_g.access(dp_dxf.row,  dp_dxf.col)]  += factor * (*(dp_dxf.value));
    }
    for (const auto& dp_dp : hes.dp_dp) {
        values[hes_h.access(dp_dp.row,   dp_dp.col)]   += factor * (*(dp_dp.value));
    }
}

void GDOP::callback_evaluation() {
    problem->full->callbackEval(curr_x.raw() + off_x, curr_x.raw() + off_xu_total);
    problem->boundary->callbackEval(curr_x.raw(), curr_x.raw() + off_last_xu, curr_x.raw() + off_xu_total);
    evaluation_state.eval_f = true;
    evaluation_state.eval_g = true;
}

void GDOP::callback_jacobian() {
    if (!(evaluation_state.eval_f && evaluation_state.eval_g)) {
        callback_evaluation();
    }
    problem->full->callbackJac(curr_x.raw() + off_x, curr_x.raw() + off_xu_total);
    problem->boundary->callbackJac(curr_x.raw(), curr_x.raw() + off_last_xu, curr_x.raw() + off_xu_total);
    evaluation_state.grad_f = true;
    evaluation_state.jac_g = true;
}

void GDOP::callback_hessian() {
    if (!(evaluation_state.eval_f && evaluation_state.eval_g)) {
        callback_evaluation();
    }
    if (!(evaluation_state.grad_f && evaluation_state.jac_g)) {
        callback_jacobian();
    }
    problem->full->callbackHes(curr_x.raw() + off_x, curr_x.raw() + off_xu_total);
    problem->boundary->callbackHes(curr_x.raw(), curr_x.raw() + off_last_xu, curr_x.raw() + off_xu_total);
    evaluation_state.hes_lag = true;
}
