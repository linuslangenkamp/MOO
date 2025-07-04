#ifndef OPT_GDOP_PROBLEM_H
#define OPT_GDOP_PROBLEM_H

#include <optional>
#include <memory>
#include <vector>

#include <base/nlp_structs.h>
#include <base/mesh.h>


class FullSweep {
public:
    // Evaluation of L(), f(), g() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T + first and second derivatives
    // Idea: Call Fullsweep.setValues(), Fullsweep.callback_eval(), callback_jac(), callHess() -> Just iterate over COO
    // function evals / diffs are on callback interfaced side

    FullSweep(FixedVector<FunctionLFG>&& lfg_in, std::unique_ptr<AugmentedHessianLFG> aug_hes, std::unique_ptr<AugmentedParameterHessian> aug_pp_hes,
     Collocation& collocation, Mesh& mesh, FixedVector<Bounds>&& g_bounds, bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size)
    : lfg(std::move(lfg_in)), aug_hes(std::move(aug_hes)), aug_pp_hes(std::move(aug_pp_hes)), collocation(collocation), mesh(mesh), g_bounds(std::move(g_bounds)),
      has_lagrange(has_lagrange), f_index_start((int) has_lagrange), f_index_end(f_index_start + f_size), g_index_start(f_index_end),
      g_index_end(g_index_start + g_size), f_size(f_size), g_size(g_size), fg_size(f_size + g_size), x_size(x_size), u_size(u_size),
      p_size(p_size), eval_size(lfg.size()), aug_hes_size(this->aug_hes->nnz()) {
        for (const auto& func : lfg) {
            jac_size += func.jac.nnz();
        }
        eval_buffer = FixedVector<f64>(mesh.node_count * eval_size);
        jac_buffer = FixedVector<f64>(mesh.node_count * jac_size);
        aug_hes_buffer = FixedVector<f64>(mesh.node_count * aug_hes_size);
    };

    virtual ~FullSweep() = default;

    FixedVector<FunctionLFG> lfg;
    std::unique_ptr<AugmentedHessianLFG> aug_hes;
    std::unique_ptr<AugmentedParameterHessian> aug_pp_hes;
    Collocation& collocation;
    Mesh& mesh;

    FixedVector<Bounds> g_bounds;

    bool has_lagrange = false;
    int f_index_start = 0;
    int f_index_end = 0;
    int g_index_start = 0;
    int g_index_end = 0;
    int f_size = 0;
    int g_size = 0;
    int fg_size = 0;

    int x_size = 0;
    int u_size = 0;
    int p_size = 0;

    // sizes of 1 buffer chuck
    int eval_size = 0;
    int jac_size = 0;
    int aug_hes_size = 0;

    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> aug_hes_buffer;

    /* TODO: add this buffer for parallel parameters, make this threaded; #threads of these buffers; sum them at the end */
    FixedVector<f64> aug_pp_hes_buffer; // make it like list<FixedVector<f64>>, each thread sum to own buffer (just size p * p each)

    // TODO: add good explanation of the interface

    // fill eval_buffer accoring to sparsity structure
    virtual void callback_eval(const f64* xu_nlp, const f64* p) = 0;

    // fill jac_buffer accoring to sparsity structure
    virtual void callback_jac(const f64* xu_nlp, const f64* p) = 0;

    /* fill aug_hes_buffer and aug_pp_hes_buffer accoring to sparsity structure
     * aug_hes_buffer is $\lambda^T * \nabla² (f, g) + lfactor * \nabla² L$ all except w.r.t. pp
     * aug_pp_hes_buffer is $\lambda^T * \nabla²_{pp} (f, g) + lfactor * \nabla²_{pp} L$
     * lambdas are exact multipliers (no transform needed) to each block [f, g]_{ij}
     * lagrange_factors are exact factor for lagrange terms in interval i, nodes j */
    virtual void callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, f64* lambda) = 0;

    void print_jacobian_sparsity_pattern() {
        std::cout << "\n=== LFG Jacobian Sparsity ===\n================================\n";
        for (size_t i = 0; i < lfg.size(); ++i) {
            std::cout << "FunctionLFG[" << i << "] - ";
            if (has_lagrange && i == 0) {
                std::cout << "L - Lagrange term:\n";
            } else if (i >= static_cast<size_t>(f_index_start) && i < static_cast<size_t>(f_index_end)) {
                std::cout << "f[" << (i - f_index_start) << "] - Dynamic Equation:\n";
            } else if (i >= static_cast<size_t>(g_index_start) && i < static_cast<size_t>(g_index_end)) {
                std::cout << "g[" << (i - g_index_start) << "] - Path Constraint:\n";
            } 

            std::cout << "  dx sparsity pattern:\n";
            for (const auto& entry : lfg[i].jac.dx) {
                std::cout << "    x_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
            }

            std::cout << "  du sparsity pattern:\n";
            for (const auto& entry : lfg[i].jac.du) {
                std::cout << "    u_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
            }

            std::cout << "  dp sparsity pattern:\n";
            for (const auto& entry : lfg[i].jac.dp) {
                std::cout << "    p_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
            }

            std::cout << "-----------------------------\n";
        }
        std::cout << "================================\n";
    }
};

class BoundarySweep {
public:
    BoundarySweep(FixedVector<FunctionMR>&& mr_in, std::unique_ptr<AugmentedHessianMR> aug_hes, Mesh& mesh,
                  FixedVector<Bounds>&& r_bounds, bool has_mayer, int r_size, int x_size, int p_size)
    : mr(std::move(mr_in)), aug_hes(std::move(aug_hes)), mesh(mesh), r_bounds(std::move(r_bounds)), has_mayer(has_mayer),
      r_index_start((int) has_mayer), r_index_end(r_index_start + r_size), r_size(r_size), x_size(x_size), p_size(p_size) {
        int jac_buffer_size = 0;
        for (const auto& func : mr) {
            jac_buffer_size += func.jac.nnz();
        }
        eval_buffer = FixedVector<f64>(r_index_end);
        jac_buffer = FixedVector<f64>(jac_buffer_size);
        aug_hes_buffer = FixedVector<f64>(this->aug_hes->nnz());
    };

    virtual ~BoundarySweep() = default;

    // M, R :: assert x0 Size == xf Size
    FixedVector<FunctionMR> mr;
    std::unique_ptr<AugmentedHessianMR> aug_hes;
    Mesh& mesh;

    FixedVector<Bounds> r_bounds;

    bool has_mayer;
    int r_index_start;
    int r_index_end;
    int r_size; // assert; check with fixed initial states, these should not be contained here!

    int x_size;
    int p_size;

    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> aug_hes_buffer;

    virtual void callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

    virtual void callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

   /* lambdas are exact multipliers (no transform needed) to [r]
    * mayer_factor is eact multiplier (no transform needed) of M */
    virtual void callback_aug_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p, const f64 mayer_factor, f64* lambda) = 0;

    void print_jacobian_sparsity_pattern() {
        std::cout << "\n=== MR Jacobian Sparsity ===\n================================\n";
        for (size_t i = 0; i < mr.size(); ++i) {
            std::cout << "FunctionMR[" << i << "] - ";
            if (has_mayer && i == 0) {
                std::cout << "M - Mayer term:\n";
            } else if (i >= static_cast<size_t>(r_index_start) && i < static_cast<size_t>(r_index_end)) {
                std::cout << "r[" << (i - r_index_start) << "] - Boundary Constraint:\n";
            } else {
                std::cout << "(unknown):\n";
            }

            std::cout << "  dx0 sparsity pattern:\n";
            for (const auto& entry : mr[i].jac.dx0) {
                std::cout << "    x0_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
            }

            std::cout << "  dxf sparsity pattern:\n";
            for (const auto& entry : mr[i].jac.dxf) {
                std::cout << "    xf_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
            }

            std::cout << "  dp sparsity pattern:\n";
            for (const auto& entry : mr[i].jac.dp) {
                std::cout << "    p_idx = " << entry.col << " (jac_buf_index = " << entry.buf_index << ")\n";
            }
            std::cout << "-----------------------------\n";
        }
        std::cout << "================================\n";
    }
};

class Problem {
public:
    Problem(std::unique_ptr<FullSweep>&& full, std::unique_ptr<BoundarySweep>&& boundary, Mesh& mesh, FixedVector<Bounds>&& x_bounds,
               FixedVector<Bounds>&& u_bounds, FixedVector<Bounds>&& p_bounds, FixedVector<std::optional<f64>>&& x0_fixed,
               FixedVector<std::optional<f64>>&& xf_fixed)
    : full(std::move(full)), boundary(std::move(boundary)), mesh(mesh), x_bounds(std::move(x_bounds)), u_bounds(std::move(u_bounds)),
      p_bounds(std::move(p_bounds)), x0_fixed(std::move(x0_fixed)), xf_fixed(std::move(xf_fixed)),
      x_size(this->full->x_size), u_size(this->full->u_size), p_size(this->full->p_size) {
    };

    std::unique_ptr<FullSweep> full;
    std::unique_ptr<BoundarySweep> boundary;

    Mesh& mesh;

    FixedVector<Bounds> x_bounds;
    FixedVector<Bounds> u_bounds;
    FixedVector<Bounds> p_bounds;

    FixedVector<std::optional<f64>> x0_fixed; // set value if a state has a fixed initial value, remove the constraint from r()!
    FixedVector<std::optional<f64>> xf_fixed; // set value if a state has a fixed final value, remove the constraint from r()!

    int x_size;
    int u_size;
    int p_size;

    // FIXME: TODO: get rid of int where possible: below could actually overflow with decent hardware!
    //              for now it might be sufficient to just static cast to size_t for the calculations, since
    //              in nearly all other calculations with indices these are not that large

    /* note entry != index in lfg, but rather full->lfg[*].buf_index */
    inline f64 lfg_eval(int entry, int interval_i, int node_j) {
        return full->eval_buffer[entry + full->eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_L(int interval_i, int node_j) {
        assert(full->has_lagrange);
        return full->eval_buffer[full->lfg[0].buf_index + full->eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_f(int f_index, int interval_i, int node_j) {
        return full->eval_buffer[full->lfg[full->f_index_start + f_index].buf_index + full->eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_g(int g_index, int interval_i, int node_j) {
        return full->eval_buffer[full->lfg[full->g_index_start + g_index].buf_index + full->eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_jac(int entry, int interval_i, int node_j) {
        return full->jac_buffer[entry + full->jac_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_aug_hes(int entry, int interval_i, int node_j) {
        return full->aug_hes_buffer[entry + full->aug_hes_size * mesh.acc_nodes[interval_i][node_j]];
    }
    /* TODO: add and make threaded */
    inline f64 lfg_aug_pp_hes(int entry) {
        return full->aug_pp_hes_buffer[entry];
    }

    /* note entry != index in mr, but rather boundary->mr[*].buf_index */
    inline f64 mr_eval(int entry) {
        return full->eval_buffer[entry];
    }

    inline f64 mr_eval_M() {
        assert(boundary->has_mayer);
        return boundary->eval_buffer[boundary->mr[0].buf_index];
    }

    inline f64 mr_eval_r(int r_index) {
        return boundary->eval_buffer[boundary->mr[boundary->r_index_start + r_index].buf_index];
    }

    inline f64 mr_jac(int entry) {
        return boundary->jac_buffer[entry];
    }

    inline f64 mr_aug_hes(int entry) {
        return boundary->aug_hes_buffer[entry];
    }
};

#endif // OPT_GDOP_PROBLEM_H
