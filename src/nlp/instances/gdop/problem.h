#ifndef OPT_GDOP_PROBLEM_H
#define OPT_GDOP_PROBLEM_H

#include <optional>
#include <memory>
#include <vector>

#include <base/nlp_structs.h>
#include <base/mesh.h>

namespace GDOP {


struct ProblemConstants {
    // sizes
    const int x_size;
    const int u_size;
    const int p_size;
    const int xu_size;

    // objective structure
    const bool has_mayer;
    const bool has_lagrange;

    // function dims
    const int f_size;
    const int g_size;
    const int r_size;
    const int fg_size;

    // index ranges in buffers
    const int f_index_start;
    const int f_index_end;
    const int g_index_start;
    const int g_index_end;
    const int r_index_start;
    const int r_index_end;

    // variable bounds
    const FixedVector<Bounds> x_bounds;
    const FixedVector<Bounds> u_bounds;
    const FixedVector<Bounds> p_bounds;

    // fixed initial and final states
    const FixedVector<std::optional<f64>> x0_fixed;
    const FixedVector<std::optional<f64>> xf_fixed;

    // bounds
    const FixedVector<Bounds> r_bounds;
    const FixedVector<Bounds> g_bounds;

    // Mesh and Collocation references
    const Mesh& mesh;
    const Collocation& collocation;

    ProblemConstants(
        int x_size,
        int u_size,
        int p_size,
        bool has_mayer,
        bool has_lagrange,
        FixedVector<Bounds>&& x_bounds,
        FixedVector<Bounds>&& u_bounds,
        FixedVector<Bounds>&& p_bounds,
        FixedVector<std::optional<f64>> x0_fixed,
        FixedVector<std::optional<f64>> xf_fixed,
        FixedVector<Bounds>&& r_bounds,
        FixedVector<Bounds>&& g_bounds,
        const Mesh& mesh,
        const Collocation& collocation)
    : x_size(x_size),
      u_size(u_size),
      p_size(p_size),
      xu_size(x_size + u_size),
      has_mayer(has_mayer),
      has_lagrange(has_lagrange),
      f_size(x_size),
      g_size(static_cast<int>(g_bounds.size())),
      r_size(r_bounds.size()),
      fg_size(f_size + g_size),
      f_index_start(static_cast<int>(has_lagrange)),
      f_index_end(f_index_start + f_size),
      g_index_start(f_index_end),
      g_index_end(g_index_start + g_size),
      r_index_start(has_mayer ? 1 : 0),
      r_index_end(r_index_start + r_size),
      x_bounds(std::move(x_bounds)),
      u_bounds(std::move(u_bounds)),
      p_bounds(std::move(p_bounds)),
      x0_fixed(std::move(x0_fixed)),
      xf_fixed(std::move(xf_fixed)),
      r_bounds(std::move(r_bounds)),
      g_bounds(std::move(g_bounds)),
      mesh(mesh),
      collocation(collocation)
    {}
};

class FullSweep {
public:
    // Evaluation of L(), f(), g() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T + first and second derivatives
    // Idea: Call Fullsweep.setValues(), Fullsweep.callback_eval(), callback_jac(), callHess() -> Just iterate over COO
    // function evals / diffs are on callback interfaced side

    FullSweep(FixedVector<FunctionLFG>&& lfg_in,
              std::unique_ptr<AugmentedHessianLFG> aug_hes,
              std::unique_ptr<AugmentedParameterHessian> aug_pp_hes,
              const ProblemConstants& pc)
    : lfg(std::move(lfg_in)),
      aug_hes(std::move(aug_hes)),
      aug_pp_hes(std::move(aug_pp_hes)),
      pc(pc),
      eval_size(lfg.size()),
      jac_size(compute_jac_nnz_vec_fn_lfg(lfg)),
      aug_hes_size(this->aug_hes->nnz())
      {
        resize_buffers();
    };

    virtual ~FullSweep() = default;
    // TODO: add good explanation of the interface

    FixedVector<FunctionLFG> lfg;
    std::unique_ptr<AugmentedHessianLFG> aug_hes;
    std::unique_ptr<AugmentedParameterHessian> aug_pp_hes;
    const ProblemConstants& pc;

    // sizes of 1 buffer chuck
    int eval_size = 0;
    int jac_size = 0;
    int aug_hes_size = 0;

    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> aug_hes_buffer;

    /* TODO: add this buffer for parallel parameters, make this threaded; #threads of these buffers; sum them at the end */
    FixedVector<f64> aug_pp_hes_buffer; // make it like list<FixedVector<f64>>, each thread sum to own buffer (just size p * p each)

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

    inline const f64* get_xu_ij(const f64* xu_nlp, int i, int j) {
        return xu_nlp + pc.xu_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_lambda_ij(f64* lambda, int i, int j) {
        return lambda + pc.fg_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_eval_buffer(int i, int j) {
        return eval_buffer.raw() + eval_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_jac_buffer(int i, int j) {
        return jac_buffer.raw() + jac_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_aug_hes_buffer(int i, int j) {
        return aug_hes_buffer.raw() + aug_hes_size * pc.mesh.acc_nodes[i][j];
    }

    void resize_buffers();
    void print_jacobian_sparsity_pattern();
};

class BoundarySweep {
public:
    BoundarySweep(FixedVector<FunctionMR>&& mr_in,
                  std::unique_ptr<AugmentedHessianMR> aug_hes,
                  const ProblemConstants& pc)
    : mr(std::move(mr_in)),
      aug_hes(std::move(aug_hes)),
      pc(pc)
    {
        resize_buffers();
    };

    virtual ~BoundarySweep() = default;

    FixedVector<FunctionMR> mr;
    std::unique_ptr<AugmentedHessianMR> aug_hes;
    const ProblemConstants& pc;

    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> aug_hes_buffer;

    virtual void callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

    virtual void callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

   /* lambdas are exact multipliers (no transform needed) to [r]
    * mayer_factor is eact multiplier (no transform needed) of M */
    virtual void callback_aug_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p, const f64 mayer_factor, f64* lambda) = 0;

    void resize_buffers();
    void print_jacobian_sparsity_pattern();
};

class Problem {
public:
    Problem(std::unique_ptr<FullSweep>&& full, std::unique_ptr<BoundarySweep>&& boundary, std::unique_ptr<const ProblemConstants>&& pc)
    : full(std::move(full)), boundary(std::move(boundary)), pc(std::move(pc)) {
    };

    std::unique_ptr<FullSweep> full;
    std::unique_ptr<BoundarySweep> boundary;
    std::unique_ptr<const ProblemConstants> pc;

    void resize_buffers();

    // FIXME: TODO: get rid of int where possible: below could actually overflow with decent hardware!
    //              for now it might be sufficient to just static cast to size_t for the calculations, since
    //              in nearly all other calculations with indices these are not that large

    /* note entry != index in lfg, but rather full->lfg[*].buf_index */
    inline f64 lfg_eval(int entry, int interval_i, int node_j) {
        return full->eval_buffer[entry + full->eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_L(int interval_i, int node_j) {
        assert(full->pc.has_lagrange);
        return full->eval_buffer[full->lfg[0].buf_index + full->eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_f(int f_index, int interval_i, int node_j) {
        return full->eval_buffer[full->lfg[pc->f_index_start + f_index].buf_index + full->eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_g(int g_index, int interval_i, int node_j) {
        return full->eval_buffer[full->lfg[pc->g_index_start + g_index].buf_index + full->eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_jac(int entry, int interval_i, int node_j) {
        return full->jac_buffer[entry + full->jac_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_aug_hes(int entry, int interval_i, int node_j) {
        return full->aug_hes_buffer[entry + full->aug_hes_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    /* TODO: add and make threaded */
    inline f64 lfg_aug_pp_hes(int entry) {
        return full->aug_pp_hes_buffer[entry];
    }

    inline f64 mr_eval_M() {
        assert(boundary->pc.has_mayer);
        return boundary->eval_buffer[boundary->mr[0].buf_index];
    }

    inline f64 mr_eval_r(int r_index) {
        return boundary->eval_buffer[boundary->mr[pc->r_index_start + r_index].buf_index];
    }

    inline f64 mr_jac(int entry) {
        return boundary->jac_buffer[entry];
    }

    inline f64 mr_aug_hes(int entry) {
        return boundary->aug_hes_buffer[entry];
    }
};

} // namespace GDOP

#endif // OPT_GDOP_PROBLEM_H
