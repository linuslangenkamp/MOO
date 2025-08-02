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
    : x_size(static_cast<int>(x_bounds.size())),
      u_size(static_cast<int>(u_bounds.size())),
      p_size(static_cast<int>(p_bounds.size())),
      xu_size(x_size + u_size),
      has_mayer(has_mayer),
      has_lagrange(has_lagrange),
      f_size(x_size),
      g_size(static_cast<int>(g_bounds.size())),
      r_size(static_cast<int>(r_bounds.size())),
      fg_size(f_size + g_size),
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

// ================== Continuous, Dynamic - Lfg Block ==================

struct BlockLFG {
    std::unique_ptr<FunctionLFG> L; // Lagrange Term
    FixedVector<FunctionLFG> f;     // Dynamic Equations
    FixedVector<FunctionLFG> g;     // Path Constraints

    BlockLFG(bool lagrange_exists,
             int size_f,
             int size_g)
    : L(lagrange_exists ? std::make_unique<FunctionLFG>() : nullptr),
      f(FixedVector<FunctionLFG>(size_f)),
      g(FixedVector<FunctionLFG>(size_g)) {}


    inline int size() {
        return (L ? 1 : 0) + f.int_size() + g.int_size();
    };

    int compute_jac_nnz();
};

struct FullSweepBuffers {
    // sizes of 1 buffer chuck
    const int eval_size = 0;
    const int jac_size = 0;
    const int aug_hes_size = 0;

    FixedVector<f64> eval;
    FixedVector<f64> jac;
    FixedVector<f64> aug_hes;

    /* TODO: add this buffer for parallel parameters, make this threaded; #threads of these buffers; sum them at the end */
    FixedVector<f64> aug_pp_hes; // make it like array<FixedVector<f64>>, each thread sum to own buffer (just size p * p each)

    FullSweepBuffers(BlockLFG& lfg,
                     AugmentedHessianLFG& aug_hes,
                     const ProblemConstants& pc) :
        eval_size(lfg.size()),
        jac_size(lfg.compute_jac_nnz()),
        aug_hes_size(aug_hes.nnz())
    {
        resize(pc.mesh);
    }

    void resize(const Mesh& mesh);
};

class FullSweep {
    friend class Problem;

public:
    FullSweep(BlockLFG&& lfg_in,
              std::unique_ptr<AugmentedHessianLFG> aug_hes_in,
              std::unique_ptr<AugmentedParameterHessian> aug_pp_hes_in,
              const ProblemConstants& pc_in)
    : lfg(std::move(lfg_in)),
      aug_hes(std::move(aug_hes_in)),
      aug_pp_hes(std::move(aug_pp_hes_in)),
      pc(pc_in),
      buffers(lfg, *aug_hes, pc) {}

    virtual ~FullSweep() = default;

    // eval + jacobian structure (includes sparsity + mapping to buffer indices - location to write to)
    BlockLFG lfg;

    // augmented Hessian structures (excluding parameters / parameters only)
    std::unique_ptr<AugmentedHessianLFG> aug_hes;
    std::unique_ptr<AugmentedParameterHessian> aug_pp_hes;

    // stores all the relevant constants such as dimensions, bounds, offsets, mesh and collocation
    const ProblemConstants& pc;

    // fill eval_buffer accoring to sparsity structure
    virtual void callback_eval(const f64* xu_nlp, const f64* p) = 0;

    // fill jac_buffer accoring to sparsity structure
    virtual void callback_jac(const f64* xu_nlp, const f64* p) = 0;

    /* fill aug_hes_buffer and aug_pp_hes_buffer accoring to sparsity structure
     * aug_hes_buffer is $\lambda^T * \nabla² (f, g) + lfactor * \nabla² L$ all except w.r.t. pp
     * aug_pp_hes_buffer is $\lambda^T * \nabla²_{pp} (f, g) + lfactor * \nabla²_{pp} L$
     * lambdas are exact multipliers (no transform needed) to each block [f, g]_{ij}
     * lagrange_factors are exact factor for lagrange terms in interval i, nodes j */
    virtual void callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, const f64* lambda) = 0;

    inline const f64* get_xu_ij(const f64* xu_nlp, int i, int j) {
        return xu_nlp + pc.xu_size * pc.mesh.acc_nodes[i][j];
    }

    inline const f64* get_lambda_ij(const f64* lambda, int i, int j) {
        return lambda + pc.fg_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_eval_buffer(int i, int j) {
        return buffers.eval.raw() + buffers.eval_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_jac_buffer(int i, int j) {
        return buffers.jac.raw() + buffers.jac_size * pc.mesh.acc_nodes[i][j];
    }

    inline f64* get_aug_hes_buffer(int i, int j) {
        return buffers.aug_hes.raw() + buffers.aug_hes_size * pc.mesh.acc_nodes[i][j];
    }

    void print_jacobian_sparsity_pattern();

private:
    // buffers to write to in callbacks
    FullSweepBuffers buffers;
};

// ================== Boundary - Mr Block ==================

struct BlockMR {
    std::unique_ptr<FunctionMR> M; // Mayer Term
    FixedVector<FunctionMR> r;     // Boundary Constraints

    BlockMR(bool mayer_exists,
            int size_r)
    : M(mayer_exists ? std::make_unique<FunctionMR>() : nullptr),
      r(FixedVector<FunctionMR>(size_r)) {}

    inline int size() {
        return (M ? 1 : 0) + r.int_size();
    };

    int compute_jac_nnz();
};

struct BoundarySweepBuffers {
    FixedVector<f64> eval;
    FixedVector<f64> jac;
    FixedVector<f64> aug_hes;

    BoundarySweepBuffers(BlockMR& mr,
                         AugmentedHessianMR& aug_hes,
                         const ProblemConstants& pc)
    : eval(FixedVector<f64>(mr.size())),
      jac(FixedVector<f64>(mr.compute_jac_nnz())),
      aug_hes(FixedVector<f64>(aug_hes.nnz())) {}
};

class BoundarySweep {
    friend class Problem;

public:
    BoundarySweep(BlockMR&& mr_in,
                  std::unique_ptr<AugmentedHessianMR> aug_hes_in,
                  const ProblemConstants& pc_in)
    : mr(std::move(mr_in)),
      aug_hes(std::move(aug_hes_in)),
      pc(pc_in),
      buffers(mr, *aug_hes, pc) {}

    virtual ~BoundarySweep() = default;

    // eval + jacobian structure (includes sparsity + mapping to buffer indices - location to write to)
    BlockMR mr;

     // augmented Hessian structure
    std::unique_ptr<AugmentedHessianMR> aug_hes;

    // stores all the relevant constants such as dimensions, bounds, offsets, mesh and collocation
    const ProblemConstants& pc;

    virtual void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) = 0;

    virtual void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) = 0;

   /* lambdas are exact multipliers (no transform needed) to [r]
    * mayer_factor is eact multiplier (no transform needed) of M */
    virtual void callback_aug_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, const f64* lambda) = 0;

    inline f64* get_eval_buffer() {
        return buffers.eval.raw();
    }

    inline f64* get_jac_buffer() {
        return buffers.jac.raw();
    }

    inline f64* get_aug_hes_buffer() {
        return buffers.aug_hes.raw();
    }

    inline void fill_zero_aug_hes_buffer() {
        buffers.aug_hes.fill_zero();
    }

    void print_jacobian_sparsity_pattern();

private:
    // buffers to write to in callbacks
    BoundarySweepBuffers buffers;
};

class Problem {
public:
    Problem(std::unique_ptr<FullSweep>&& full, std::unique_ptr<BoundarySweep>&& boundary, std::unique_ptr<const ProblemConstants>&& pc)
    : full(std::move(full)), boundary(std::move(boundary)), pc(std::move(pc)) {
    };

    std::unique_ptr<FullSweep> full;
    std::unique_ptr<BoundarySweep> boundary;
    std::unique_ptr<const ProblemConstants> pc;

    // FIXME: TODO: get rid of int where possible: below could actually overflow with decent hardware!
    //              for now it might be sufficient to just static cast to size_t for the calculations, since
    //              in nearly all other calculations with indices these are not that large

    inline f64 lfg_eval_L(int interval_i, int node_j) {
        assert(full->pc.has_lagrange && full->lfg.L);
        return full->buffers.eval[full->lfg.L->buf_index + full->buffers.eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_f(int f_index, int interval_i, int node_j) {
        return full->buffers.eval[full->lfg.f[f_index].buf_index + full->buffers.eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_eval_g(int g_index, int interval_i, int node_j) {
        return full->buffers.eval[full->lfg.g[g_index].buf_index + full->buffers.eval_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_jac(int entry, int interval_i, int node_j) {
        return full->buffers.jac[entry + full->buffers.jac_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    inline f64 lfg_aug_hes(int entry, int interval_i, int node_j) {
        return full->buffers.aug_hes[entry + full->buffers.aug_hes_size * pc->mesh.acc_nodes[interval_i][node_j]];
    }

    /* TODO: add and make threaded */
    inline f64 lfg_aug_pp_hes(int entry) {
        return full->buffers.aug_pp_hes[entry];
    }

    inline f64 mr_eval_M() {
        assert(boundary->pc.has_mayer && boundary->mr.M);
        return boundary->buffers.eval[boundary->mr.M->buf_index];
    }

    inline f64 mr_eval_r(int r_index) {
        return boundary->buffers.eval[boundary->mr.r[r_index].buf_index];
    }

    inline f64 mr_jac(int entry) {
        return boundary->buffers.jac[entry];
    }

    inline f64 mr_aug_hes(int entry) {
        return boundary->buffers.aug_hes[entry];
    }

    inline void resize_buffers(const Mesh& mesh) {
        full->buffers.resize(mesh);
    }
};

} // namespace GDOP

#endif // OPT_GDOP_PROBLEM_H
