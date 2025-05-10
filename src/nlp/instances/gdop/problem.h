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

    FullSweep(FixedVector<FunctionLFG>&& lfg_in, Mesh& mesh, FixedVector<Bounds>&& g_bounds, bool has_lagrange,
              int f_size, int g_size, int x_size, int u_size, int p_size)
    : lfg(std::move(lfg_in)), mesh(mesh), g_bounds(std::move(g_bounds)), has_lagrange(has_lagrange), f_index_start((int) has_lagrange),
      f_index_end(f_index_start + f_size), g_index_start(f_index_end), g_index_end(g_index_start + g_size), f_size(f_size),
      g_size(g_size), fg_size(f_size + g_size), x_size(x_size), u_size(u_size), p_size(p_size), eval_size(lfg.size()) {
        for (const auto& func : lfg) {
            jac_size += func.jac.nnz();
            hes_size += func.hes.nnz();
        }
        eval_buffer = FixedVector<F64>(mesh.node_count * eval_size);
        jac_buffer = FixedVector<F64>(mesh.node_count * jac_size);
        hes_buffer = FixedVector<F64>(mesh.node_count * hes_size);
    };

    FixedVector<FunctionLFG> lfg;
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
    int hes_size = 0;

    FixedVector<F64> eval_buffer;
    FixedVector<F64> jac_buffer;
    FixedVector<F64> hes_buffer; // can lead to severe memory consumption / order of a few GB, so maybe rethink this for large scale problems

    // fill eval_buffer, jac_buffer, hes_buffer
    virtual void callback_eval(const F64* xu_nlp, const F64* p) = 0;

    virtual void callback_jac(const F64* xu_nlp, const F64* p) = 0;

    virtual void callback_hes(const F64* xu_nlp, const F64* p) = 0;
};

class BoundarySweep {
public:

    BoundarySweep(FixedVector<FunctionMR>&& mr_in, Mesh& mesh, FixedVector<Bounds>&& r_bounds, bool has_mayer,
                  int r_size, int x_size, int p_size)
    : mr(std::move(mr_in)), mesh(mesh), r_bounds(std::move(r_bounds)), has_mayer(has_mayer), r_index_start((int) has_mayer),
      r_index_end(r_index_start + r_size), r_size(r_size), x_size(x_size), p_size(p_size) {
        int jac_buffer_size = 0;
        int hes_buffer_size = 0;
        for (const auto& func : mr) {
            jac_buffer_size += func.jac.nnz();
            hes_buffer_size += func.hes.nnz();
        }
        eval_buffer = FixedVector<F64>(r_index_end);
        jac_buffer = FixedVector<F64>(jac_buffer_size);
        hes_buffer = FixedVector<F64>(hes_buffer_size);
    };

    // M, R :: assert x0 Size == xf Size
    FixedVector<FunctionMR> mr;
    Mesh& mesh;

    FixedVector<Bounds> r_bounds;

    bool has_mayer;
    int r_index_start;
    int r_index_end;
    int r_size; // assert; check with fixed initial states, these should not be contained here!

    int x_size;
    int p_size;

    FixedVector<F64> eval_buffer;
    FixedVector<F64> jac_buffer;
    FixedVector<F64> hes_buffer;

    virtual void callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) = 0;

    virtual void callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) = 0;

    virtual void callback_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p) = 0;
};

class Problem {
public:
    Problem(FullSweep& full, BoundarySweep& boundary, Mesh& mesh, FixedVector<Bounds>&& x_bounds,
               FixedVector<Bounds>&& u_bounds, FixedVector<Bounds>&& p_bounds, FixedVector<std::optional<F64>>&& x0_fixed,
               FixedVector<std::optional<F64>>&& xf_fixed)
    : full(full), boundary(boundary), mesh(mesh), x_bounds(std::move(x_bounds)), u_bounds(std::move(u_bounds)),
      p_bounds(std::move(p_bounds)), x0_fixed(std::move(x0_fixed)), xf_fixed(std::move(xf_fixed)),
      x_size(this->full.x_size), u_size(this->full.u_size), p_size(this->full.p_size) {
    };
    
    FullSweep& full;
    BoundarySweep& boundary;

    Mesh& mesh;

    FixedVector<Bounds> x_bounds;
    FixedVector<Bounds> u_bounds;
    FixedVector<Bounds> p_bounds;

    FixedVector<std::optional<F64>> x0_fixed; // set value if a state has a fixed initial value, remove the constraint from r()!
    FixedVector<std::optional<F64>> xf_fixed; // set value if a state has a fixed final value, remove the constraint from r()!

    int x_size;
    int u_size;
    int p_size;

   // TODO: get rid of int where possible: below could actually overflow with decent hardware!

    /* note entry != index in lfg, but rather full.lfg[*].buf_index */
    inline F64 lfg_eval(int entry, int interval_i, int node_j) {
        return full.eval_buffer[entry + full.eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline F64 lfg_eval_L(int interval_i, int node_j) {
        assert(full.has_lagrange);
        return full.eval_buffer[full.lfg[0].buf_index + full.eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline F64 lfg_eval_f(int f_index, int interval_i, int node_j) {
        return full.eval_buffer[full.lfg[full.f_index_start + f_index].buf_index + full.eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline F64 lfg_eval_g(int g_index, int interval_i, int node_j) {
        return full.eval_buffer[full.lfg[full.g_index_start + g_index].buf_index + full.eval_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline F64 lfg_jac(int entry, int interval_i, int node_j) {
        return full.jac_buffer[entry + full.jac_size * mesh.acc_nodes[interval_i][node_j]];
    }

    inline F64 lfg_hes(int entry, int interval_i, int node_j) {
        return full.hes_buffer[entry + full.hes_size * mesh.acc_nodes[interval_i][node_j]];
    }

    /* note entry != index in mr, but rather boundary.mr[*].buf_index */
    inline F64 mr_eval(int entry) {
        return full.eval_buffer[entry];
    }

    inline F64 mr_eval_M() {
        assert(boundary.has_mayer);
        return boundary.eval_buffer[boundary.mr[0].buf_index];
    }

    inline F64 mr_eval_r(int r_index) {
        return boundary.eval_buffer[boundary.mr[boundary.r_index_start + r_index].buf_index];
    }

    inline F64 mr_jac(int entry) {
        return boundary.jac_buffer[entry];
    }
    
    inline F64 mr_hes(int entry) {
        return boundary.hes_buffer[entry];
    }
};

#endif // OPT_GDOP_PROBLEM_H
