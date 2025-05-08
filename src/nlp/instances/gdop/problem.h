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
        eval_buffer = FixedVector<f64>(mesh.node_count * eval_size);
        jac_buffer = FixedVector<f64>(mesh.node_count * jac_size);
        hes_buffer = FixedVector<f64>(mesh.node_count * hes_size);
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

    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> hes_buffer; // can lead to severe memory consumption / order of a few GB, so maybe rethink this for large scale problems

    // fill eval_buffer, jac_buffer, hes_buffer
    virtual void callback_eval(const f64* xu_nlp, const f64* p) = 0;

    virtual void callback_jac(const f64* xu_nlp, const f64* p) = 0;

    virtual void callback_hes(const f64* xu_nlp, const f64* p) = 0;
    
    inline f64 get_eval_l(const int offset) {
        return *(lfg[0].eval + offset * eval_size);
    }

    inline f64 get_eval_f(const int f_index, const int offset) {
        return *(lfg[f_index_start + f_index].eval + offset * eval_size);
    }

    inline f64 get_eval_g(const int g_index, const int offset) {
        return *(lfg[g_index_start + g_index].eval + offset * eval_size);
    }
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
        eval_buffer = FixedVector<f64>(r_index_end);
        jac_buffer = FixedVector<f64>(jac_buffer_size);
        hes_buffer = FixedVector<f64>(hes_buffer_size);
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

    FixedVector<f64> eval_buffer;
    FixedVector<f64> jac_buffer;
    FixedVector<f64> hes_buffer;

    virtual void callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

    virtual void callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

    virtual void callback_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p) = 0;

    inline f64 get_eval_m() {
        return *(mr[0].eval);
    }

    inline f64 get_eval_r(const int r_index) {
        return *(mr[r_index_start + r_index].eval);
    }
};

class Problem {
public:
    Problem(FullSweep& full, BoundarySweep& boundary, FixedVector<Bounds>&& x_bounds,
               FixedVector<Bounds>&& u_bounds, FixedVector<Bounds>&& p_bounds, FixedVector<std::optional<f64>>&& x0_fixed,
               FixedVector<std::optional<f64>>&& xf_fixed)
    : full(full), boundary(boundary), x_bounds(std::move(x_bounds)), u_bounds(std::move(u_bounds)),
      p_bounds(std::move(p_bounds)), x0_fixed(std::move(x0_fixed)), xf_fixed(std::move(xf_fixed)),
      x_size(this->full.x_size), u_size(this->full.u_size), p_size(this->full.p_size) {
    };
    
    FullSweep& full;
    BoundarySweep& boundary;

    FixedVector<Bounds> x_bounds;
    FixedVector<Bounds> u_bounds;
    FixedVector<Bounds> p_bounds;

    FixedVector<std::optional<f64>> x0_fixed; // set value if a state has a fixed initial value, remove the constraint from r()!
    FixedVector<std::optional<f64>> xf_fixed; // set value if a state has a fixed final value, remove the constraint from r()!

    int x_size;
    int u_size;
    int p_size;
};

#endif // OPT_GDOP_PROBLEM_H
