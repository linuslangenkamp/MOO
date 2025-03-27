#ifndef OPT_GDOP_PROBLEM_H
#define OPT_GDOP_PROBLEM_H

#include <optional>
#include <memory>
#include <vector>

#include "../base/nlp_structs.h"
#include "../base/mesh.h"


class FullSweep {
public:
    // Evaluation of L(), f(), g() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T + first and second derivatives
    // Idea: Call Fullsweep.setValues(), Fullsweep.callbackEval(), callbackJac(), callHess() -> Just iterate over COO
    // function evals / diffs are on callback interfaced side

    FullSweep(FixedVector<FunctionLFG>&& lfg_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds, bool has_lagrange,
              int f_size, int g_size, int x_size, int u_size, int p_size)
    : lfg(std::move(lfg_in)), mesh(mesh), f_size(f_size), g_size(g_size), fg_size(f_size + g_size), has_lagrange(has_lagrange), f_index_start((int) has_lagrange),
      f_index_end(f_index_start + f_size), g_index_start(f_index_end), g_index_end(g_index_start + g_size), x_size(x_size),
      u_size(u_size), p_size(p_size), g_bounds(g_bounds), eval_size(lfg.size()) {
        for (const auto& func : lfg) {
            jac_size += func.jac.nnz();
            hes_size += func.hes.nnz();
        }
        eval_buffer = FixedVector<double>(mesh->node_count * eval_size);
        jac_buffer = FixedVector<double>(mesh->node_count * jac_size);
        hes_buffer = FixedVector<double>(mesh->node_count * hes_size);
    };

    FixedVector<FunctionLFG> lfg;
    std::shared_ptr<Mesh> mesh;

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

    FixedVector<double> eval_buffer;
    FixedVector<double> jac_buffer;
    FixedVector<double> hes_buffer; // can lead to several memory consumption / order of a few GB, so maybe rethink this for large scale problems

    // fill eval_buffer, jac_buffer, hes_buffer
    virtual void callbackEval(const double* xu_nlp, const double* p) = 0;

    virtual void callbackJac(const double* xu_nlp, const double* p) = 0;

    virtual void callbackHes(const double* xu_nlp, const double* p) = 0;
    
    inline double getEvalL(const int offset) {
        return *(lfg[0].eval + offset * eval_size);
    }

    inline double getEvalF(const int f_index, const int offset) {
        return *(lfg[f_index_start + f_index].eval + offset * eval_size);
    }

    inline double getEvalG(const int g_index, const int offset) {
        return *(lfg[g_index_start + g_index].eval + offset * eval_size);
    }
};

class BoundarySweep {
public:

    BoundarySweep(FixedVector<FunctionMR>&& mr_in, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds, bool has_mayer,
                  int r_size, int x_size, int p_size)
    : mr(std::move(mr_in)), mesh(mesh), r_size(r_size), has_mayer(has_mayer), r_index_start((int) has_mayer),
      r_index_end(r_index_start + r_size), x_size(x_size), p_size(p_size), r_bounds(r_bounds) {
        int jac_buffer_size = 0;
        int hes_buffer_size = 0;
        for (const auto& func : mr) {
            jac_buffer_size += func.jac.nnz();
            hes_buffer_size += func.hes.nnz();
        }
        eval_buffer = FixedVector<double>(r_index_end);
        jac_buffer = FixedVector<double>(jac_buffer_size);
        hes_buffer = FixedVector<double>(hes_buffer_size);
    };

    // M, R :: assert x0 Size == xf Size
    FixedVector<FunctionMR> mr;
    std::shared_ptr<Mesh> mesh;

    FixedVector<Bounds> r_bounds;

    bool has_mayer = false;
    int r_index_start = 0;
    int r_index_end = 0;
    int r_size = 0; // assert; check with fixed initial states, these should not be contained here!

    int x_size = 0;
    int p_size = 0;

    FixedVector<double> eval_buffer;
    FixedVector<double> jac_buffer;
    FixedVector<double> hes_buffer;

    virtual void callbackEval(const double* x0_nlp, const double* xf_nlp, const double* p) = 0;

    virtual void callbackJac(const double* x0_nlp, const double* xf_nlp, const double* p) = 0;

    virtual void callbackHes(const double* x0_nlp, const double* xf_nlp, const double* p) = 0;

    inline double getEvalM() {
        return *(mr[0].eval);
    }

    inline double getEvalR(const int r_index) {
        return *(mr[r_index_start + r_index].eval);
    }
};

class Problem {
public:
    Problem(std::unique_ptr<FullSweep> full, std::unique_ptr<BoundarySweep> boundary, FixedVector<Bounds>& x_bounds,
               FixedVector<Bounds>& u_bounds, FixedVector<Bounds>& p_bounds, FixedVector<std::optional<double>>& x0_fixed,
               FixedVector<std::optional<double>>& xf_fixed)
    : full(std::move(full)), boundary(std::move(boundary)), x_size(this->full->x_size), u_size(this->full->u_size), p_size(this->full->p_size),
      x_bounds(std::move(x_bounds)), u_bounds(std::move(u_bounds)), p_bounds(std::move(p_bounds)),
      x0_fixed(std::move(x0_fixed)), xf_fixed(std::move(xf_fixed)) {
    };
    
    std::unique_ptr<FullSweep> full;
    std::unique_ptr<BoundarySweep> boundary;

    FixedVector<Bounds> x_bounds;
    FixedVector<Bounds> u_bounds;
    FixedVector<Bounds> p_bounds;

    FixedVector<std::optional<double>> x0_fixed; // set value if a state has a fixed initial value, remove the constraint from r()!
    FixedVector<std::optional<double>> xf_fixed; // set value if a state has a fixed final value, remove the constraint from r()!

    int x_size;
    int u_size;
    int p_size;
};

#endif // OPT_GDOP_PROBLEM_H
