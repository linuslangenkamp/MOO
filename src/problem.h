#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H

#include <optional>
#include <memory>
#include <vector>

/* FullSweep: Evaluation of L(), f(), g() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T
 * + first and second derivatives
 */

// LFG - generic global function f(x, u, p, t)
// used for Lagrange term (L), dynamic (F), path (G)

struct JacobianSparsity {
    int index;
    double* value;
};

struct HessianSparsity {
    int index1;
    int index2;
    double* value;
};

struct Bounds {
    double lb;
    double ub;
};

struct JacobianLFG {
    // coordinate format jacobian for LFGH functions
    std::vector<JacobianSparsity> dx;
    std::vector<JacobianSparsity> du;
    std::vector<JacobianSparsity> dp;
};

struct HessianLFG {
    // coordinate format hessian for LFG functions
    std::vector<HessianSparsity> dx_dx;
    std::vector<HessianSparsity> du_dx;
    std::vector<HessianSparsity> du_du;
    std::vector<HessianSparsity> dp_dx;
    std::vector<HessianSparsity> dp_du;
    std::vector<HessianSparsity> dp_dp;
};

struct FunctionLFG {
    double* eval;
    JacobianLFG jac;
    HessianLFG hes;
};

// MR - generic boundary function r(x(t0), x(tf), p, t0, tf)
// used for Mayer term (M), boundary constraints (R)
 
struct JacobianMR {
    std::vector<JacobianSparsity> dx0;
    std::vector<JacobianSparsity> dxf;
    std::vector<JacobianSparsity> dp;
};

struct HessianMR {
    std::vector<HessianSparsity> dx0_dx0;
    
    std::vector<HessianSparsity> dxf_dx0;
    std::vector<HessianSparsity> dxf_dxf;

    std::vector<HessianSparsity> dp_dx0;
    std::vector<HessianSparsity> dp_dxf;
    std::vector<HessianSparsity> dp_dp;
};

struct FunctionMR {
    double* eval;
    JacobianMR jac;
    HessianMR hes;
};

struct FullSweep {
    // Idea: Call Fullsweep.setValues(), Fullsweep.callEval(), callJac(), callHess() -> Just iterate over COO
    // function evals / diffs are on openmodelica side, 

    std::vector<FunctionLFG> lfg;

    std::vector<Bounds> g_bounds; // g^L <= g(x, u, p, t) <= g^U :: path constraints

    //  F, G indices in lfg vector
    int f_index_start, f_index_end;
    int g_index_start, g_index_end;
    int f_size;
    int g_size;
    int fg_size;
    bool has_lagrange;

    // sizes of 1 chunk in the data array - which has to be defined somewhere, maybe in NLP
    int eval_size; // eval_size == g_index_end == number of functions in lfg
    int jac_size;
    int hes_size;

    // create this based on some input data for the Fullsweep. Create the Grad
    FullSweep() {
        // now just use as: fillInputData() -> iterate over function COO's
    };

    inline void fillInputData(double* x, double* u, double* p, double* t) {
        // TODO: think about this part
        // input unscaled values for x, u, p, t at some (i, j)
        // fill the OM buffer (should be scoped / parallelizable / localized)
    };

    inline void callEval(double* evalData) {
        // TODO: think about this part
        // get data array from solver
        // calculate eval 
   }

    inline void callJac(double* jacData) {
        // TODO: think about this part
    }

    inline void callHess(double* hesData) {
        // TODO: think about this part
    }

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

struct BoundarySweep {
    // M, R :: assert x0 Size == xf Size
    std::vector<FunctionMR> mr;

    bool has_mayer;
    int r_index_start, r_index_end;
    int r_size; // assert; check with fixed initial states, these should not be contained here!
    
    std::vector<Bounds> r_bounds; // r^L <= r(x0, xf, p) <= r^U :: boundary constraints

    inline void fillInputData(double* x0, double* xf, double* p) {
        // TODO: think about this part
        // unscale
        // fill the OM buffer
    };
    
    inline void callEval(double* evalData) {
        // TODO: think about this part
        // get data array from solver
        // calculate eval 
   }

    inline void callJac(double* jacData) {
        // TODO: think about this part
    }

    inline void callHess(double* hesData) {
        // TODO: think about this part
    }

    inline double getEvalM() {
        return *(mr[0].eval);
    }

    inline double getEvalR(const int r_index) {
        return *(mr[r_index_start + r_index].eval);
    }
};

struct Problem {
    FullSweep full;
    BoundarySweep boundary;

    std::vector<Bounds> x_bounds;
    std::vector<Bounds> u_bounds;
    std::vector<Bounds> p_bounds;

    std::vector<std::optional<double>> x0_fixed; // set value if a state has a fixed initial value, remove the constraint from r()!
    std::vector<std::optional<double>> xf_fixed; // set value if a state has a fixed final value, remove the constraint from r()!

    int x_size;
    int u_size;
    int p_size;
};

#endif  // OPT_PROBLEM_H
