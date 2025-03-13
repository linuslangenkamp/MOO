#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H

#include "util.h"

/* FullSweep: Evaluation of L(), f(), g() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T
 * + first and second derivatives
 */

// LFG - generic global function f(x, u, p, t)
// used for Lagrange term (L), dynamic (F), path (G)

struct JacobianLFG {
    // coordinate format jacobian for LFGH functions
    std::vector<std::tuple<uint32_t, double*>> dx;
    std::vector<std::tuple<uint32_t, double*>> du;
    std::vector<std::tuple<uint32_t, double*>> dp;
};

struct HessianLFG {
    // coordinate format hessian for LFG functions
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dx_dx;

    std::vector<std::tuple<uint32_t, uint32_t, double*>> du_dx;
    std::vector<std::tuple<uint32_t, uint32_t, double*>> du_du;

    std::vector<std::tuple<uint32_t, uint32_t, double*>> dp_dx;
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dp_du;
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dp_dp;
};

struct FunctionLFG {
    double* eval;
    JacobianLFG jacCOO;
    HessianLFG hessCOO;
};

// MR - generic boundary function r(x(t0), x(tf), p, t0, tf)
// used for Mayer term (M), boundary constraints (R)
 
struct JacobianMR {
    std::vector<std::tuple<uint32_t, double*>> dx0;
    std::vector<std::tuple<uint32_t, double*>> dxf;
    std::vector<std::tuple<uint32_t, double*>> dp;
};

struct HessianMR {
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dx0_dx0;
    
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dxf_dx0;
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dxf_dxf;

    std::vector<std::tuple<uint32_t, uint32_t, double*>> dp_dx0;
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dp_dxf;
    std::vector<std::tuple<uint32_t, uint32_t, double*>> dp_dp;
};


struct FunctionMR {
    double* eval;
    JacobianMR jacCOO;
    HessianMR hessCOO;
};

// TODO: where do we add nominals? or is the problem scaled automatically in OpenModelica?
// TODO: add a scaling / nominal wrapper, which updates the func-evals depending on nominal values
// TOOD: add scaling for variables
struct FullSweep {
    // Idea: Call Fullsweep.setValues(), Fullsweep.callEval(), callJac(), callHess() -> Just iterate over COO
    // function evals / diffs are on openmodelica side, 

    std::vector<FunctionLFG> functionsLFG;

    std::vector<std::tuple<double, double>> gBounds; // g^L <= g(x, u, p, t) <= g^U :: path constraints

    // L, F, G indices in functionLFGH vector
    uint32_t LIndexStart, LIndexEnd;
    uint32_t FIndexStart, FIndexEnd;
    uint32_t GIndexStart, GIndexEnd;

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
};

struct BoundarySweep {
    // M, R :: assert x0 Size == xf Size
    std::vector<FunctionMR> functionsMR;
    
    uint32_t MIndexStart, MIndexEnd;
    uint32_t RIndexStart, RIndexEnd;

    std::vector<std::tuple<double, double>> rBounds; // r^L <= r(x0, xf, p) dt <= r^U :: boundary constraints

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
};

struct Problem {
    FullSweep fullSweep;
    BoundarySweep boundarySweep;

    std::vector<std::tuple<double, double>> xBounds, uBounds;
    std::vector<std::tuple<double, double>> pBounds;

    uint32_t xSize;
    uint32_t uSize;
    uint32_t pSize;
};

#endif  // OPT_PROBLEM_H
