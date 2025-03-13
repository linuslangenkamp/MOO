#ifndef OPT_PROBLEM_H
#define OPT_PROBLEM_H

#include "util.h"

/* FullSweep: Evaluation of L(), f(), g(), h() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T
 * + first and second derivatives
 */

// LFGH - generic global function f(x, u, p, t)
// used for Lagrange term (L), dynamic (F), path (G), integral constraints (H)

struct GradientLFGH {
    // coordinate format gradient for LFGH functions
    std::vector<std::tuple<int, gNumber*>> dx;
    std::vector<std::tuple<int, gNumber*>> du;
    std::vector<std::tuple<int, gNumber*>> dp;
    gNumber* dt;
};

struct HessianLFGH {
    // coordinate format hessian for LFGH functions
    std::vector<std::tuple<int, int, gNumber*>> dx_dx;

    std::vector<std::tuple<int, int, gNumber*>> du_dx;
    std::vector<std::tuple<int, int, gNumber*>> du_du;

    // ? bool: contains p?
    std::vector<std::tuple<int, int, gNumber*>> dp_dx;
    std::vector<std::tuple<int, int, gNumber*>> dp_du;
    std::vector<std::tuple<int, int, gNumber*>> dp_dp;

    // ? bool: contains t?
    std::vector<std::tuple<int, gNumber*>> dt_dx;
    std::vector<std::tuple<int, gNumber*>> dt_du;
    std::vector<std::tuple<int, gNumber*>> dt_dp;
    gNumber* dt_dt;
};

struct FunctionLFGH {
    gNumber* eval;
    GradientLFGH gradCOO;
    HessianLFGH hessCOO;
};
// MR - generic boundary function r(x(t0), x(tf), p, t0, tf)
// used for Mayer term (M), boundary constraints (R)
 
struct GradientMR {
    std::vector<std::tuple<int, gNumber*>> dx0;
    std::vector<std::tuple<int, gNumber*>> dxf;
    std::vector<std::tuple<int, gNumber*>> dp;
    gNumber* dt0;
    gNumber* dtf;
};

struct HessianMR {
    std::vector<std::tuple<int, int, gNumber*>> dx0_dx0;
    
    std::vector<std::tuple<int, int, gNumber*>> dxf_dx0;
    std::vector<std::tuple<int, int, gNumber*>> dxf_dxf;

    std::vector<std::tuple<int, int, gNumber*>> dp_dx0;
    std::vector<std::tuple<int, int, gNumber*>> dp_dxf;
    std::vector<std::tuple<int, int, gNumber*>> dp_dp;

    std::vector<std::tuple<int, gNumber*>> dt0_dx0;
    std::vector<std::tuple<int, gNumber*>> dt0_dxf;
    std::vector<std::tuple<int, gNumber*>> dt0_dp;
    gNumber* dt0_dt0;

    std::vector<std::tuple<int, gNumber*>> dtf_dx0;
    std::vector<std::tuple<int, gNumber*>> dtf_dxf;
    std::vector<std::tuple<int, gNumber*>> dtf_dp;
    gNumber* dtf_dt0;
    gNumber* dtf_dtf;
};


struct FunctionMR {
    gNumber* eval;
    GradientMR gradCOO;
    HessianMR hessCOO;
};

// E - generic even function e(x^{pre}(tf), x^{suc}(t0), p)
// used for event constraints / linkages
 
struct GradientE {
    // coordinate format gradient for event functions
    // TODO: maybe extend this to also incorporate x0pre, xfpre, x0suc, xfsuc? 
    // <=> r(x0pre, p) = x0pre - p == 0 => include p here, just overhead of 1 p and more simple structures
    std::vector<std::tuple<int, gNumber*>> dxfPre;
    std::vector<std::tuple<int, gNumber*>> dx0Suc;
    std::vector<std::tuple<int, gNumber*>> dp;
};

struct HessianE {
    // coordinate format hessian for event functions
    std::vector<std::tuple<int, int, gNumber*>> dxfPre_dxfPre;

    std::vector<std::tuple<int, int, gNumber*>> dx0Suc_dxfPre;
    std::vector<std::tuple<int, int, gNumber*>> dx0Suc_dx0Suc;

    std::vector<std::tuple<int, int, gNumber*>> dp_dxfPre;
    std::vector<std::tuple<int, int, gNumber*>> dp_dx0Suc;
    std::vector<std::tuple<int, int, gNumber*>> dp_dp;
};

struct FunctionE {
    gNumber* eval;
    GradientE gradCOO;
    HessianE hessCOO;
};

// TODO: where do we add nominals? or is the problem scaled automatically in OpenModelica?
// TODO: add a scaling / nominal wrapper, which updates the func-evals depending on nominal values
// TOOD: add scaling for variables
struct FullSweep {
    // Idea: Call Fullsweep.setValues(), Fullsweep.callEval(), callGrad(), callHess() -> Just iterate over COO
    // function evals / diffs are on openmodelica side, 

    // definitely make this *FunctionLFGH, so they can be reused
    std::vector<std::shared_ptr<FunctionLFGH>> functionsLFGH;

    std::vector<std::tuple<gNumber, gNumber>> gBounds; // g^L <= g(x, u, p, t) <= g^U :: path constraints
    std::vector<std::tuple<gNumber, gNumber>> hBounds; // h^L <= \int_{t_0}^{t_f} h(x, u, p, t) dt <= h^U :: integral constraints

    int xSize;
    int uSize;
    int pSize;
    int tSize;

    // L, F, G, H indices in functionLFGH vector
    int LIndexStart, LIndexEnd;
    int FIndexStart, FIndexEnd;
    int GIndexStart, GIndexEnd;
    int HIndexStart, HIndexEnd;

    // create this based on some input data for the Fullsweep. Create the Grad
    FullSweep() {
        // now just use as: fillInputData() -> iterate over function COO's
    };

    inline void fillInputData(gNumber* x, gNumber* u, gNumber* p, gNumber* t) {
        // TODO: think about this part
        // unscale
        // fill the OM buffer
    };

    inline void callEval() {
        // TODO: think about this part
        // scale
   }

    inline void callGrad() {
        // TODO: think about this part
        // scale
    }

    inline void callHess() {
        // TODO: think about this part
        // scale
    }
};

struct BoundarySweep {
    // M, R :: guard x0 Size == xf Size
    std::vector<std::shared_ptr<FunctionMR>> functionsMR;
    
    int xSize;
    int pSize;
    int t0Size;
    int tfSize;

    int MIndexStart, MIndexEnd;
    int RIndexStart, RIndexEnd;

    std::vector<std::tuple<gNumber, gNumber>> rBounds; // r^L <= r(x0, xf, p, t0, tf) dt <= r^U :: boundary constraints

    inline void fillInputData(gNumber* x0, gNumber* xf, gNumber* p, gNumber* t0, gNumber* tf) {
        // TODO: think about this part
        // unscale
        // fill the OM buffer
    };
    
    inline void callEval() {
        // TODO: think about this part
        // scale
   }

    inline void callGrad() {
        // TODO: think about this part
        // scale
    }

    inline void callHess() {
        // TODO: think about this part
        // scale
    }
};

// Linkage stuff
struct PhaseConnection {
    // t_{0, suc} - t_{f, pre}

    // no need for fancy COO
    const gNumber dt0Suc = 1;
    const gNumber dtfPre = -1;

    inline gNumber eval(const gNumber* tfPre, const gNumber* t0Suc) {
        return *t0Suc - *tfPre;
    }
};

struct Linkage {
    std::vector<std::shared_ptr<FunctionE>> functionsE;       // event functions
    std::vector<std::tuple<gNumber, gNumber>> eBounds;        // e^L <= e(x_{f, pre}, x_{0, suc}, p) <= e^U :: event constraints

    // optional, if the time horion is fixed and t0Suc = tfPre => just set hard in box constraints
    std::shared_ptr<PhaseConnection> connection;              // phase connection
    std::shared_ptr<std::tuple<gNumber, gNumber>> connBounds; // t^L <= t_{0, suc} - t_{f, pre}      <= t^U :: phase connections

    // TODO: do i really need these sizes everywhere? or just set to phase
    int xfPreSize, x0SucSize, pSize;
    int phasePre, phaseSuc;
};

struct PhaseVariables {
    std::vector<std::tuple<gNumber, gNumber>> xBounds, uBounds;
    std::tuple<gNumber, gNumber> t0Bounds, tfBounds;
};

struct GlobalParameters {
    std::vector<std::tuple<gNumber, gNumber>> pBounds;
};

struct Phase {
    std::shared_ptr<FullSweep> fullSweep;
    std::shared_ptr<BoundarySweep> boundarySweep;
    std::shared_ptr<PhaseVariables> variables;
};

struct Model {
    std::vector<std::shared_ptr<Phase>> phases;
    std::vector<std::shared_ptr<Linkage>> linkages;
    std::shared_ptr<GlobalParameters> parameters;
};


// TODO: think about the ptrs and where they are advantagous in the structures?
// Where do we need to duplicate Phases / Linkages?
// also duplicate functions, how to remember or from OM side
// dont create new COO arrays, if we know that these are duplicates

#endif  // OPT_PROBLEM_H
