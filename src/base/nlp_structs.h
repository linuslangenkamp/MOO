#ifndef OPT_NLP_STRUCTS_H
#define OPT_NLP_STRUCTS_H

#include <vector>

struct Bounds {
    double lb;
    double ub;
};

struct JacobianSparsity {
    int col;
    double* value;
};

struct HessianSparsity {
    int row;
    int col;
    double* value;
};

// LFG - generic global function f(x, u, p, t)
// used for Lagrange term (L), dynamic (F), path (G) in GDOP

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

// MR - semi-generic boundary function r(x(t0), x(tf), p)
// used for Mayer term (M), boundary constraints (R) in GDOP

struct JacobianMR {
    // coordinate format jacobian for MR functions
    std::vector<JacobianSparsity> dx0;
    std::vector<JacobianSparsity> dxf;
    std::vector<JacobianSparsity> dp;
};

struct HessianMR {
    // coordinate format hessian for MR functions
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

#endif // OPT_NLP_STRUCTS_H
