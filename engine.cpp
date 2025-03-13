#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <set>
#include <chrono>
#include <symengine/llvm_double.h>
#include <symengine/parser.h>
#include "engine.h"

// g++ -o engine engine.cpp engine.h -O2 -lsymengine -lgmp -lLLVM-18 -g  && ./engine 

using namespace SymEngine;

// indices for LFGH expressions
constexpr int L_INDEX_START = 0;
constexpr int L_INDEX_END = 1;
constexpr int F_INDEX_START = 1;
constexpr int F_INDEX_END = 2;
constexpr int G_INDEX_START = 2;
constexpr int G_INDEX_END = 3;
constexpr int H_INDEX_START = 3;
constexpr int H_INDEX_END = 4;

// indices for MR expressions
constexpr int M_INDEX_START = 0;
constexpr int M_INDEX_END = 1;
constexpr int R_INDEX_START = 1;
constexpr int R_INDEX_END = 2;

// SymEngine -> LLVM -> binary options
constexpr bool SYM_LLVM_USE_CSE = true;
constexpr int SYM_LLVM_OPT_LEVEL = 3;

inline std::unique_ptr<SymEngine::LLVMDoubleVisitor> compileVectorExprLLVM(
    const std::vector<RCP<const Basic>>& vecexpr, const std::vector<RCP<const Basic>>& vars, bool use_cse, int opt_level) {
    auto llvm_func = std::make_unique<SymEngine::LLVMDoubleVisitor>();
    llvm_func->init(vars, vecexpr, use_cse, opt_level); // TODO: investigate use_cse, opt_level?
    return llvm_func;
}

// TODO: JUST NAME IT AS HASH OF THE EXPRESSION VECTOR + INPUT VECTOR --> HASH AGAIN AND IF == -> READ FROM FILE
void saveBinaryLLVM(std::unique_ptr<SymEngine::LLVMDoubleVisitor>& llvmVisitor, std::string& filepath) {
    auto binary = llvmVisitor->dumps();
    std::ofstream file(filepath, std::ios::binary);
    if (file.is_open()) {
        file.write(binary.data(), binary.size());
        file.close();
    }
    else {
        std::cerr << "Error: Failed to write binary to " << filepath << "." << std::endl;
    }
}

std::unique_ptr<SymEngine::LLVMDoubleVisitor> readBinaryLLVM(std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary | std::ios::ate);
    if (!file.is_open()) {
        std::cerr << "Error: Failed to open file " << filepath << " to read binary data." << std::endl;
    }

    size_t fileSize = file.tellg();
    file.seekg(0, std::ios::beg);
    std::string binaryData(fileSize, '\0');
    file.read(&binaryData[0], fileSize);
    file.close();

    if (file.fail()) {
        std::cerr << "Error: Failed to read binary data from " << filepath << "." << std::endl;
    }

    std::unique_ptr<SymEngine::LLVMDoubleVisitor> llvm_func = std::make_unique<LLVMDoubleVisitor>();
    llvm_func->loads(binaryData);

    return llvm_func;
}

inline bool isZero(RCP<const Basic> expr) { return (!(is_zero(*expr) == tribool::tritrue)); }
    
struct PhaseVariables {
    std::vector<RCP<const Basic>> x, u, x0, xf;
    RCP<const Basic> t, t0, tf;
    std::vector<std::tuple<double, double>> xBounds, uBounds;
    std::tuple<double, double> t0Bounds, tfBounds;
    PhaseVariables(std::vector<RCP<const Basic>> x = {}, std::vector<RCP<const Basic>> u = {}, std::vector<RCP<const Basic>> x0 = {},
              std::vector<RCP<const Basic>> xf = {}, RCP<const Basic> t = {},  RCP<const Basic> t0 = {}, RCP<const Basic> tf = {},
              std::vector<std::tuple<double, double>> xBounds = {}, std::vector<std::tuple<double, double>> uBounds = {}, std::vector<std::tuple<double, double>> pBounds = {},
              std::tuple<double, double> t0Bounds = {}, std::tuple<double, double> tfBounds = {})
            : x(x), u(u), x0(x0), xf(xf), t(t), t0(t0), tf(tf), xBounds(xBounds), uBounds(uBounds), t0Bounds(t0Bounds), tfBounds(tfBounds) {
            }
};

struct GlobalParameters {
    std::vector<RCP<const Basic>> p;
    std::vector<std::tuple<double, double>> pBounds;

    GlobalParameters(std::vector<RCP<const Basic>> p = {}, std::vector<std::tuple<double, double>> pBounds = {})
    : p(p), pBounds(pBounds) {}
};

// symbol in symbols vector -> f(x, u, p, t, mp), but mp is fixed for all new inputs
// altough being compiled as variable and thus input -> reading from file with new parameters is possible!
struct ModelParameters {
    std::vector<RCP<const Basic>> mp;
    std::vector<double> mpValues;

    ModelParameters(std::vector<RCP<const Basic>> mp = {}, std::vector<double> mpValues = {})
    : mp(mp), mpValues(mpValues) {
    }
};

struct PhaseFunctions {
    RCP<const Basic> M, L;
    std::vector<RCP<const Basic>> f, g, h, r;
    std::vector<std::tuple<double, double>> gBounds, hBounds, rBounds;
    PhaseFunctions(RCP<const Basic> M = {}, RCP<const Basic> L = {}, std::vector<RCP<const Basic>> f = {}, std::vector<RCP<const Basic>> g = {},
              std::vector<RCP<const Basic>> h = {}, std::vector<RCP<const Basic>> r = {}, std::vector<std::tuple<double, double>> gBounds = {},
              std::vector<std::tuple<double, double>> hBounds = {}, std::vector<std::tuple<double, double>> rBounds = {})
            : M(M), L(L), f(f), g(g), h(h), r(r), gBounds(gBounds), hBounds(hBounds), rBounds(rBounds) {
            }
};

/* FullSweep: Evaluation of L(), f(), g(), h() for a given z_{i,j} = (x_{i,j}, u_{i,j}, p, t_{i,j})^T
 * + first and second derivatives
 */

// LFGH - generic global function f(x, u, p, t)
// used for Lagrange term (L), dynamic (F), path (G), integral constraints (H)

struct GradientLFGH {
    // coordinate format gradient for LFGH functions
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dx;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> du;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dp;
    std::tuple<double*, RCP<const Basic>> dt;

    std::vector<RCP<const Basic>> toVector() {
        std::vector<RCP<const Basic>> vectorGrad;
        for (auto& coo : dx) { vectorGrad.push_back(std::get<2>(coo)); } 
        for (auto& coo : du) { vectorGrad.push_back(std::get<2>(coo)); } 
        for (auto& coo : dp) { vectorGrad.push_back(std::get<2>(coo)); } 
        if (!(std::get<1>(dt)).is_null()) {
            vectorGrad.push_back(std::get<1>(dt));
        }
        return vectorGrad;
    }
};

struct HessianLFGH {
    // coordinate format hessian for LFGH functions
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dx_dx;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> du_dx;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> du_du;

    // ? bool: contains p?
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dx;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_du;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dp;

    // ? bool: contains t?
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dt_dx;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dt_du;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dt_dp;

    std::tuple<double*, RCP<const Basic>> dt_dt;

    std::vector<RCP<const Basic>> toVector() {
        std::vector<RCP<const Basic>> vectorHess;
        for (auto& coo : dx_dx) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : du_dx) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : du_du) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dx) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_du) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dp) { vectorHess.push_back(std::get<3>(coo)); } 

        for (auto& coo : dt_dx) { vectorHess.push_back(std::get<2>(coo)); } 
        for (auto& coo : dt_du) { vectorHess.push_back(std::get<2>(coo)); } 
        for (auto& coo : dt_dp) { vectorHess.push_back(std::get<2>(coo)); } 
        if (!std::get<1>(dt_dt).is_null()) {
            vectorHess.push_back(std::get<1>(dt_dt));
        }
        return vectorHess;
    }
};

// MR - generic boundary function r(x(t0), x(tf), p, t0, tf)
// used for Mayer term (M), boundary constraints (R)
 
struct GradientMR {
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dx0;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dxf;

    std::vector<std::tuple<int, double*, RCP<const Basic>>> dp;

    std::tuple<double*, RCP<const Basic>> dt0;
    std::tuple<double*, RCP<const Basic>> dtf;

    std::vector<RCP<const Basic>> toVector() {
        std::vector<RCP<const Basic>> vectorGrad;
        for (auto& coo : dx0) { vectorGrad.push_back(std::get<2>(coo)); } 
        for (auto& coo : dxf) { vectorGrad.push_back(std::get<2>(coo)); } 
        for (auto& coo : dp) { vectorGrad.push_back(std::get<2>(coo)); } 
        if (!(std::get<1>(dt0)).is_null()) {
            vectorGrad.push_back(std::get<1>(dt0));
        }
        if (!(std::get<1>(dtf)).is_null()) {
            vectorGrad.push_back(std::get<1>(dtf));
        }
        return vectorGrad;
    }
};

struct HessianMR {
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dx0_dx0;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dxf_dx0;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dxf_dxf;

    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dx0;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dxf;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dp;

    std::vector<std::tuple<int, double*, RCP<const Basic>>> dt0_dx0;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dt0_dxf;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dt0_dp;
    std::tuple<double*, RCP<const Basic>> dt0_dt0;

    std::vector<std::tuple<int, double*, RCP<const Basic>>> dtf_dx0;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dtf_dxf;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dtf_dp;
    std::tuple<double*, RCP<const Basic>> dtf_dt0;
    std::tuple<double*, RCP<const Basic>> dtf_dtf;

    std::vector<RCP<const Basic>> toVector() {
        std::vector<RCP<const Basic>> vectorHess;
        for (auto& coo : dx0_dx0) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dxf_dx0) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dxf_dxf) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dx0) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dxf) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dp) { vectorHess.push_back(std::get<3>(coo)); } 

        for (auto& coo : dt0_dx0) { vectorHess.push_back(std::get<2>(coo)); } 
        for (auto& coo : dt0_dxf) { vectorHess.push_back(std::get<2>(coo)); } 
        for (auto& coo : dt0_dp) { vectorHess.push_back(std::get<2>(coo)); } 
        if (!std::get<1>(dt0_dt0).is_null()) {
            vectorHess.push_back(std::get<1>(dt0_dt0));
        }

        for (auto& coo : dtf_dx0) { vectorHess.push_back(std::get<2>(coo)); } 
        for (auto& coo : dtf_dxf) { vectorHess.push_back(std::get<2>(coo)); } 
        for (auto& coo : dtf_dp) { vectorHess.push_back(std::get<2>(coo)); } 
        if (!std::get<1>(dtf_dt0).is_null()) {
            vectorHess.push_back(std::get<1>(dtf_dt0));
        }
        if (!std::get<1>(dtf_dtf).is_null()) {
            vectorHess.push_back(std::get<1>(dtf_dtf));
        }
        return vectorHess;
    }
};

// E - generic even function e(x^{pre}(tf), x^{suc}(t0), p)
// used for event constraints / linkages
 
struct GradientE {
    // coordinate format gradient for event / linkage functions
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dxfPre;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dx0Suc;
    std::vector<std::tuple<int, double*, RCP<const Basic>>> dp;

    std::vector<RCP<const Basic>> toVector() {
        std::vector<RCP<const Basic>> vectorGrad;
        for (auto& coo : dxfPre) { vectorGrad.push_back(std::get<2>(coo)); } 
        for (auto& coo : dx0Suc) { vectorGrad.push_back(std::get<2>(coo)); } 
        for (auto& coo : dp) { vectorGrad.push_back(std::get<2>(coo)); } 
        return vectorGrad;
    }
};

struct HessianE {
    // coordinate format hessian for event / linkage functions
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dxfPre_dxfPre;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dx0Suc_dxfPre;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dx0Suc_dx0Suc;

    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dxfPre;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dx0Suc;
    std::vector<std::tuple<int, int, double*, RCP<const Basic>>> dp_dp;

    std::vector<RCP<const Basic>> toVector() {
        std::vector<RCP<const Basic>> vectorHess;
        for (auto& coo : dxfPre_dxfPre) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dx0Suc_dxfPre) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dx0Suc_dx0Suc) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dxfPre) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dx0Suc) { vectorHess.push_back(std::get<3>(coo)); } 
        for (auto& coo : dp_dp) { vectorHess.push_back(std::get<3>(coo)); } 
        return vectorHess;
    }
};

// gradient for standard non-time variables, i.e. x, u, p
void calcGradVar(const RCP<const Basic>& symExpr, const std::vector<RCP<const Basic>>& symVars, std::vector<std::tuple<int, double*, RCP<const Basic>>>& gradientCOO, int& nnzGrad) {
    for (int i = 0; i < symVars.size(); i++) {         
        if (free_symbols(*symExpr).count(symVars[i])) {
            auto symDer = symExpr->diff(rcp_dynamic_cast<const Symbol>(symVars[i]));
            if (isZero(symDer)) {
                gradientCOO.push_back(std::tuple<int, double*, RCP<const Basic>>{i, nullptr, symDer});
                nnzGrad++; 
            }
        }
    }
}

// hessian for pairs of standard non-time variables, i.e. x, u, p
void calcHessVarVar(const std::vector<std::tuple<int, double*, RCP<const Basic>>>& gradientCOO, const std::vector<RCP<const Basic>>& symVars, 
                    std::vector<std::tuple<int, int, double*, RCP<const Basic>>>& hessianCOO, int& nnzHess, const bool sameVariable) {
    int gradIndex = 0;
    for (const auto& coo : gradientCOO) {
        auto row = std::get<0>(coo);
        auto symDer1 = std::get<2>(coo);
        for (int col = 0; col < (sameVariable ? row + 1 : symVars.size()); col++) {
            if (free_symbols(*symDer1).count(symVars[col])) {
                auto symDer2 = symDer1->diff(rcp_dynamic_cast<const Symbol>(symVars[col]));
                if (isZero(symDer2)) {
                    hessianCOO.push_back(std::tuple<int, int, double*, RCP<const Basic>>{row, col, nullptr, symDer2});
                    nnzHess++;
                }
            }
        }
    }
}

// hessian for standard non-time variables, i.e. x, u, p and w.r.t. time
void calcHessTimeVar(const std::tuple<double*, RCP<const Basic>>& gradientCOO, const std::vector<RCP<const Basic>>& symVars,
                     std::vector<std::tuple<int, double*, RCP<const Basic>>>& hessianCOO, int& nnzHess) {
    for (int i = 0; i < symVars.size(); i++) {         
        if (free_symbols(*std::get<1>(gradientCOO)).count(symVars[i])) {
            auto symDer2 = std::get<1>(gradientCOO)->diff(rcp_dynamic_cast<const Symbol>(symVars[i]));
            if (isZero(symDer2)) {
                hessianCOO.push_back(std::tuple<int, double*, RCP<const Basic>>{i, nullptr, symDer2});
                nnzHess++; 
            }
        }
    }
}

struct FunctionLFGH {
    RCP<const Basic> symExpr;
    GradientLFGH gradCOO;
    HessianLFGH hessCOO;

    FunctionLFGH(const RCP<const Basic> expr) : symExpr(expr) {
    }

    void calcGrad(const PhaseVariables& variables, const GlobalParameters& parameters, int& nnzGrad) {
        calcGradVar(symExpr, variables.x, gradCOO.dx, nnzGrad);
        calcGradVar(symExpr, variables.u, gradCOO.du, nnzGrad);
        calcGradVar(symExpr, parameters.p, gradCOO.dp, nnzGrad);

        if (!variables.t.is_null() && free_symbols(*symExpr).count(variables.t)) {
            auto symDer = symExpr->diff(rcp_dynamic_cast<const Symbol>(variables.t));
            if (isZero(symDer)) {
                gradCOO.dt = std::tuple<double*, RCP<const Basic>>{nullptr, symDer};
                nnzGrad++;
            }    
        }
    }

    void setGradPtr(int& offset, double* start) {
        for (auto& coo : gradCOO.dx) { std::get<1>(coo) = &start[offset++]; } 
        for (auto& coo : gradCOO.du) { std::get<1>(coo) = &start[offset++]; } 
        for (auto& coo : gradCOO.dp) { std::get<1>(coo) = &start[offset++]; } 

        if (!std::get<1>(gradCOO.dt).is_null()) {
            std::get<0>(gradCOO.dt) = &start[offset++]; 
        }
    }

    void calcHess(const PhaseVariables& variables, const GlobalParameters& parameters, int& nnzHess) {
        calcHessVarVar(gradCOO.dx, variables.x, hessCOO.dx_dx, nnzHess, true);
        calcHessVarVar(gradCOO.du, variables.x, hessCOO.du_dx, nnzHess, false);
        calcHessVarVar(gradCOO.du, variables.u, hessCOO.du_du, nnzHess, true);
        calcHessVarVar(gradCOO.dp, variables.x, hessCOO.dp_dx, nnzHess, false);
        calcHessVarVar(gradCOO.dp, variables.u, hessCOO.dp_du, nnzHess, false);
        calcHessVarVar(gradCOO.dp, parameters.p, hessCOO.dp_dp, nnzHess, true);
        
        if (std::get<0>(gradCOO.dt) != NULL) {
            calcHessTimeVar(gradCOO.dt, variables.x, hessCOO.dt_dx, nnzHess);
            calcHessTimeVar(gradCOO.dt, variables.u, hessCOO.dt_du, nnzHess);
            calcHessTimeVar(gradCOO.dt, parameters.p, hessCOO.dt_dp, nnzHess);

            if (!variables.t.is_null() && free_symbols(*std::get<1>(gradCOO.dt)).count(variables.t)) {
                auto symDer2 = std::get<1>(gradCOO.dt)->diff(rcp_dynamic_cast<const Symbol>(variables.t));
                if (isZero(symDer2)) {
                    hessCOO.dt_dt = std::tuple<double*, RCP<const Basic>>{nullptr, symDer2};
                    nnzHess++;
                }    
            }
        }
    }

    void setHessPtr(int& offset, double* start) {
        for (auto& coo : hessCOO.dx_dx) { std::get<2>(coo) = &start[offset++]; } 
        for (auto& coo : hessCOO.du_dx) { std::get<2>(coo) = &start[offset++]; }
        for (auto& coo : hessCOO.du_du) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dx) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_du) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dp) { std::get<2>(coo) = &start[offset++]; } 

        // index 1 for (index, double*), since t is implicit
        for (auto& coo : hessCOO.dt_dx) { std::get<1>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dt_du) { std::get<1>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dt_dp) { std::get<1>(coo) = &start[offset++]; }  

        if (!std::get<1>(hessCOO.dt_dt).is_null()) {
            std::get<0>(hessCOO.dt_dt) = &start[offset++];
        }
    }
};

struct FullSweep {
    // L, F, G, H
    std::vector<std::unique_ptr<FunctionLFGH>> functionsLFGH; 
    int ranges[5]; // {0, dimL, dimL + dimF, dimL + dimF + dimG, dimL + dimF + dimG + dimH}
    int xuSize;
    int pSize;
    int tSize;

    // evaluation with compiled code
    std::unique_ptr<double[]> inputData; // data input to llvm callbacks, must be filled before every call
    std::unique_ptr<double[]> evalData, gradData, hessData;  // data filled from llvm callbacks
    std::unique_ptr<SymEngine::LLVMDoubleVisitor> llvmEval, llvmGrad, llvmHess; // llvm callbacks

    FullSweep(const PhaseVariables& variables, const GlobalParameters& parameters, const ModelParameters& mparameters, const PhaseFunctions& functions) {
        initVariableSizes(variables, parameters);
        initFuncIndices(functions);
        initFunctions(functions);
        initDerivatives(variables, parameters);
        compileLLVM(variables, parameters, mparameters);
        // now just use as: fillInputData() -> iterate over function COO's
    };

    void initVariableSizes(const PhaseVariables& variables, const GlobalParameters& parameters) {
        xuSize = variables.x.size() + variables.u.size();
        pSize = parameters.p.size();
        tSize = variables.t.is_null() ? 0 : 1;
    }

    void initFuncIndices(const PhaseFunctions& functions) {
        ranges[0] = 0;
        ranges[1] = ranges[0] + !functions.L.is_null();
        ranges[2] = ranges[1] + functions.f.size();
        ranges[3] = ranges[2] + functions.g.size();
        ranges[4] = ranges[3] + functions.h.size();

        // memory for output of function evaluation
        evalData = std::make_unique<double[]>(ranges[4]);
    }
    
    void initFunctions(const PhaseFunctions& functions) {
        if (!functions.L.is_null()) {
            functionsLFGH.push_back(std::make_unique<FunctionLFGH>(functions.L));
        }
        for (auto& dynamic : functions.f) {
            functionsLFGH.push_back(std::make_unique<FunctionLFGH>(dynamic));
        }
        for (auto& path : functions.g) {
            functionsLFGH.push_back(std::make_unique<FunctionLFGH>(path));
        }
        for (auto& integral : functions.h) {
            functionsLFGH.push_back(std::make_unique<FunctionLFGH>(integral));
        }
    }

    void initDerivatives(const PhaseVariables& variables, const GlobalParameters& parameters) {  
        int nnzGrad = 0;
        for (auto& function : functionsLFGH) {
            function->calcGrad(variables, parameters, nnzGrad);
        }
        gradData = std::make_unique<double[]>(nnzGrad);
        int indexGrad = 0;
        for (auto& function : functionsLFGH) {
            function->setGradPtr(indexGrad, gradData.get());
        }
                
        int nnzHess = 0;
        for (auto& function : functionsLFGH) {
            function->calcHess(variables, parameters, nnzHess);
        }
        hessData = std::make_unique<double[]>(nnzHess);
        int indexHess = 0;
        for (auto& function : functionsLFGH) {
            function->setHessPtr(indexHess, hessData.get());
        }
    }

    void compileLLVM(const PhaseVariables& variables, const GlobalParameters& parameters, const ModelParameters& mparameters) {
        // init variable vector, memory for inputs
        std::vector<RCP<const Basic>> symVariables;
        std::vector<RCP<const Basic>> symEvalVector;
        std::vector<RCP<const Basic>> symGradVector;
        std::vector<RCP<const Basic>> symHessVector;

        for (auto& xvar: variables.x) {
            symVariables.push_back(xvar);
        }
        for (auto& uvar: variables.u) {
            symVariables.push_back(uvar);
        }
        for (auto& pvar: parameters.p) {
            symVariables.push_back(pvar);
        }

        if (!variables.t.is_null()) {
            symVariables.push_back(variables.t);
        }
        size_t variablesSize = symVariables.size();

        for (auto& mpvar: mparameters.mp) {
            symVariables.push_back(mpvar);
        }
        
        inputData = std::make_unique<double[]>(symVariables.size());
        for (size_t i = 0; i < mparameters.mp.size(); i++) {
            inputData[variablesSize + i] = mparameters.mpValues[i];
        }

        for (auto& function : functionsLFGH) {
            symEvalVector.push_back(function->symExpr);
            auto vectorGrad = function->gradCOO.toVector();
            auto vectorHess = function->hessCOO.toVector();
            symGradVector.insert(symGradVector.end(), vectorGrad.begin(), vectorGrad.end());
            symHessVector.insert(symHessVector.end(), vectorHess.begin(), vectorHess.end());
        }
        llvmEval = compileVectorExprLLVM(symEvalVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
        llvmGrad = compileVectorExprLLVM(symGradVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
        llvmHess = compileVectorExprLLVM(symHessVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
    }

    inline void fillInputData(double* xu, double* p, double* t) {
        // fills the input data array with new values (x, u, p, t, (mp)) - mp values are given as constructed
        for (int i = 0; i < xuSize; i++) {
            inputData[i] = xu[i];
        }
        for (int i = 0; i < pSize; i++) {
            inputData[xuSize + i] = p[i];
        }
        if (tSize > 0) {
            inputData[xuSize + pSize] = *t;
        }
    };

    inline void callEval() {
        llvmEval->call(evalData.get(), inputData.get());
    }

    inline void callGrad() {
        llvmGrad->call(gradData.get(), inputData.get());
    }

    inline void callHess() {
        llvmHess->call(hessData.get(), inputData.get());
    }
};

struct FunctionMR {
    RCP<const Basic> symExpr;
    GradientMR gradCOO;
    HessianMR hessCOO;

    FunctionMR(const RCP<const Basic> expr) : symExpr(expr) {
    }

    void calcGrad(const PhaseVariables& variables, const GlobalParameters& parameters, int& nnzGrad) {
        calcGradVar(symExpr, variables.x0, gradCOO.dx0, nnzGrad);
        calcGradVar(symExpr, variables.xf, gradCOO.dxf, nnzGrad);
        calcGradVar(symExpr, parameters.p, gradCOO.dp, nnzGrad);

        if (!variables.t0.is_null() && free_symbols(*symExpr).count(variables.t0)) {
            auto symDer = symExpr->diff(rcp_dynamic_cast<const Symbol>(variables.t0));
            if (isZero(symDer)) {
                gradCOO.dt0 = std::tuple<double*, RCP<const Basic>>{nullptr, symDer};
                nnzGrad++;
            }    
        }
        if (!variables.tf.is_null() && free_symbols(*symExpr).count(variables.tf)) {
            auto symDer = symExpr->diff(rcp_dynamic_cast<const Symbol>(variables.tf));
            if (isZero(symDer)) {
                gradCOO.dtf = std::tuple<double*, RCP<const Basic>>{nullptr, symDer};
                nnzGrad++;
            }    
        }
    }

    void setGradPtr(int& offset, double* start) {
        for (auto& coo : gradCOO.dx0) { std::get<1>(coo) = &start[offset++]; } 
        for (auto& coo : gradCOO.dxf) { std::get<1>(coo) = &start[offset++]; } 
        for (auto& coo : gradCOO.dp) { std::get<1>(coo) = &start[offset++]; } 

        if (!std::get<1>(gradCOO.dt0).is_null()) {
            std::get<0>(gradCOO.dt0) = &start[offset++];
        }
        if (!std::get<1>(gradCOO.dtf).is_null()) {
            std::get<0>(gradCOO.dtf) = &start[offset++];
        }
    }

    void calcHess(const PhaseVariables& variables, const GlobalParameters& parameters, int& nnzHess) {
        calcHessVarVar(gradCOO.dx0, variables.x0, hessCOO.dx0_dx0, nnzHess, true);
        calcHessVarVar(gradCOO.dxf, variables.x0, hessCOO.dxf_dx0, nnzHess, false);
        calcHessVarVar(gradCOO.dxf, variables.xf, hessCOO.dxf_dxf, nnzHess, true);
        calcHessVarVar(gradCOO.dp, variables.x0, hessCOO.dp_dx0, nnzHess, false);
        calcHessVarVar(gradCOO.dp, variables.xf, hessCOO.dp_dxf, nnzHess, false);
        calcHessVarVar(gradCOO.dp, parameters.p, hessCOO.dp_dp, nnzHess, true);
        
        if (std::get<0>(gradCOO.dt0) != NULL) {
            calcHessTimeVar(gradCOO.dt0, variables.x0, hessCOO.dt0_dx0, nnzHess);
            calcHessTimeVar(gradCOO.dt0, variables.xf, hessCOO.dt0_dxf, nnzHess);
            calcHessTimeVar(gradCOO.dt0, parameters.p, hessCOO.dt0_dp, nnzHess);

            if (!variables.t0.is_null() && free_symbols(*std::get<1>(gradCOO.dt0)).count(variables.t0)) {
                auto symDer2 = std::get<1>(gradCOO.dt0)->diff(rcp_dynamic_cast<const Symbol>(variables.t0));
                if (isZero(symDer2)) {
                    hessCOO.dt0_dt0 = std::tuple<double*, RCP<const Basic>>{nullptr, symDer2};
                    nnzHess++;
                }    
            }
        }
        if (std::get<0>(gradCOO.dtf) != NULL) {
            calcHessTimeVar(gradCOO.dtf, variables.x0, hessCOO.dtf_dx0, nnzHess);
            calcHessTimeVar(gradCOO.dtf, variables.xf, hessCOO.dtf_dxf, nnzHess);
            calcHessTimeVar(gradCOO.dtf, parameters.p, hessCOO.dtf_dp, nnzHess);

            if (!variables.tf.is_null() && free_symbols(*std::get<1>(gradCOO.dt0)).count(variables.tf)) {
                auto symDer2 = std::get<1>(gradCOO.dt0)->diff(rcp_dynamic_cast<const Symbol>(variables.tf));
                if (isZero(symDer2)) {
                    hessCOO.dtf_dt0 = std::tuple<double*, RCP<const Basic>>{nullptr, symDer2};
                    nnzHess++;
                }    
            }
        }
        if (!variables.tf.is_null() && free_symbols(*std::get<1>(gradCOO.dtf)).count(variables.tf)) {
            auto symDer2 = std::get<1>(gradCOO.dtf)->diff(rcp_dynamic_cast<const Symbol>(variables.tf));
            if (isZero(symDer2)) {
                hessCOO.dtf_dtf = std::tuple<double*, RCP<const Basic>>{nullptr, symDer2};
                nnzHess++;
            }    
        }
    }

    void setHessPtr(int& offset, double* start) {
        for (auto& coo : hessCOO.dx0_dx0) { std::get<2>(coo) = &start[offset++]; } 
        for (auto& coo : hessCOO.dxf_dx0) { std::get<2>(coo) = &start[offset++]; }
        for (auto& coo : hessCOO.dxf_dxf) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dx0) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dxf) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dp) { std::get<2>(coo) = &start[offset++]; } 

        for (auto& coo : hessCOO.dt0_dx0) { std::get<1>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dt0_dxf) { std::get<1>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dt0_dp) { std::get<1>(coo) = &start[offset++]; }  

        if (!std::get<1>(hessCOO.dt0_dt0).is_null()) {
            std::get<0>(hessCOO.dt0_dt0) = &start[offset++];
        }

        for (auto& coo : hessCOO.dtf_dx0) { std::get<1>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dtf_dxf) { std::get<1>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dtf_dp) { std::get<1>(coo) = &start[offset++]; }  

        if (!std::get<1>(hessCOO.dtf_dt0).is_null()) {
            std::get<0>(hessCOO.dtf_dt0) = &start[offset++];;
        }
        if (!std::get<1>(hessCOO.dtf_dtf).is_null()) {
            std::get<0>(hessCOO.dtf_dtf) = &start[offset++];
        }
    }
};

struct BoundarySweep {
    // M, R :: ensure x0 Size == xf Size
    std::vector<std::unique_ptr<FunctionMR>> functionsMR; 
    int ranges[3]; // {0, dimM, dimM + dimR}
    int xSize, pSize, t0Size, tfSize;
    int timeIndex; // 2*xSize + pSize;

    std::unique_ptr<double[]> inputData; // data input to llvm callbacks, must be filled before every call
    std::unique_ptr<double[]> evalData, gradData, hessData; // data filled from llvm callbacks
    std::unique_ptr<SymEngine::LLVMDoubleVisitor> llvmEval, llvmGrad, llvmHess; // llvm callbacks
    BoundarySweep(const PhaseVariables& variables, const GlobalParameters& parameters, const ModelParameters& mparameters, const PhaseFunctions& functions) {
        initVariableSizes(variables, parameters);
        initFuncIndices(functions);
        initFunctions(functions);
        initDerivatives(variables, parameters);
        compileLLVM(variables, parameters, mparameters);
    };

    void initVariableSizes(const PhaseVariables& variables, const GlobalParameters& parameters) {
        xSize = variables.x0.size();
        pSize = parameters.p.size();
        t0Size = variables.t0.is_null() ? 0 : 1;
        tfSize = variables.tf.is_null() ? 0 : 1;
        timeIndex = 2 * xSize + pSize;
    }

    void initFuncIndices(const PhaseFunctions& functions) {
        ranges[0] = 0;
        ranges[1] = ranges[0] + !functions.M.is_null();
        ranges[2] = ranges[1] + functions.r.size();

        // memory for output of function evaluation
        evalData = std::make_unique<double[]>(ranges[2]);
    }

    void initFunctions(const PhaseFunctions& functions) {
        if (!functions.M.is_null()) {
            functionsMR.push_back(std::make_unique<FunctionMR>(functions.M));
        }
        for (auto& boundary : functions.r) {
            functionsMR.push_back(std::make_unique<FunctionMR>(boundary));
        }
    }

    void initDerivatives(const PhaseVariables& variables, const GlobalParameters& parameters) {  
        int nnzGrad = 0;
        for (auto& function : functionsMR) {
            function->calcGrad(variables, parameters, nnzGrad);
        }
        gradData = std::make_unique<double[]>(nnzGrad);
        int indexGrad = 0;
        for (auto& function : functionsMR) {
            function->setGradPtr(indexGrad, gradData.get());
        }

        int nnzHess = 0;
        for (auto& function : functionsMR) {
            function->calcHess(variables, parameters, nnzHess);
        }
        hessData = std::make_unique<double[]>(nnzHess);
        int indexHess = 0;
        for (auto& function : functionsMR) {
            function->setHessPtr(indexHess, hessData.get());
        }
    }

    void compileLLVM(const PhaseVariables& variables, const GlobalParameters& parameters, const ModelParameters& mparameters) {
        // init variable vector, memory for inputs
        std::vector<RCP<const Basic>> symVariables;
        std::vector<RCP<const Basic>> symEvalVector;
        std::vector<RCP<const Basic>> symGradVector;
        std::vector<RCP<const Basic>> symHessVector;

        for (auto& x0var: variables.x0) {
            symVariables.push_back(x0var);
        }
        for (auto& xfvar: variables.xf) {
            symVariables.push_back(xfvar);
        }
        for (auto& pvar: parameters.p) {
            symVariables.push_back(pvar);
        }
        if (!variables.t0.is_null()) {
            symVariables.push_back(variables.t0);
        }
        if (!variables.tf.is_null()) {
            symVariables.push_back(variables.tf);
        }

        size_t variablesSize = symVariables.size();

        for (auto& mpvar: mparameters.mp) {
            symVariables.push_back(mpvar);
        }
        
        inputData = std::make_unique<double[]>(symVariables.size());
        for (size_t i = 0; i < mparameters.mp.size(); i++) {
            inputData[variablesSize + i] = mparameters.mpValues[i];
        }

        for (auto& function : functionsMR) {
            symEvalVector.push_back(function->symExpr);
            auto vectorGrad = function->gradCOO.toVector();
            auto vectorHess = function->hessCOO.toVector();
            symGradVector.insert(symGradVector.end(), vectorGrad.begin(), vectorGrad.end());
            symHessVector.insert(symHessVector.end(), vectorHess.begin(), vectorHess.end());
        }
        llvmEval = compileVectorExprLLVM(symEvalVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
        llvmGrad = compileVectorExprLLVM(symGradVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
        llvmHess = compileVectorExprLLVM(symHessVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
    }

    inline void fillInputData(double* x0, double* xf, double* p, double* t0, double* tf) {
        // fills the input data array with new values (x0, xf, p, t0, tf, (mp)) - mp values are given as constructed
        for (int i = 0; i < xSize; i++) {
            inputData[i] = x0[i];
            inputData[xSize + i] = xf[i];
        }

        for (int i = 0; i < pSize; i++) {
            inputData[2*xSize + i] = p[i];
        }

        if (t0Size > 0 && tfSize > 0) {
            inputData[timeIndex] = *t0;
            inputData[timeIndex + 1] = *tf;
        }
        else if (t0Size > 0) {
            inputData[timeIndex] = *t0;
        }
        else if (tfSize > 0) {
            inputData[timeIndex] = *tf;
        }
    };

    inline void callEval() {
        llvmEval->call(evalData.get(), inputData.get());
    }

    inline void callGrad() {
        llvmGrad->call(gradData.get(), inputData.get());
    }

    inline void callHess() {
        llvmHess->call(hessData.get(), inputData.get());
    }
};

struct Phase {
    FullSweep fullSweep;
    BoundarySweep boundarySweep;
    PhaseVariables variables;
    PhaseFunctions functions;

    Phase(const PhaseVariables& variables, const GlobalParameters& parameters, const ModelParameters& mparameters, const PhaseFunctions& functions)
    : variables(variables), functions(functions), fullSweep(variables, parameters, mparameters, functions),
      boundarySweep(variables, parameters, mparameters, functions) {
    }
};

struct FunctionE {
    RCP<const Basic> symExpr;
    GradientE gradCOO;
    HessianE hessCOO;

    FunctionE(const RCP<const Basic> expr) : symExpr(expr) {
    }

    void calcGrad(const PhaseVariables& variablesPre, const PhaseVariables& variablesSuc, const GlobalParameters& parameters, int& nnzGrad) {
        calcGradVar(symExpr, variablesPre.xf, gradCOO.dxfPre, nnzGrad);
        calcGradVar(symExpr, variablesSuc.x0, gradCOO.dx0Suc, nnzGrad);
        calcGradVar(symExpr, parameters.p, gradCOO.dp, nnzGrad);
    }

    void setGradPtr(int& offset, double* start) {
        for (auto& coo : gradCOO.dxfPre) { std::get<1>(coo) = &start[offset++]; } 
        for (auto& coo : gradCOO.dx0Suc) { std::get<1>(coo) = &start[offset++]; } 
        for (auto& coo : gradCOO.dp) { std::get<1>(coo) = &start[offset++]; } 
    }

    void calcHess(const PhaseVariables& variablesPre, const PhaseVariables& variablesSuc, const GlobalParameters& parameters, int& nnzHess) {
        calcHessVarVar(gradCOO.dxfPre, variablesPre.xf, hessCOO.dxfPre_dxfPre, nnzHess, true);
        calcHessVarVar(gradCOO.dx0Suc, variablesPre.xf, hessCOO.dx0Suc_dxfPre, nnzHess, false);
        calcHessVarVar(gradCOO.dx0Suc, variablesSuc.x0, hessCOO.dx0Suc_dx0Suc, nnzHess, true);
        calcHessVarVar(gradCOO.dp, variablesPre.xf, hessCOO.dp_dxfPre, nnzHess, false);
        calcHessVarVar(gradCOO.dp, variablesSuc.x0, hessCOO.dp_dx0Suc, nnzHess, false);
        calcHessVarVar(gradCOO.dp, parameters.p, hessCOO.dp_dp, nnzHess, true);
    }

    void setHessPtr(int& offset, double* start) {
        for (auto& coo : hessCOO.dxfPre_dxfPre) { std::get<2>(coo) = &start[offset++]; } 
        for (auto& coo : hessCOO.dx0Suc_dxfPre) { std::get<2>(coo) = &start[offset++]; }
        for (auto& coo : hessCOO.dx0Suc_dx0Suc) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dxfPre) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dx0Suc) { std::get<2>(coo) = &start[offset++]; }  
        for (auto& coo : hessCOO.dp_dp) { std::get<2>(coo) = &start[offset++]; } 
    }
};

struct PhaseConnection {
    // t_{0, suc} - t_{f, pre}
    const double dt0Suc = 1;
    const double dtfPre = -1;
    const int phasePre, phaseSuc;

    PhaseConnection(int phasePre, int PhaseSuc) : phasePre(phasePre), phaseSuc(phaseSuc) {}
 
    inline double eval(const double* tfPre, const double* t0Suc) {
        return *t0Suc - *tfPre;
    }
};

struct Linkage {
    std::vector<std::unique_ptr<FunctionE>> functionsE; // event functions
    std::vector<std::tuple<double, double>> eBounds;    // e^L <= e(x_{f, pre}, x_{0, suc}, p) <= e^U :: event constraints
    std::unique_ptr<PhaseConnection> connection;        // phase connection
    std::tuple<double, double> connBounds;              // t^L <= t_{0, suc} - t_{f, pre}      <= t^U :: phase connections
    int x0SucSize, xfPreSize, pSize;
    int phasePre, phaseSuc;

    // evaluation with compiled code
    std::unique_ptr<double[]> inputData; // data input to llvm callbacks, must be filled before every call
    std::unique_ptr<double[]> evalData, gradData, hessData;  // data filled from llvm callbacks
    std::unique_ptr<SymEngine::LLVMDoubleVisitor> llvmEval, llvmGrad, llvmHess; // llvm callbacks

    Linkage(int phasePre, int phaseSuc, const PhaseVariables& variablesPre = {}, const PhaseVariables& variablesSuc = {}, const GlobalParameters& parameters = {},
            const ModelParameters& mparameters = {}, const std::vector<RCP<const Basic>>& functions = {}, const std::vector<std::tuple<double, double>>& eBounds = {},
            const std::tuple<double, double>& connBounds = {})
            : eBounds(eBounds), connBounds(connBounds), x0SucSize(variablesSuc.x0.size()), xfPreSize(variablesPre.xf.size()), pSize(parameters.p.size()),
            phasePre(phasePre), phaseSuc(phaseSuc), evalData(std::make_unique<double[]>(functions.size())) {
        initFunctions(functions);
        initDerivatives(variablesPre, variablesSuc, parameters);
        compileLLVM(variablesPre, variablesSuc, parameters, mparameters);
    };

    void initFunctions(const std::vector<RCP<const Basic>>& functions) {
        for (auto& event : functions) {
            functionsE.push_back(std::make_unique<FunctionE>(event));
        }
        connection = std::make_unique<PhaseConnection>(phasePre, phaseSuc);
    }

    void initDerivatives(const PhaseVariables& variablesPre, const PhaseVariables& variablesSuc, const GlobalParameters& parameters) {  
        int nnzGrad = 0;
        for (auto& function : functionsE) {
            function->calcGrad(variablesPre, variablesSuc, parameters, nnzGrad);
        }
        gradData = std::make_unique<double[]>(nnzGrad);
        int indexGrad = 0;
        for (auto& function : functionsE) {
            function->setGradPtr(indexGrad, gradData.get());
        }

        int nnzHess = 0;
        for (auto& function : functionsE) {
            function->calcHess(variablesPre, variablesSuc, parameters, nnzHess);
        }
        hessData = std::make_unique<double[]>(nnzHess);
        int indexHess = 0;
        for (auto& function : functionsE) {
            function->setHessPtr(indexHess, hessData.get());
        }
    }

    void compileLLVM(const PhaseVariables& variablesPre, const PhaseVariables& variablesSuc, const GlobalParameters& parameters,
                     const ModelParameters& mparameters) {
        // init variable vector, memory for inputs
        std::vector<RCP<const Basic>> symVariables;
        std::vector<RCP<const Basic>> symEvalVector;
        std::vector<RCP<const Basic>> symGradVector;
        std::vector<RCP<const Basic>> symHessVector;

        for (auto& xfPreVar: variablesPre.xf) {
            symVariables.push_back(xfPreVar);
        }
        for (auto& x0Suc: variablesSuc.x0) {
            symVariables.push_back(x0Suc);
        }
        for (auto& pvar: parameters.p) {
            symVariables.push_back(pvar);
        }

        size_t variablesSize = symVariables.size();

        for (auto& mpvar: mparameters.mp) {
            symVariables.push_back(mpvar);
        }
        
        inputData = std::make_unique<double[]>(symVariables.size());
        for (size_t i = 0; i < mparameters.mp.size(); i++) {
            inputData[variablesSize + i] = mparameters.mpValues[i];
        }
        
        for (auto& function : functionsE) {
            symEvalVector.push_back(function->symExpr);
            auto vectorGrad = function->gradCOO.toVector();
            auto vectorHess = function->hessCOO.toVector();
            symGradVector.insert(symGradVector.end(), vectorGrad.begin(), vectorGrad.end());
            symHessVector.insert(symHessVector.end(), vectorHess.begin(), vectorHess.end());
        }
        llvmEval = compileVectorExprLLVM(symEvalVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
        llvmGrad = compileVectorExprLLVM(symGradVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
        llvmHess = compileVectorExprLLVM(symHessVector, symVariables, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
    }

    inline void fillInputData(double* xfPre, double* x0Suc, double* p) {
        // fills the input data array with new values (xfPre, x0Suc, p, (mp)) - mp values are given as constructed
        for (int i = 0; i < xfPreSize; i++) {
            inputData[i] = xfPre[i];
        }

        for (int i = 0; i < x0SucSize; i++) {
            inputData[xfPreSize + i] = x0Suc[i];
        }

        for (int i = 0; i < pSize; i++) {
            inputData[xfPreSize + x0SucSize + i] = p[i];
        }

    };

    inline void callEval() {
        llvmEval->call(evalData.get(), inputData.get());
    }

    inline void callGrad() {
        llvmGrad->call(gradData.get(), inputData.get());
    }

    inline void callHess() {
        llvmHess->call(hessData.get(), inputData.get());
    }
};

struct Model {
    std::vector<std::unique_ptr<Phase>> phases;
    std::vector<std::unique_ptr<Linkage>> linkages;
    std::unique_ptr<GlobalParameters> parameters;
};

std::vector<RCP<const Basic>> symGrad(const RCP<const Basic>& expr) {
    std::vector<RCP<const Basic>> grad;
    // sort free symbols, look up what type (x, u, p, t) -> adj + symbolic diff -> only take non-zeros
    for (auto var : free_symbols(*expr)) {
        auto symDer = expr->diff(SymEngine::rcp_dynamic_cast<const Symbol>(var));
        grad.push_back(symDer);
    }
    return grad;
}

void evalDiffGDOPT(double *out, const double* p) {
    const double x0 = sin(p[3]);
    const double x1 = exp(x0*p[2]);
    const double x2 = x1/p[4];
    out[0] = 2*p[1]*p[0];
    out[1] = pow(p[0], 2);
    out[2] = x0*x2;
    out[3] = x2*p[2]*cos(p[3]);
    out[4] = -x1/pow(p[4], 2);
}

std::array<std::vector<double>, 6> evalDiff2(const double *x, const double *u, const double *p, double* t)  {
    const double x0 = exp(u[0]);
    const double x1 = pow(u[0], 2);
    const double x2 = 1 + x1;
    const double x3 = pow(x2, -2);
    const double x4 = 2*x3;
    const double x5 = x0/x2;
    const double x6 = x5 - x0*x4*u[0];
    const double x7 = x0*x[0];
    return {std::vector<double>{}, {x6}, {-x4*x7 + x5*x[0] + 8*x1*x7/pow(x2, 3) - 4*x3*x7*u[0]}, {}, {}, {}};
}

FullSweep attitudeControlInit() {
    std::vector<RCP<const Basic>> x =  {symbol("x0"), symbol("x1"), symbol("x2"), symbol("x3"), symbol("x4"), symbol("x5"), symbol("x6"), symbol("x7"), symbol("x8")};
    std::vector<RCP<const Basic>> u =  {symbol("u0"), symbol("u1"), symbol("u2")};
    RCP<const Basic> f1 = parse("4.12885955460403e-08*(-u0 + (593.469491054971*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 0.233455284475168*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(-x4 + x3*x5)*(-133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(0.466910568950337*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 368.600313490232*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - (x1*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) - x2*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7))) - 2.15137087061754e-10*(-u1 + (1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))*(-593.469491054971*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) - 66.5087324606991*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))) + 2*(x3 + x5*x4)*(-0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(133.017464921398*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 108.746865122769*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) - (-x0*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) + x2*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6))) + 9.25440118196439e-09*(-u2 + (133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(x3 + x5*x4)*(-3.7365818746602*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 737.200626980464*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(-217.493730245539*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 3.7365818746602*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) - (x0*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7) - x1*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6)))");
    RCP<const Basic> f2 = parse("-2.15137087061754e-10*(-u0 + (593.469491054971*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 0.233455284475168*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(-x4 + x3*x5)*(-133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(0.466910568950337*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 368.600313490232*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - (x1*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) - x2*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7))) + 1.05114398568516e-08*(-u1 + (1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))*(-593.469491054971*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) - 66.5087324606991*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))) + 2*(x3 + x5*x4)*(-0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(133.017464921398*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 108.746865122769*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) - (-x0*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) + x2*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6))) - 5.64896642555011e-11*(-u2 + (133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(x3 + x5*x4)*(-3.7365818746602*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 737.200626980464*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(-217.493730245539*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 3.7365818746602*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) - (x0*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7) - x1*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6)))");
    RCP<const Basic> f3 = parse("9.25440118196439e-09*(-u0 + (593.469491054971*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 0.233455284475168*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(-x4 + x3*x5)*(-133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(0.466910568950337*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 368.600313490232*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - (x1*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) - x2*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7))) - 5.64896642555011e-11*(-u1 + (1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))*(-593.469491054971*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) - 66.5087324606991*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))) + 2*(x3 + x5*x4)*(-0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(133.017464921398*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 108.746865122769*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) - (-x0*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) + x2*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6))) + 1.513006699675e-08*(-u2 + (133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(x3 + x5*x4)*(-3.7365818746602*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 737.200626980464*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(-217.493730245539*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 3.7365818746602*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) - (x0*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7) - x1*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6)))");
    RCP<const Basic> f4 = parse("0.5*(1 + x3**2)*(x0 + 0.00227276774*(x5 + x3*x4)/(1 + x3**2 + x4**2 + x5**2)) + 0.5*(x1 + 0.00113638387*(1 + 2*(-x3**2 - x5**2)/(1 + x3**2 + x4**2 + x5**2)))*(-x5 + x3*x4) + 0.5*(x4 + x3*x5)*(x2 + 0.00227276774*(-x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2))");
    RCP<const Basic> f5 = parse("0.5*(x1 + 0.00113638387*(1 + 2*(-x3**2 - x5**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + x4**2) + 0.5*(-x3 + x5*x4)*(x2 + 0.00227276774*(-x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2)) + 0.5*(x5 + x3*x4)*(x0 + 0.00227276774*(x5 + x3*x4)/(1 + x3**2 + x4**2 + x5**2))");
    RCP<const Basic> f6 = parse("0.5*(1 + x5**2)*(x2 + 0.00227276774*(-x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2)) + 0.5*(x1 + 0.00113638387*(1 + 2*(-x3**2 - x5**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4) + 0.5*(-x4 + x3*x5)*(x0 + 0.00227276774*(x5 + x3*x4)/(1 + x3**2 + x4**2 + x5**2))");
    RCP<const Basic> f7 = parse("u0");
    RCP<const Basic> f8 = parse("u1");
    RCP<const Basic> f9 = parse("u2");

    RCP<const Basic> g1 = parse("x6**2 + x7**2 + x8**2");

    RCP<const Basic> r1 = parse("-2.15137087061754e-10*((1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))*(-593.469491054971*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) - 66.5087324606991*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))) + 2*(x3 + x5*x4)*(-0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(133.017464921398*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 108.746865122769*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) - (-x0*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) + x2*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6))) + 4.12885955460403e-08*((593.469491054971*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 0.233455284475168*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(-x4 + x3*x5)*(-133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(0.466910568950337*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 368.600313490232*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - (x1*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) - x2*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7))) + 9.25440118196439e-09*((133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(x3 + x5*x4)*(-3.7365818746602*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 737.200626980464*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(-217.493730245539*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 3.7365818746602*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) - (x0*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7) - x1*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6)))");
    RCP<const Basic> r2 = parse("1.05114398568516e-08*((1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))*(-593.469491054971*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) - 66.5087324606991*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))) + 2*(x3 + x5*x4)*(-0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(133.017464921398*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 108.746865122769*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) - (-x0*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) + x2*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6))) - 2.15137087061754e-10*((593.469491054971*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 0.233455284475168*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(-x4 + x3*x5)*(-133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(0.466910568950337*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 368.600313490232*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - (x1*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) - x2*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7))) - 5.64896642555011e-11*((133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(x3 + x5*x4)*(-3.7365818746602*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 737.200626980464*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(-217.493730245539*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 3.7365818746602*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) - (x0*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7) - x1*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6)))");
    RCP<const Basic> r3 = parse("-5.64896642555011e-11*((1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))*(-593.469491054971*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) - 66.5087324606991*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2))) + 2*(x3 + x5*x4)*(-0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(133.017464921398*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2) + 108.746865122769*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) - (-x0*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) + x2*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6))) + 9.25440118196439e-09*((593.469491054971*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 0.233455284475168*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(-x4 + x3*x5)*(-133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 1.8682909373301*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))/(1 + x3**2 + x4**2 + x5**2) + 2*(0.466910568950337*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - 368.600313490232*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) - (x1*(-17167509.4448*x0 + 60260.4448*x1 + 76594401.336*x2 + x8) - x2*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7))) + 1.513006699675e-08*((133.017464921398*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 0.466910568950337*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))*(1 + 2*(-x3**2 - x4**2)/(1 + x3**2 + x4**2 + x5**2)) + 2*(x3 + x5*x4)*(-3.7365818746602*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 737.200626980464*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) + 2*(-x4 + x3*x5)*(-217.493730245539*(x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2) + 3.7365818746602*(-x4 + x3*x5)/(1 + x3**2 + x4**2 + x5**2))/(1 + x3**2 + x4**2 + x5**2) - (x0*(482250.9936*x0 + 95144639.344*x1 + 60260.4448*x2 + x7) - x1*(28070191.1616*x0 + 482250.9936*x1 - 17167509.4448*x2 + x6)))");
    RCP<const Basic> r4 = parse("0.5*(1 + x3**2)*(x0 + 0.00227276774*(x5 + x3*x4)/(1 + x3**2 + x4**2 + x5**2)) + 0.5*(x1 + 0.00113638387*(1 + 2*(-x3**2 - x5**2)/(1 + x3**2 + x4**2 + x5**2)))*(-x5 + x3*x4) + 0.5*(x4 + x3*x5)*(x2 + 0.00227276774*(-x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2)) ");
    RCP<const Basic> r5 = parse("0.5*(x1 + 0.00113638387*(1 + 2*(-x3**2 - x5**2)/(1 + x3**2 + x4**2 + x5**2)))*(1 + x4**2) + 0.5*(-x3 + x5*x4)*(x2 + 0.00227276774*(-x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2)) + 0.5*(x5 + x3*x4)*(x0 + 0.00227276774*(x5 + x3*x4)/(1 + x3**2 + x4**2 + x5**2))");
    RCP<const Basic> r6 = parse("0.5*(1 + x5**2)*(x2 + 0.00227276774*(-x3 + x5*x4)/(1 + x3**2 + x4**2 + x5**2)) + 0.5*(x1 + 0.00113638387*(1 + 2*(-x3**2 - x5**2)/(1 + x3**2 + x4**2 + x5**2)))*(x3 + x5*x4) + 0.5*(-x4 + x3*x5)*(x0 + 0.00227276774*(x5 + x3*x4)/(1 + x3**2 + x4**2 + x5**2)) ");
    RCP<const Basic> r7 = parse("x6");
    RCP<const Basic> r8 = parse("x7");
    RCP<const Basic> r9 = parse("x8");

    RCP<const Basic> L = parse("u0**2 + u1**2 + u2**2");

    auto dynamicEquations = {f1, f2, f3, f4, f5, f6, f7, f8, f9};
    auto pathEquations = {g1};
    auto finalEquations = {r1, r2, r3, r4, r5, r6, r7, r8, r9};
    auto LagrangeTerm = L;

    PhaseVariables variables(x, u, {}, x);
    GlobalParameters parameters;
    ModelParameters mparameters;
    PhaseFunctions functions({}, L, dynamicEquations, pathEquations, {}, finalEquations);

    FullSweep fs(variables, parameters, mparameters, functions);
    BoundarySweep bs(variables, parameters, mparameters, functions);
    return fs;
}

void attitudeControlHessCall(FullSweep& fs, int iterations) {
    double xuvals[12] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    for (int i = 0; i < iterations; i++) {
        fs.fillInputData(xuvals, NULL, NULL);
        fs.callHess();
    }
}

std::array<std::vector<double>, 6> gdoptF1(const double *x, const double *u, const double *p, double t) {
        const double x0 = pow(x[5], 2);
        const double x1 = pow(x[4], 2);
        const double x2 = pow(x[3], 2);
        const double x3 = 1 + x0 + x1 + x2;
        const double x4 = pow(x3, -1);
        const double x5 = 434.987460491077*x4;
        const double x6 = -x5;
        const double x7 = pow(x3, -3);
        const double x8 = x2*x7;
        const double x9 = x[3]*x[5];
        const double x10 = x9 - x[4];
        const double x11 = 1064.13971937119*x10;
        const double x12 = pow(x3, -2);
        const double x13 = x2*x12;
        const double x14 = -x1 - x2;
        const double x15 = x7*x14;
        const double x16 = x2*x15;
        const double x17 = 532.069859685593*x12;
        const double x18 = -x9*x17;
        const double x19 = x12*x10;
        const double x20 = 266.034929842796*x19;
        const double x21 = x14*x12;
        const double x22 = 434.987460491077*x21;
        const double x23 = -x20 - x22;
        const double x24 = x18 + x23;
        const double x25 = 4.30274174123508e-10*x4;
        const double x26 = x25*x10;
        const double x27 = 1.4807041891143e-07*x8;
        const double x28 = x[5]*x[4];
        const double x29 = x28 + x[3];
        const double x30 = 737.200626980464*x4;
        const double x31 = 3.7365818746602*x4;
        const double x32 = x30*x10 - x31*x29;
        const double x33 = x32*x29;
        const double x34 = 266.034929842796*x4;
        const double x35 = 4747.75592843977*x10;
        const double x36 = 2373.87796421989*x12;
        const double x37 = 1186.93898210994*x19;
        const double x38 = 266.034929842796*x21;
        const double x39 = x37 + x38;
        const double x40 = x39 + x9*x36;
        const double x41 = 1 + 2*x4*x14;
        const double x42 = 2.15137087061754e-10*x41;
        const double x43 = 4*x21;
        const double x44 = 4*x4;
        const double x45 = -x43 - x44;
        const double x46 = 16*x13 + 16*x16 + x45;
        const double x47 = 0.466910568950337*x4;
        const double x48 = x47*x10;
        const double x49 = 133.017464921398*x4;
        const double x50 = x49*x29;
        const double x51 = x48 + x50;
        const double x52 = 9.25440118196439e-09*x51;
        const double x53 = 5897.60501584371*x10;
        const double x54 = 29.8926549972816*x29;
        const double x55 = x12*x[3];
        const double x56 = 2948.80250792186*x12;
        const double x57 = -x9*x56;
        const double x58 = 1474.40125396093*x19;
        const double x59 = -x58;
        const double x60 = 7.47316374932039*x12;
        const double x61 = x60*x29;
        const double x62 = x57 + x59 + x61;
        const double x63 = 1.85088023639288e-08*x4;
        const double x64 = x63*x29;
        const double x65 = -x31;
        const double x66 = x60*x[3];
        const double x67 = x65 + x30*x[5] - x58*x[3] + x66*x29;
        const double x68 = x29*x12;
        const double x69 = x68*x[3];
        const double x70 = 266.034929842796*x12;
        const double x71 = x70*x29;
        const double x72 = x71*x[3];
        const double x73 = x47*x[5];
        const double x74 = 0.933821137900673*x19;
        const double x75 = x74*x[3];
        const double x76 = x49 - x72 + x73 - x75;
        const double x77 = -x43*x[3] - x44*x[3];
        const double x78 = x66*x14;
        const double x79 = 7.47316374932039*x4;
        const double x80 = x79*x[3];
        const double x81 = -x73 + x75 - x78 - x80;
        const double x82 = 8.60548348247015e-10*x4;
        const double x83 = 593.469491054971*x4;
        const double x84 = -66.5087324606991*x41 - x83*x10;
        const double x85 = 2.15137087061754e-10*x84;
        const double x86 = 1.8682909373301*x41;
        const double x87 = -x48 + x86;
        const double x88 = x87*x29;
        const double x89 = 3.44219339298806e-09*x88;
        const double x90 = 3.73528455160269*x29;
        const double x91 = 1474.40125396093*x4;
        const double x92 = 5897.60501584371*x15;
        const double x93 = 0.933821137900673*x68;
        const double x94 = 1474.40125396093*x21;
        const double x95 = -x93 + x94;
        const double x96 = 8.25771910920806e-08*x4;
        const double x97 = x96*x29;
        const double x98 = x47 + x91*x[3] - x93*x[3] + x94*x[3];
        const double x99 = -x74;
        const double x100 = 3.73528455160269*x10;
        const double x101 = x8*x100;
        const double x102 = 1.86764227580135*x12;
        const double x103 = x9*x102;
        const double x104 = -x103;
        const double x105 = 1064.13971937119*x29;
        const double x106 = x8*x105;
        const double x107 = -x71;
        const double x108 = 532.069859685593*x55;
        const double x109 = 9.25440118196439e-09*x41;
        const double x110 = -x49;
        const double x111 = x110 + x72 + x78 + x80;
        const double x112 = x19*x[3];
        const double x113 = 1.72109669649403e-09*x87;
        const double x114 = -x79;
        const double x115 = 29.8926549972816*x13;
        const double x116 = 29.8926549972816*x16;
        const double x117 = x60*x14;
        const double x118 = -x117 + x74;
        const double x119 = x103 + x118;
        const double x120 = x25*x29;
        const double x121 = 217.493730245539*x4;
        const double x122 = -x29*x121 + x31*x10;
        const double x123 = x10*x122;
        const double x124 = -x50 - x86;
        const double x125 = x10*x124;
        const double x126 = 6.60617528736645e-07*x125;
        const double x127 = x49*x[5];
        const double x128 = x127 - x20*x[3] - x22*x[3] - x5*x[3];
        const double x129 = -368.600313490232*x41 + x47*x29;
        const double x130 = 3.30308764368322e-07*x129;
        const double x131 = 108.746865122769*x41 + x49*x10;
        const double x132 = x10*x131;
        const double x133 = 3.44219339298806e-09*x132;
        const double x134 = -x121;
        const double x135 = 434.987460491077*x68;
        const double x136 = x31*x[5];
        const double x137 = x134 + x136 + x135*x[3] - x66*x10;
        const double x138 = 3.70176047278576e-08*x4;
        const double x139 = x138*x[5];
        const double x140 = x29*x129;
        const double x141 = 6.60617528736645e-07*x140;
        const double x142 = 4747.75592843977*x29;
        const double x143 = 0.933821137900673*x4;
        const double x144 = 1186.93898210994*x68;
        const double x145 = 0.933821137900673*x21;
        const double x146 = -x144 + x145;
        const double x147 = 4.12885955460403e-08*x41;
        const double x148 = x82*x[5];
        const double x149 = x83 + x143*x[3] - x144*x[3] + x145*x[3];
        const double x150 = 1.65154382184161e-07*x4;
        const double x151 = -0.233455284475168*x41 + x83*x29;
        const double x152 = 4.12885955460403e-08*x151;
        const double x153 = 29.8926549972816*x10;
        const double x154 = 1739.94984196431*x29;
        const double x155 = 14.9463274986408*x12;
        const double x156 = x9*x155;
        const double x157 = -x156;
        const double x158 = -x60*x10;
        const double x159 = x135 + x157 + x158;
        const double x160 = x63*x10;
        const double x161 = x150*x[5];
        const double x162 = x83*x[5];
        const double x163 = -x162 + x34*x[3] + x37*x[3] + x38*x[3];
        const double x164 = x117 + x71;
        const double x165 = x96*x10;
        const double x166 = 7.40352094557151e-08*x32;
        const double x167 = 8.60548348247015e-10*x131;
        const double x168 = 8.60548348247015e-10*x87;
        const double x169 = 1.65154382184161e-07*x124;
        const double x170 = 1.65154382184161e-07*x129;
        const double x171 = 3.70176047278576e-08*x122;
        const double x172 = 3.70176047278576e-08*x32;
        const double x173 = x19*x167 - x19*x169 - x19*x171 + x68*x168 - x68*x170 - x68*x172;
        const double x174 = x9*x12;
        const double x175 = 7.40352094557151e-08*x174;
        const double x176 = 1.72109669649403e-09*x174;
        const double x177 = -x122*x175 - 3.30308764368322e-07*x124*x174 + x176*x131;
        const double x178 = x60*x[4];
        const double x179 = x14*x178;
        const double x180 = x79*x[4];
        const double x181 = x12*x[4];
        const double x182 = 0.933821137900673*x181;
        const double x183 = x10*x182;
        const double x184 = -x179 - x180 + x183 + x47;
        const double x185 = x25*x184;
        const double x186 = 3.70176047278576e-08*x181;
        const double x187 = x29*x186;
        const double x188 = 1474.40125396093*x181;
        const double x189 = x73 + x14*x188 - x29*x182 + x91*x[4];
        const double x190 = x96*x189;
        const double x191 = x63*x67;
        const double x192 = 1186.93898210994*x181;
        const double x193 = x162 + x14*x182 + x143*x[4] - x29*x192;
        const double x194 = 1.65154382184161e-07*x189;
        const double x195 = x[3]*x[4];
        const double x196 = x7*x195;
        const double x197 = x15*x195;
        const double x198 = x181*x[3];
        const double x199 = 16*x197 + 16*x198;
        const double x200 = x10*x186;
        const double x201 = 0.933821137900673*x12;
        const double x202 = -4*x14*x181 - x44*x[4];
        const double x203 = 1739.94984196431*x15;
        const double x204 = x63*x137;
        const double x205 = x28*x12;
        const double x206 = -x136 - x30 - x10*x188 + x61*x[4];
        const double x207 = 3.70176047278576e-08*x206;
        const double x208 = 1.65154382184161e-07*x98;
        const double x209 = x29*x181;
        const double x210 = 434.987460491077*x181;
        const double x211 = 266.034929842796*x181;
        const double x212 = x110 - x10*x211 - x14*x210 - x5*x[4];
        const double x213 = 8.60548348247015e-10*x212;
        const double x214 = 1186.93898210994*x12;
        const double x215 = x10*x181;
        const double x216 = 8.60548348247015e-10*x215;
        const double x217 = 1.65154382184161e-07*x215;
        const double x218 = x96*x111;
        const double x219 = x29*x211;
        const double x220 = -x127 + x179 + x180 + x219;
        const double x221 = 1.65154382184161e-07*x220;
        const double x222 = x63*x206;
        const double x223 = x98*x96;
        const double x224 = x83 + x10*x192 + x14*x211 + x34*x[4];
        const double x225 = 2.15137087061754e-10*x77;
        const double x226 = 1474.40125396093*x12;
        const double x227 = x25*x212;
        const double x228 = 1.4807041891143e-07*x123;
        const double x229 = 2.15137087061754e-10*x202;
        const double x230 = x100*x196;
        const double x231 = 0.933821137900673*x55;
        const double x232 = x28*x201;
        const double x233 = x7*x105;
        const double x234 = x233*x195;
        const double x235 = x9*x70;
        const double x236 = 3.73528455160269*x15;
        const double x237 = 1.4807041891143e-07*x33;
        const double x238 = x7*x154;
        const double x239 = 434.987460491077*x12;
        const double x240 = x81*x25;
        const double x241 = x65 - x10*x178 - x121*x[5] + x29*x210;
        const double x242 = 3.70176047278576e-08*x241;
        const double x243 = x63*x241;
        const double x244 = -x47;
        const double x245 = x127 - x183 - x219 + x244;
        const double x246 = 9.25440118196439e-09*x77;
        const double x247 = 29.8926549972816*x198;
        const double x248 = 29.8926549972816*x197;
        const double x249 = x96*x220;
        const double x250 = 9.25440118196439e-09*x202;
        const double x251 = x25*x128;
        const double x252 = 8.60548348247015e-10*x81;
        const double x253 = 8.60548348247015e-10*x184;
        const double x254 = x1*x7;
        const double x255 = x254*x100;
        const double x256 = 1.86764227580135*x181;
        const double x257 = x1*x12;
        const double x258 = 29.8926549972816*x257;
        const double x259 = x1*x15;
        const double x260 = 29.8926549972816*x259;
        const double x261 = 16*x257 + 16*x259 + x45;
        const double x262 = 1.72109669649403e-09*x181;
        const double x263 = x28*x155;
        const double x264 = 869.974920982154*x12;
        const double x265 = x28*x264;
        const double x266 = 3.30308764368322e-07*x181;
        const double x267 = x254*x105;
        const double x268 = x28*x17;
        const double x269 = x164 + x268;
        const double x270 = 1.4807041891143e-07*x254;
        const double x271 = 7.40352094557151e-08*x181;
        const double x272 = -x268;
        const double x273 = x107 + x272 + x99;
        const double x274 = -x28*x102;
        const double x275 = x274 + x95;
        const double x276 = x146 - x36*x28;
        const double x277 = x173 + x205*x113 - x205*x130 - x205*x166;
        const double x278 = 1.4807041891143e-07*x51;
        const double x279 = x9*x15;
        const double x280 = 3.70176047278576e-08*x21;
        const double x281 = x280*x[5];
        const double x282 = x60*x[5];
        const double x283 = x14*x282;
        const double x284 = x47*x[3];
        const double x285 = x74*x[5];
        const double x286 = -x283 - x284 + x285;
        const double x287 = 8.60548348247015e-10*x286;
        const double x288 = -x20*x[5] - x22*x[5] + x49*x[3];
        const double x289 = x25*x288;
        const double x290 = x7*x9;
        const double x291 = x0*x70;
        const double x292 = -x291 + x49;
        const double x293 = 3.30308764368322e-07*x151;
        const double x294 = 6.60617528736645e-07*x151;
        const double x295 = -x10*x282 - x121*x[4] + x135*x[5] + x31*x[3];
        const double x296 = x295*x[3];
        const double x297 = x9*x233;
        const double x298 = x70*x[5];
        const double x299 = x211*x[3];
        const double x300 = 29.8926549972816*x15;
        const double x301 = x9*x300;
        const double x302 = x19*x[5];
        const double x303 = x25*x286;
        const double x304 = x12*x[5];
        const double x305 = x47*x[4] - x93*x[5] + x94*x[5];
        const double x306 = x96*x305;
        const double x307 = x63*x295;
        const double x308 = x68*x[5];
        const double x309 = x0*x214;
        const double x310 = 1064.13971937119*x15;
        const double x311 = -x144*x[5] + x145*x[5] + x83*x[4];
        const double x312 = 4.12885955460403e-08*x311;
        const double x313 = -x214*x[5];
        const double x314 = x192*x[3];
        const double x315 = 1.65154382184161e-07*x305;
        const double x316 = x0*x12;
        const double x317 = x201*x[5];
        const double x318 = -x317;
        const double x319 = x182*x[3];
        const double x320 = -x319;
        const double x321 = x0*x201;
        const double x322 = x290*x100;
        const double x323 = 0.933821137900673*x13;
        const double x324 = 3.44219339298806e-09*x84;
        const double x325 = 8.60548348247015e-10*x21;
        const double x326 = x325*x[5];
        const double x327 = -x299;
        const double x328 = -x321 + x47;
        const double x329 = x0*x60;
        const double x330 = 1.65154382184161e-07*x21;
        const double x331 = x330*x[5];
        const double x332 = x49*x[4];
        const double x333 = x71*x[5];
        const double x334 = x284 - x285 + x332 - x333;
        const double x335 = x283 - x332 + x333;
        const double x336 = x96*x335;
        const double x337 = x335*x[3];
        const double x338 = x288*x[3];
        const double x339 = x30*x[3] - x31*x[4] - x58*x[5] + x61*x[5];
        const double x340 = x63*x339;
        const double x341 = x37*x[5] + x38*x[5] - x83*x[3];
        const double x342 = x66*x[4];
        const double x343 = x28*x15;
        const double x344 = x7*x28;
        const double x345 = x344*x100;
        const double x346 = x28*x300;
        const double x347 = x344*x105;
        const double x348 = 266.034929842796*x257;
        const double x349 = x0*x15;
        const double x350 = x0*x7;
        const double x351 = x350*x100;
        const double x352 = x350*x105;
        const double x353 = 1.4807041891143e-07*x350;
        const double x354 = x21*x[5];
        const double x355 = 6.60617528736645e-07*x350;
        const double x356 = 29.8926549972816*x349;
        const double x357 = 3.44219339298806e-09*x350;
		return {std::vector<double>{-0.00153915238223185, 0.088075537203803, 0.00394975006480731, 0.00891443466595092, 0.60714200477175, -0.00241059768257546, x173 + x177 + x109*(x101 + x104 + x106 + x107 - x108 + x99) - 3.30308764368322e-07*x111*x112 + x111*x161 + 1.72109669649403e-09*x112*x128 - 7.40352094557151e-08*x112*x137 - x120*(-x101 + x114 + x115 + x116 + x119) - x128*x148 + x137*x139 + x147*(-3.73528455160269*x13 + x143 + x146 - 3.73528455160269*x16 - 2373.87796421989*x55 + x8*x142) + x160*(x159 + 869.974920982154*x55 + x8*x153 - x8*x154) + x165*(-x106 + x108 - x115 - x116 + x164 + x79) - x26*(1739.94984196431*x13 + 1739.94984196431*x16 + x24 + x6 + x8*x11) + x27*x123 + x33*x27 - x42*(-1064.13971937119*x13 - 1064.13971937119*x16 + x34 + x40 - x8*x35) + x46*x152 + x52*x46 + x55*x113 - x55*x130 - x55*x166 + x64*(14.9463274986408*x55 + x62 + x8*x53 - x8*x54) + x67*x138 - 7.40352094557151e-08*x67*x69 - 3.30308764368322e-07*x69*x98 + 8.25771910920806e-08*x77*x149 - 4.30274174123508e-10*x77*x163 + 1.85088023639288e-08*x77*x76 + x8*x126 - x8*x133 + x8*x141 - x8*x89 + 1.72109669649403e-09*x81*x69 - x81*x82 - x85*x46 + x97*(-5897.60501584371*x13 - 1.86764227580135*x55 + x91 + x95 - x2*x92 + x8*x90) + x98*x150, -x185 + x190 - x204 - x218 + x222 + x251 - x120*(-x230 - x231 + x232 + x247 + x248) + x126*x196 + x147*(-x192 - 3.73528455160269*x198 + x196*x142 - x236*x195 - x9*x214) + x160*(x210 + x66 + x196*x153 - x238*x195 - x60*x28 + x9*x239) + x165*(x211 - x234 + x235 - x247 - x248) + x168*x174 + x168*x181 - x170*x174 - x170*x181 - x172*x181 - x174*x172 + x191*x[5] - x196*x133 + x196*x141 + x199*x152 - x200*x137 + 4.12885955460403e-08*x202*x149 + x205*x167 - x205*x169 - x205*x171 - x208*x209 + x209*x252 + x213*x112 + x216*x128 - x217*x111 - x221*x112 + x223*x[5] - x224*x225 - x227*x[5] + x228*x196 - x229*x163 + x237*x196 - x240*x[5] - x242*x112 + x243*x[5] + x245*x246 + x249*x[5] - x26*(1739.94984196431*x198 + x11*x196 + x203*x195 - x70*x28 + x70*x[3]) - x42*(-1064.13971937119*x197 - 1064.13971937119*x198 - x214*x[3] + x28*x214 - x35*x196) + x52*x199 - x55*x167 + x55*x169 + x55*x171 + x64*(x178 + 1474.40125396093*x55 - x28*x226 + x53*x196 - x54*x196 + x9*x60) - x67*x187 - x69*x194 - x69*x207 + x69*x253 + x76*x250 + 4.12885955460403e-08*x77*x193 - x85*x199 - x89*x196 + x97*(-x182 - 5897.60501584371*x198 - x9*x201 + x90*x196 - x92*x195) + (-x211 + x230 + x231 - x232 + x234 - x235)*x109, x277 + x147*(x143 - 3.73528455160269*x257 - 3.73528455160269*x259 + x276 + x254*x142) + x160*(x135 + x158 + 14.9463274986408*x181 + x265 - x1*x238 + x254*x153) + x161*x189 - x184*x148 + 8.25771910920806e-08*x202*x193 - 4.30274174123508e-10*x202*x224 + 1.85088023639288e-08*x202*x245 + x206*x139 - x220*x150 - x241*x138 + x254*x126 - x254*x133 + x254*x141 - x26*(532.069859685593*x181 + x23 + 1739.94984196431*x257 + 1739.94984196431*x259 + x6 + x11*x254) + x261*x152 - x262*x131 + x266*x124 + x270*x123 + x271*x122 + x33*x270 - x42*(-2373.87796421989*x181 - 1064.13971937119*x257 - 1064.13971937119*x259 + x34 + x39 - x35*x254) + x52*x261 + x64*(2948.80250792186*x181 + x263 + x59 + x61 + x53*x254 - x54*x254) + x82*x212 - x85*x261 - x89*x254 + x97*(-5897.60501584371*x257 - 5897.60501584371*x259 + x275 + x91 + x90*x254) + (x255 + x256 + x267 + x273)*x109 + (-x258 - x260 - x267 + x269 + x79)*x165 - (x114 + x118 - x255 - x256 + x258 + x260)*x120 + x10*x212*x262 - x10*x266*x220 - x10*x271*x241 - x29*x206*x271 + x29*x262*x184 - x29*x266*x189, -x303 + x306 + x340 + x13*x167 - x13*x169 - x13*x171 + x147*(x104 + x313 - x314 + x290*x142 - x9*x236) + x160*(x31 - x329 - x2*x60 + x210*x[3] + x239*x[5] + x290*x153 - x9*x238) + x165*(x157 - x297 + x298 + x299 - x301) + x168*x198 - x170*x198 - x172*x198 - 3.70176047278576e-08*x19*x296 - 1.65154382184161e-07*x19*x337 + 8.60548348247015e-10*x19*x338 + x191*x[4] + x204*x[3] - x208*x308 + x218*x[3] + x223*x[4] - x225*x341 + x237*x290 - x240*x[4] + x246*x334 - x25*x131 - x251*x[3] + x252*x308 - x26*(-266.034929842796*x13 + x292 + x11*x290 + x9*x203 + x9*x264) + x279*x278 + x279*x294 - x279*x324 - x289*x[5] + x290*x126 - x290*x133 + x290*x141 + x290*x228 + x293*x174 - 1.65154382184161e-07*x302*x111 + 8.60548348247015e-10*x302*x128 - 3.70176047278576e-08*x302*x137 + x304*x168 - x304*x170 - x304*x172 + x307*x[5] + x316*x167 - x316*x169 - x316*x171 + x326*x163 - x331*x149 + x336*x[5] - x42*(1186.93898210994*x13 + x18 + x309 - x83 - x35*x290 - x9*x310) + x51*x175 + x63*x122 + x64*(-1474.40125396093*x13 + x282 + x30 + x342 - x0*x226 + x53*x290 - x54*x290) - 3.70176047278576e-08*x67*x308 + x69*x287 - x69*x315 - 3.70176047278576e-08*x69*x339 - x76*x281 + x77*x312 - x84*x176 - x89*x290 + x96*x124 + x97*(x318 + x320 + x57 - x9*x92 + x90*x290) - (x156 + x244 + x301 + x321 - x322 + x323)*x120 + (x297 - x298 + x322 - x323 + x327 + x328)*x109, x289 - x307 - x336 - x120*(x263 + x318 + x319 - x345 + x346) + x147*(-1186.93898210994*x257 + x274 - x309 + x83 - x28*x236 + x344*x142) + x160*(x134 + 434.987460491077*x257 + x282 - x342 + x0*x239 - x28*x238 + x344*x153) + x167*x198 - x169*x198 - x171*x198 - x185*x[4] + x190*x[4] - x200*x295 + x202*x312 + x205*x293 - x207*x308 + x209*x287 - x209*x315 + x213*x302 + x216*x288 - x217*x335 - x221*x302 + x222*x[4] + x224*x326 - x227*x[3] + x228*x344 - x229*x341 + x237*x344 - x242*x302 + x243*x[3] + x249*x[3] + x250*x334 + x253*x308 + x257*x168 - x257*x170 - x257*x172 - x26*(x265 + x298 + x327 + x11*x344 + x28*x203) + x278*x343 - x281*x245 + x294*x343 - x303*x[5] - x304*x167 + x304*x169 + x304*x171 + x306*x[5] - x308*x194 + x316*x168 - x316*x170 - x316*x172 - x324*x343 - x331*x193 - x339*x187 + x340*x[5] + x344*x126 - x344*x133 + x344*x141 - x42*(x272 + x313 + x314 - x28*x310 - x35*x344) + 7.40352094557151e-08*x51*x205 + x63*x32 + x64*(-1474.40125396093*x198 + x329 + x65 + x1*x60 + x226*x[5] + x53*x344 - x54*x344) - 1.72109669649403e-09*x84*x205 - x87*x25 - x89*x344 + x96*x129 + x97*(-0.933821137900673*x257 + x328 - x56*x28 + x90*x344 - x92*x28) + (x110 - x263 + x291 - x346 - x347 + x348)*x165 + (x292 + x317 + x320 + x345 + x347 - x348)*x109, x177 + x277 - x120*(x119 - x351 + x356) + x147*(x276 - 3.73528455160269*x349 + x350*x142) + x160*(x159 + x265 - x0*x238 + x350*x153) + x165*(x269 - x352 - x356) - x26*(x24 + 1739.94984196431*x349 + x11*x350) + x278*x349 + 1.72109669649403e-09*x286*x308 + 1.72109669649403e-09*x288*x302 + x294*x349 - 7.40352094557151e-08*x295*x302 + x296*x138 - 3.30308764368322e-07*x302*x335 - 3.30308764368322e-07*x308*x305 - 7.40352094557151e-08*x308*x339 - x324*x349 + x33*x353 - x330*x151 + x337*x150 + x353*x123 - 3.30308764368322e-07*x354*x311 - 7.40352094557151e-08*x354*x334 + 1.72109669649403e-09*x354*x341 + x355*x125 + x355*x140 - x357*x132 - x42*(-1064.13971937119*x349 + x40 - x35*x350) - x51*x280 + x64*(x263 + x62 + x53*x350 - x54*x350) - x82*x338 + x84*x325 - x88*x357 + x97*(x275 - x0*x92 + x90*x350) + (x104 + x273 + x351 + x352)*x109 + x305*x150*x[4] + x339*x138*x[4] - x82*x286*x[4], 9.25440118196439e-09, 2.15137087061754e-10, -9.25440118196439e-09, 4.12885955460403e-08, -2.15137087061754e-10, -4.12885955460403e-08}, {}, {}, {}, {}, {}};
	}

    std::array<std::vector<double>, 6> gdoptF2(const double *x, const double *u, const double *p, double t) {
        const double x0 = pow(x[3], 2);
        const double x1 = pow(x[5], 2);
        const double x2 = pow(x[4], 2);
        const double x3 = 1 + x0 + x1 + x2;
        const double x4 = pow(x3, -3);
        const double x5 = x0*x4;
        const double x6 = x[3]*x[5];
        const double x7 = x6 - x[4];
        const double x8 = 1064.13971937119*x7;
        const double x9 = pow(x3, -1);
        const double x10 = 434.987460491077*x9;
        const double x11 = -x10;
        const double x12 = pow(x3, -2);
        const double x13 = x0*x12;
        const double x14 = -x0 - x2;
        const double x15 = x4*x14;
        const double x16 = 1739.94984196431*x15;
        const double x17 = 532.069859685593*x12;
        const double x18 = -x6*x17;
        const double x19 = x7*x12;
        const double x20 = 266.034929842796*x19;
        const double x21 = 434.987460491077*x12;
        const double x22 = x21*x14;
        const double x23 = -x20 - x22;
        const double x24 = x18 + x23;
        const double x25 = 2.10228797137032e-08*x9;
        const double x26 = x7*x25;
        const double x27 = x[5]*x[4];
        const double x28 = x27 + x[3];
        const double x29 = 737.200626980464*x9;
        const double x30 = 3.7365818746602*x9;
        const double x31 = -x30*x28 + x7*x29;
        const double x32 = x31*x28;
        const double x33 = 9.03834628088017e-10*x32;
        const double x34 = 266.034929842796*x9;
        const double x35 = 4747.75592843977*x5;
        const double x36 = x0*x15;
        const double x37 = 2373.87796421989*x12;
        const double x38 = 1186.93898210994*x19;
        const double x39 = x14*x12;
        const double x40 = 266.034929842796*x39;
        const double x41 = x38 + x40;
        const double x42 = x41 + x6*x37;
        const double x43 = 1 + 2*x9*x14;
        const double x44 = 1.05114398568516e-08*x43;
        const double x45 = 4*x39;
        const double x46 = 4*x9;
        const double x47 = -x45 - x46;
        const double x48 = 16*x13 + 16*x36 + x47;
        const double x49 = 0.466910568950337*x9;
        const double x50 = x7*x49;
        const double x51 = 133.017464921398*x9;
        const double x52 = x51*x28;
        const double x53 = x50 + x52;
        const double x54 = 5.64896642555011e-11*x53;
        const double x55 = 1474.40125396093*x19;
        const double x56 = -x55;
        const double x57 = 5897.60501584371*x7;
        const double x58 = 2948.80250792186*x12;
        const double x59 = -x6*x58;
        const double x60 = 29.8926549972816*x28;
        const double x61 = 7.47316374932039*x12;
        const double x62 = x61*x28;
        const double x63 = x12*x[3];
        const double x64 = 1.12979328511002e-10*x9;
        const double x65 = x64*x28;
        const double x66 = -x30;
        const double x67 = x66 + x29*x[5] - x55*x[3] + x62*x[3];
        const double x68 = x28*x12;
        const double x69 = x68*x[3];
        const double x70 = 266.034929842796*x68;
        const double x71 = x70*x[3];
        const double x72 = x49*x[5];
        const double x73 = 0.933821137900673*x19;
        const double x74 = x73*x[3];
        const double x75 = x51 - x71 + x72 - x74;
        const double x76 = -x45*x[3] - x46*x[3];
        const double x77 = x61*x[3];
        const double x78 = x77*x14;
        const double x79 = 7.47316374932039*x9;
        const double x80 = x79*x[3];
        const double x81 = -x72 + x74 - x78 - x80;
        const double x82 = 4.20457594274063e-08*x9;
        const double x83 = 593.469491054971*x9;
        const double x84 = -66.5087324606991*x43 - x7*x83;
        const double x85 = 1.05114398568516e-08*x84;
        const double x86 = 1.8682909373301*x43;
        const double x87 = -x50 + x86;
        const double x88 = 1.68183037709625e-07*x87;
        const double x89 = x88*x28;
        const double x90 = 3.73528455160269*x28;
        const double x91 = 1.86764227580135*x12;
        const double x92 = 1474.40125396093*x9;
        const double x93 = 5897.60501584371*x15;
        const double x94 = 5897.60501584371*x12;
        const double x95 = 0.933821137900673*x68;
        const double x96 = 1474.40125396093*x12;
        const double x97 = -x95 + x96*x14;
        const double x98 = 4.30274174123508e-10*x9;
        const double x99 = x98*x28;
        const double x100 = x96*x[3];
        const double x101 = x49 + x14*x100 + x92*x[3] - x95*x[3];
        const double x102 = 3.73528455160269*x7;
        const double x103 = x5*x102;
        const double x104 = 1064.13971937119*x28;
        const double x105 = x5*x104;
        const double x106 = 532.069859685593*x63;
        const double x107 = x6*x91;
        const double x108 = -x107;
        const double x109 = -x70 - x73;
        const double x110 = x108 + x109;
        const double x111 = 5.64896642555011e-11*x43;
        const double x112 = -x51;
        const double x113 = x112 + x71 + x78 + x80;
        const double x114 = x19*x[3];
        const double x115 = 8.40915188548126e-08*x87;
        const double x116 = -x79;
        const double x117 = 29.8926549972816*x13;
        const double x118 = 29.8926549972816*x15;
        const double x119 = x0*x118;
        const double x120 = x61*x14;
        const double x121 = -x120 + x73;
        const double x122 = x107 + x121;
        const double x123 = x25*x28;
        const double x124 = 217.493730245539*x9;
        const double x125 = -x28*x124 + x7*x30;
        const double x126 = x7*x125;
        const double x127 = 9.03834628088017e-10*x126;
        const double x128 = -x52 - x86;
        const double x129 = 3.44219339298806e-09*x128;
        const double x130 = x7*x129;
        const double x131 = x51*x[5];
        const double x132 = x131 - x10*x[3] - x20*x[3] - x22*x[3];
        const double x133 = -368.600313490232*x43 + x49*x28;
        const double x134 = 1.72109669649403e-09*x133;
        const double x135 = 108.746865122769*x43 + x7*x51;
        const double x136 = 1.68183037709625e-07*x135;
        const double x137 = x7*x136;
        const double x138 = -x124;
        const double x139 = x21*x28;
        const double x140 = x30*x[5];
        const double x141 = x7*x61;
        const double x142 = x138 + x140 + x139*x[3] - x141*x[3];
        const double x143 = 2.25958657022004e-10*x9;
        const double x144 = x143*x[5];
        const double x145 = 3.44219339298806e-09*x133;
        const double x146 = x28*x145;
        const double x147 = 0.933821137900673*x9;
        const double x148 = 1186.93898210994*x68;
        const double x149 = 0.933821137900673*x39;
        const double x150 = -x148 + x149;
        const double x151 = 2.15137087061754e-10*x43;
        const double x152 = x82*x[5];
        const double x153 = x83 + x147*x[3] - x148*x[3] + x149*x[3];
        const double x154 = 8.60548348247015e-10*x9;
        const double x155 = -0.233455284475168*x43 + x83*x28;
        const double x156 = 2.15137087061754e-10*x155;
        const double x157 = -x141;
        const double x158 = 29.8926549972816*x7;
        const double x159 = 14.9463274986408*x12;
        const double x160 = x6*x159;
        const double x161 = -x160;
        const double x162 = 1739.94984196431*x28;
        const double x163 = x4*x162;
        const double x164 = x7*x64;
        const double x165 = x154*x[5];
        const double x166 = x83*x[5];
        const double x167 = -x166 + x34*x[3] + x38*x[3] + x40*x[3];
        const double x168 = x120 + x70;
        const double x169 = x7*x98;
        const double x170 = 4.51917314044009e-10*x31;
        const double x171 = 4.20457594274063e-08*x135;
        const double x172 = 4.20457594274063e-08*x87;
        const double x173 = 8.60548348247015e-10*x128;
        const double x174 = 8.60548348247015e-10*x133;
        const double x175 = 2.25958657022004e-10*x125;
        const double x176 = 2.25958657022004e-10*x31;
        const double x177 = -x19*x171 + x19*x173 + x19*x175 - x68*x172 + x68*x174 + x68*x176;
        const double x178 = 1.72109669649403e-09*x128;
        const double x179 = x6*x12;
        const double x180 = 4.51917314044009e-10*x125;
        const double x181 = 8.40915188548126e-08*x135;
        const double x182 = x178*x179 + x179*x180 - x179*x181;
        const double x183 = x61*x[4];
        const double x184 = x14*x183;
        const double x185 = x79*x[4];
        const double x186 = x73*x[4];
        const double x187 = -x184 - x185 + x186 + x49;
        const double x188 = x25*x187;
        const double x189 = 2.25958657022004e-10*x67;
        const double x190 = x68*x[4];
        const double x191 = x12*x[4];
        const double x192 = x14*x191;
        const double x193 = 1474.40125396093*x192 + x72 + x92*x[4] - x95*x[4];
        const double x194 = x98*x193;
        const double x195 = x64*x67;
        const double x196 = 0.933821137900673*x191;
        const double x197 = x166 + x14*x196 + x147*x[4] - x148*x[4];
        const double x198 = 2.15137087061754e-10*x76;
        const double x199 = 8.60548348247015e-10*x69;
        const double x200 = x[3]*x[4];
        const double x201 = x4*x200;
        const double x202 = x15*x200;
        const double x203 = x191*x[3];
        const double x204 = 16*x202 + 16*x203;
        const double x205 = 2.25958657022004e-10*x19;
        const double x206 = x205*x[4];
        const double x207 = x4*x90;
        const double x208 = 0.933821137900673*x12;
        const double x209 = -4*x192 - x46*x[4];
        const double x210 = 2.15137087061754e-10*x209;
        const double x211 = 266.034929842796*x12;
        const double x212 = x4*x8;
        const double x213 = x64*x142;
        const double x214 = x27*x12;
        const double x215 = -x140 - x29 + x28*x183 - x55*x[4];
        const double x216 = 2.25958657022004e-10*x215;
        const double x217 = 8.60548348247015e-10*x101;
        const double x218 = x112 - x10*x[4] - x20*x[4] - x22*x[4];
        const double x219 = 4.20457594274063e-08*x218;
        const double x220 = 1186.93898210994*x12;
        const double x221 = 4747.75592843977*x7;
        const double x222 = x4*x221;
        const double x223 = 1064.13971937119*x15;
        const double x224 = x19*x[4];
        const double x225 = 4.20457594274063e-08*x132;
        const double x226 = 8.60548348247015e-10*x113;
        const double x227 = x98*x113;
        const double x228 = x70*x[4];
        const double x229 = -x131 + x184 + x185 + x228;
        const double x230 = 8.60548348247015e-10*x229;
        const double x231 = x64*x215;
        const double x232 = x98*x101;
        const double x233 = 266.034929842796*x191;
        const double x234 = x83 + x14*x233 + x34*x[4] + x38*x[4];
        const double x235 = 1.05114398568516e-08*x76;
        const double x236 = x4*x60;
        const double x237 = x25*x218;
        const double x238 = 1.05114398568516e-08*x209;
        const double x239 = x201*x102;
        const double x240 = 0.933821137900673*x63;
        const double x241 = x27*x208;
        const double x242 = x201*x104;
        const double x243 = x6*x211;
        const double x244 = 1186.93898210994*x191;
        const double x245 = 4747.75592843977*x28;
        const double x246 = x4*x245;
        const double x247 = x21*x[4];
        const double x248 = x81*x25;
        const double x249 = x66 - x124*x[5] + x28*x247 - x7*x183;
        const double x250 = 2.25958657022004e-10*x249;
        const double x251 = x64*x249;
        const double x252 = -x49;
        const double x253 = x131 - x186 - x228 + x252;
        const double x254 = 29.8926549972816*x203;
        const double x255 = x200*x118;
        const double x256 = x98*x229;
        const double x257 = x25*x132;
        const double x258 = 4.20457594274063e-08*x68;
        const double x259 = x258*x[3];
        const double x260 = x2*x4;
        const double x261 = x260*x102;
        const double x262 = x91*x[4];
        const double x263 = x2*x12;
        const double x264 = 29.8926549972816*x263;
        const double x265 = x2*x118;
        const double x266 = x2*x15;
        const double x267 = 16*x263 + 16*x266 + x47;
        const double x268 = x27*x159;
        const double x269 = x268 + x56 + x62;
        const double x270 = 869.974920982154*x12;
        const double x271 = x27*x270;
        const double x272 = x139 + x157 + x271;
        const double x273 = x7*x260;
        const double x274 = 4747.75592843977*x260;
        const double x275 = x260*x104;
        const double x276 = x27*x17;
        const double x277 = x168 + x276;
        const double x278 = 9.03834628088017e-10*x260;
        const double x279 = -x276;
        const double x280 = -x91*x27;
        const double x281 = x280 + x97;
        const double x282 = x150 - x37*x27;
        const double x283 = x177 - x214*x115 + x214*x134 + x214*x170;
        const double x284 = 4.51917314044009e-10*x53;
        const double x285 = 9.03834628088017e-10*x53;
        const double x286 = x6*x15;
        const double x287 = 2.25958657022004e-10*x39;
        const double x288 = x287*x[5];
        const double x289 = x61*x[5];
        const double x290 = x14*x289;
        const double x291 = x49*x[3];
        const double x292 = x73*x[5];
        const double x293 = -x290 - x291 + x292;
        const double x294 = -x20*x[5] - x22*x[5] + x51*x[3];
        const double x295 = x25*x294;
        const double x296 = x4*x6;
        const double x297 = x1*x211;
        const double x298 = -x297 + x51;
        const double x299 = 1.72109669649403e-09*x155;
        const double x300 = 3.44219339298806e-09*x155;
        const double x301 = x21*x[5];
        const double x302 = -x124*x[4] + x28*x301 + x30*x[3] - x7*x289;
        const double x303 = x302*x[3];
        const double x304 = x296*x104;
        const double x305 = x211*x[5];
        const double x306 = x233*x[3];
        const double x307 = x6*x118;
        const double x308 = x19*x[5];
        const double x309 = x25*x293;
        const double x310 = x12*x[5];
        const double x311 = x96*x[5];
        const double x312 = x14*x311 + x49*x[4] - x95*x[5];
        const double x313 = x98*x312;
        const double x314 = x64*x302;
        const double x315 = x68*x[5];
        const double x316 = 4.20457594274063e-08*x315;
        const double x317 = x1*x220;
        const double x318 = -x148*x[5] + x149*x[5] + x83*x[4];
        const double x319 = -x220*x[5];
        const double x320 = x244*x[3];
        const double x321 = 3.73528455160269*x15;
        const double x322 = x1*x12;
        const double x323 = x208*x[5];
        const double x324 = -x323;
        const double x325 = x196*x[3];
        const double x326 = -x325;
        const double x327 = x1*x208;
        const double x328 = x296*x102;
        const double x329 = 0.933821137900673*x13;
        const double x330 = 1.68183037709625e-07*x84;
        const double x331 = 4.20457594274063e-08*x39;
        const double x332 = x331*x[5];
        const double x333 = -x306;
        const double x334 = -x327 + x49;
        const double x335 = x1*x61;
        const double x336 = 8.40915188548126e-08*x84;
        const double x337 = 8.60548348247015e-10*x39;
        const double x338 = x337*x[5];
        const double x339 = x51*x[4];
        const double x340 = x70*x[5];
        const double x341 = x291 - x292 + x339 - x340;
        const double x342 = 5.64896642555011e-11*x341;
        const double x343 = x290 - x339 + x340;
        const double x344 = x98*x343;
        const double x345 = x343*x[3];
        const double x346 = x294*x[3];
        const double x347 = x28*x289 + x29*x[3] - x30*x[4] - x55*x[5];
        const double x348 = x64*x347;
        const double x349 = x38*x[5] + x40*x[5] - x83*x[3];
        const double x350 = x183*x[3];
        const double x351 = x27*x15;
        const double x352 = x4*x27;
        const double x353 = x352*x102;
        const double x354 = x27*x118;
        const double x355 = x352*x104;
        const double x356 = 266.034929842796*x263;
        const double x357 = x293*x[4];
        const double x358 = x312*x[4];
        const double x359 = x347*x[4];
        const double x360 = x1*x15;
        const double x361 = x1*x4;
        const double x362 = x361*x102;
        const double x363 = x361*x104;
        const double x364 = 4747.75592843977*x361;
        const double x365 = x39*x[5];
        const double x366 = x28*x361;
        const double x367 = x1*x118;
		return {std::vector<double>{-0.360856001648463, 0.000729069125149687, -2.85558803120564e-05, 0.50995897086773, -0.00809020963983677, 0.360884557528775, x177 + x182 - x101*x154 + 1.72109669649403e-09*x113*x114 - x113*x165 - 8.40915188548126e-08*x114*x132 + 4.51917314044009e-10*x114*x142 + x123*(-x103 + x116 + x117 + x119 + x122) + x132*x152 - x142*x144 - x151*(-3.73528455160269*x13 + x147 + x150 - 3.73528455160269*x36 - 2373.87796421989*x63 + x35*x28) - x164*(x139 + x157 + x161 + 869.974920982154*x63 - x0*x163 + x5*x158) - x169*(-x105 + x106 - x117 - x119 + x168 + x79) + x26*(x11 + 1739.94984196431*x13 + x24 + x0*x16 + x5*x8) + x44*(-1064.13971937119*x13 + x34 - 1064.13971937119*x36 + x42 - x7*x35) - x48*x156 - x5*x127 - x5*x130 + x5*x137 - x5*x146 - x5*x33 + x5*x89 - x54*x48 - x63*x115 + x63*x134 + x63*x170 - x65*(x56 + x59 + x62 + 14.9463274986408*x63 + x5*x57 - x5*x60) - x67*x143 + 4.51917314044009e-10*x67*x69 + 1.72109669649403e-09*x69*x101 - 1.12979328511002e-10*x75*x76 - 4.30274174123508e-10*x76*x153 + 2.10228797137032e-08*x76*x167 - 8.40915188548126e-08*x81*x69 + x81*x82 + x85*x48 - x99*(x92 + x97 - x0*x93 - x0*x94 + x5*x90 - x91*x[3]) - (x103 + x105 - x106 + x110)*x111, x188 - x194 + x213 + x227 - x231 - x257 + x123*(-x239 - x240 + x241 + x254 + x255) - x151*(-3.73528455160269*x202 - 3.73528455160269*x203 - x244 + x200*x246 - x6*x220) - x164*(x247 + x77 - x200*x163 + x201*x158 + x6*x21 - x61*x27) - x169*(x233 - x242 + x243 - x254 - x255) - x172*x191 + x174*x191 + x176*x191 - x179*x172 + x179*x174 + x179*x176 + x189*x190 - x195*x[5] - x197*x198 + x199*x193 - x201*x127 - x201*x130 + x201*x137 - x201*x146 - x204*x156 + x206*x142 - x210*x153 - x214*x171 + x214*x173 + x214*x175 + x217*x190 - x219*x114 - x224*x225 + x224*x226 + x230*x114 - x232*x[5] + x234*x235 + x237*x[5] + x238*x167 + x248*x[5] + x250*x114 - x251*x[5] - x256*x[5] - x259*x187 + x26*(1739.94984196431*x203 + x16*x200 + x211*x[3] + x212*x200 - x27*x211) - x33*x201 + x44*(-1064.13971937119*x203 - x200*x222 - x200*x223 - x220*x[3] + x27*x220) - x54*x204 + x63*x171 - x63*x173 - x63*x175 - x65*(x100 + x183 - x236*x200 + x57*x201 + x6*x61 - x96*x27) + x69*x216 - 5.64896642555011e-11*x75*x209 - 5.64896642555011e-11*x76*x253 + x85*x204 + x89*x201 - x99*(-x196 + x200*x207 - x6*x208 - x93*x200 - x94*x200) - (-x233 + x239 + x240 - x241 + x242 - x243)*x111 - x81*x258*x[4], x283 - x111*(x109 + x261 + x262 + x275 + x279) - x151*(x147 - 3.73528455160269*x263 - 3.73528455160269*x266 + x282 + x28*x274) - x164*(14.9463274986408*x191 + x272 + x260*x158 - x260*x162) - x165*x193 - x178*x191 - x180*x191 + x181*x191 + x187*x152 - 8.40915188548126e-08*x187*x190 + 1.72109669649403e-09*x190*x193 - 4.30274174123508e-10*x209*x197 - 1.12979328511002e-10*x209*x253 - x215*x144 + 4.51917314044009e-10*x215*x190 - 8.40915188548126e-08*x218*x224 + 1.72109669649403e-09*x224*x229 + 4.51917314044009e-10*x224*x249 + x229*x154 + 2.10228797137032e-08*x234*x209 + x249*x143 + x26*(x11 + 532.069859685593*x191 + x23 + 1739.94984196431*x263 + 1739.94984196431*x266 + x8*x260) - x260*x146 - x267*x156 - x273*x129 + x273*x136 - x278*x126 - x32*x278 + x44*(-2373.87796421989*x191 - 1064.13971937119*x263 - 1064.13971937119*x266 + x34 + x41 - x7*x274) - x54*x267 - x65*(2948.80250792186*x191 + x269 + x57*x260 - x60*x260) - x82*x218 + x85*x267 + x89*x260 - x99*(x281 + x92 - x2*x93 - x2*x94 + x90*x260) - (-x264 - x265 - x275 + x277 + x79)*x169 + (x116 + x121 - x261 - x262 + x264 + x265)*x123, x309 - x313 - x348 - x13*x171 + x13*x173 + x13*x175 - x151*(x108 + x319 - x320 + x296*x245 - x6*x321) - x164*(x30 + x301 - x335 - x0*x61 + x247*x[3] + x296*x158 - x6*x163) - x169*(x161 - x304 + x305 + x306 - x307) + 8.60548348247015e-10*x19*x345 - 4.20457594274063e-08*x19*x346 - x195*x[4] - x203*x172 + x203*x174 + x203*x176 + x205*x303 - x213*x[3] + x217*x315 - x225*x308 + x226*x308 - x227*x[3] - x232*x[4] + x235*x349 + x248*x[4] + x25*x135 + x257*x[3] + x26*(x298 - x0*x211 + x6*x16 + x6*x270 + x8*x296) - x284*x179 - x286*x285 - x286*x300 + x286*x330 - x293*x259 + x295*x[5] - x296*x127 - x296*x130 + x296*x137 - x296*x146 - x299*x179 + 2.25958657022004e-10*x308*x142 - x310*x172 + x310*x174 + x310*x176 + x312*x199 - x314*x[5] + x315*x189 - x318*x198 - x322*x171 + x322*x173 + x322*x175 - x33*x296 - x332*x167 + x336*x179 + x338*x153 - x344*x[5] + x44*(x18 + x317 - x83 + x0*x220 - x296*x221 - x6*x223) - x64*x125 - x65*(x289 + x29 + x350 - x0*x96 - x1*x96 + x57*x296 - x6*x236) + 2.25958657022004e-10*x69*x347 + x75*x288 - x76*x342 - x81*x316 + x89*x296 - x98*x128 - x99*(x324 + x326 + x59 - x6*x93 + x90*x296) + (x160 + x252 + x307 + x327 - x328 + x329)*x123 - (x304 - x305 + x328 - x329 + x333 + x334)*x111, -x295 + x314 + x344 + x123*(x268 + x324 + x325 - x353 + x354) - x151*(x280 - x317 + x83 - x2*x220 + x27*x246 - x27*x321) - x164*(x138 + x289 - x350 + x1*x21 + x2*x21 - x27*x163 + x352*x158) + x188*x[4] - x194*x[4] - x203*x171 + x203*x173 + x203*x175 + x206*x302 - x209*x342 - x210*x318 - x214*x284 - x214*x299 + x214*x336 + x216*x315 - x219*x308 + 8.60548348247015e-10*x224*x343 + x230*x308 - x231*x[4] - x234*x332 + x237*x[3] + x238*x349 + x250*x308 - x251*x[3] - x256*x[3] - x258*x357 + x26*(x271 + x305 + x333 + x27*x16 + x27*x212) - x263*x172 + x263*x174 + x263*x176 - x285*x351 + x288*x253 - 4.20457594274063e-08*x294*x224 + x309*x[5] + x310*x171 - x310*x173 - x310*x175 - x313*x[5] + 8.60548348247015e-10*x315*x193 - x316*x187 - x322*x172 + x322*x174 + x322*x176 - x33*x352 + x338*x197 - x348*x[5] - x351*x300 + x351*x330 - x352*x127 - x352*x130 + x352*x137 - x352*x146 + x44*(x279 + x319 + x320 - x27*x222 - x27*x223) - x64*x31 - x65*(-1474.40125396093*x203 + x311 + x335 + x66 + x2*x61 - x27*x236 + x57*x352) + 8.60548348247015e-10*x68*x358 + 2.25958657022004e-10*x68*x359 + x87*x25 + x89*x352 - x98*x133 - x99*(-0.933821137900673*x263 + x334 + x27*x207 - x58*x27 - x93*x27) - (x112 - x268 + x297 - x354 - x355 + x356)*x169 - (x298 + x323 + x326 + x353 + x355 - x356)*x111, x182 + x283 + x123*(x122 - x362 + x367) - x151*(x282 - 3.73528455160269*x360 + x28*x364) - x164*(x161 + x272 - x1*x163 + x361*x158) - x169*(x277 - x363 - x367) + x26*(x24 + x1*x16 + x8*x361) - x285*x360 - 8.40915188548126e-08*x293*x315 - 8.40915188548126e-08*x294*x308 - x303*x143 + 4.51917314044009e-10*x308*x302 + 1.72109669649403e-09*x308*x343 + 1.72109669649403e-09*x312*x315 + 4.51917314044009e-10*x315*x347 - x33*x361 + x337*x155 - x345*x154 - x358*x154 - x359*x143 - x360*x300 + x360*x330 - x361*x127 - x361*x130 + x361*x137 + 1.72109669649403e-09*x365*x318 + 4.51917314044009e-10*x365*x341 - 8.40915188548126e-08*x365*x349 - x366*x145 + x44*(-1064.13971937119*x360 + x42 - x7*x364) + x53*x287 - x65*(x269 + x59 + x57*x361 - x60*x361) + x82*x346 + x82*x357 - x84*x331 + x88*x366 - x99*(x281 - x1*x93 + x90*x361) - (x110 + x279 + x362 + x363)*x111, -5.64896642555011e-11, -1.05114398568516e-08, 5.64896642555011e-11, -2.15137087061754e-10, 1.05114398568516e-08, 2.15137087061754e-10}, {}, {}, {}, {}, {}};
	}

    std::array<std::vector<double>, 6> gdoptF3(const double *x, const double *u, const double *p, double t) {
        const double x0 = pow(x[5], 2);
        const double x1 = pow(x[4], 2);
        const double x2 = pow(x[3], 2);
        const double x3 = 1 + x0 + x1 + x2;
        const double x4 = pow(x3, -1);
        const double x5 = 434.987460491077*x4;
        const double x6 = -x5;
        const double x7 = x[3]*x[5];
        const double x8 = x7 - x[4];
        const double x9 = pow(x3, -3);
        const double x10 = x8*x9;
        const double x11 = 1064.13971937119*x10;
        const double x12 = pow(x3, -2);
        const double x13 = x2*x12;
        const double x14 = -x1 - x2;
        const double x15 = x9*x14;
        const double x16 = 1739.94984196431*x2;
        const double x17 = 532.069859685593*x12;
        const double x18 = -x7*x17;
        const double x19 = 266.034929842796*x12;
        const double x20 = 434.987460491077*x12;
        const double x21 = x20*x14;
        const double x22 = -x21 - x8*x19;
        const double x23 = x18 + x22;
        const double x24 = 1.12979328511002e-10*x4;
        const double x25 = x8*x24;
        const double x26 = 2.42081071948e-07*x2;
        const double x27 = 737.200626980464*x4;
        const double x28 = x[5]*x[4];
        const double x29 = x28 + x[3];
        const double x30 = 3.7365818746602*x4;
        const double x31 = -x30*x29 + x8*x27;
        const double x32 = x9*x29;
        const double x33 = x32*x31;
        const double x34 = 266.034929842796*x4;
        const double x35 = 4747.75592843977*x2;
        const double x36 = 1064.13971937119*x15;
        const double x37 = 2373.87796421989*x12;
        const double x38 = 1186.93898210994*x12;
        const double x39 = x14*x19 + x8*x38;
        const double x40 = x39 + x7*x37;
        const double x41 = 1 + 2*x4*x14;
        const double x42 = 5.64896642555011e-11*x41;
        const double x43 = 16*x15;
        const double x44 = x14*x12;
        const double x45 = 4*x4;
        const double x46 = -4*x44 - x45;
        const double x47 = 16*x13 + x46 + x2*x43;
        const double x48 = 0.466910568950337*x4;
        const double x49 = x8*x48;
        const double x50 = 133.017464921398*x4;
        const double x51 = x50*x29;
        const double x52 = x49 + x51;
        const double x53 = 1.513006699675e-08*x52;
        const double x54 = 5897.60501584371*x10;
        const double x55 = 29.8926549972816*x32;
        const double x56 = x12*x[3];
        const double x57 = 2948.80250792186*x12;
        const double x58 = -x7*x57;
        const double x59 = 1474.40125396093*x12;
        const double x60 = x8*x59;
        const double x61 = -x60;
        const double x62 = 7.47316374932039*x12;
        const double x63 = x62*x29;
        const double x64 = x58 + x61 + x63;
        const double x65 = 3.02601339935e-08*x4;
        const double x66 = x65*x29;
        const double x67 = -x30;
        const double x68 = 7.47316374932039*x56;
        const double x69 = 1474.40125396093*x56;
        const double x70 = x67 + x27*x[5] + x68*x29 - x8*x69;
        const double x71 = 1.21040535974e-07*x56;
        const double x72 = 266.034929842796*x56;
        const double x73 = x72*x29;
        const double x74 = x48*x[5];
        const double x75 = 0.933821137900673*x56;
        const double x76 = x8*x75;
        const double x77 = x50 - x73 + x74 - x76;
        const double x78 = 4*x14;
        const double x79 = -x45*x[3] - x78*x56;
        const double x80 = x68*x14;
        const double x81 = 7.47316374932039*x4;
        const double x82 = x81*x[3];
        const double x83 = -x74 + x76 - x80 - x82;
        const double x84 = 2.25958657022004e-10*x4;
        const double x85 = 593.469491054971*x4;
        const double x86 = -66.5087324606991*x41 - x8*x85;
        const double x87 = 5.64896642555011e-11*x86;
        const double x88 = 4.51917314044009e-10*x56;
        const double x89 = 1.8682909373301*x41;
        const double x90 = -x49 + x89;
        const double x91 = 9.03834628088017e-10*x90*x32;
        const double x92 = 3.73528455160269*x32;
        const double x93 = 1474.40125396093*x4;
        const double x94 = 5897.60501584371*x15;
        const double x95 = 0.933821137900673*x12;
        const double x96 = x95*x29;
        const double x97 = x59*x14;
        const double x98 = -x96 + x97;
        const double x99 = 1.85088023639288e-08*x4;
        const double x100 = x99*x29;
        const double x101 = 7.40352094557151e-08*x56;
        const double x102 = x48 + x69*x14 - x75*x29 + x93*x[3];
        const double x103 = x29*x102;
        const double x104 = x8*x95;
        const double x105 = -x104;
        const double x106 = 3.73528455160269*x10;
        const double x107 = x2*x106;
        const double x108 = 1.86764227580135*x12;
        const double x109 = x7*x108;
        const double x110 = -x109;
        const double x111 = 1064.13971937119*x32;
        const double x112 = x2*x111;
        const double x113 = x29*x19;
        const double x114 = -x113;
        const double x115 = 532.069859685593*x56;
        const double x116 = 1.513006699675e-08*x41;
        const double x117 = -x50;
        const double x118 = x117 + x73 + x80 + x82;
        const double x119 = x8*x118;
        const double x120 = -x81;
        const double x121 = 29.8926549972816*x13;
        const double x122 = 29.8926549972816*x15;
        const double x123 = x2*x122;
        const double x124 = 7.47316374932039*x44;
        const double x125 = x104 - x124;
        const double x126 = x109 + x125;
        const double x127 = x24*x29;
        const double x128 = 217.493730245539*x4;
        const double x129 = -x29*x128 + x8*x30;
        const double x130 = x8*x129;
        const double x131 = x9*x130;
        const double x132 = 1.4807041891143e-07*x2;
        const double x133 = -x51 - x89;
        const double x134 = x10*x133;
        const double x135 = x50*x[5];
        const double x136 = x135 - x21*x[3] - x5*x[3] - x8*x72;
        const double x137 = -368.600313490232*x41 + x48*x29;
        const double x138 = 108.746865122769*x41 + x8*x50;
        const double x139 = 9.03834628088017e-10*x138;
        const double x140 = x10*x139;
        const double x141 = -x128;
        const double x142 = x20*x29;
        const double x143 = x30*x[5];
        const double x144 = x141 + x143 + x142*x[3] - x8*x68;
        const double x145 = 6.05202679869999e-08*x4;
        const double x146 = x145*x[5];
        const double x147 = x32*x137;
        const double x148 = 0.933821137900673*x4;
        const double x149 = 3.73528455160269*x15;
        const double x150 = 0.933821137900673*x44;
        const double x151 = x150 - x38*x29;
        const double x152 = 9.25440118196439e-09*x41;
        const double x153 = x84*x[5];
        const double x154 = x38*x[3];
        const double x155 = x85 + x148*x[3] - x29*x154 + x75*x14;
        const double x156 = 3.70176047278576e-08*x4;
        const double x157 = -0.233455284475168*x41 + x85*x29;
        const double x158 = 9.25440118196439e-09*x157;
        const double x159 = 29.8926549972816*x10;
        const double x160 = 14.9463274986408*x12;
        const double x161 = x7*x160;
        const double x162 = -x161;
        const double x163 = -x8*x62;
        const double x164 = x142 + x162 + x163;
        const double x165 = x8*x65;
        const double x166 = x156*x[5];
        const double x167 = x85*x[5];
        const double x168 = -x167 + x34*x[3] + x72*x14 + x8*x154;
        const double x169 = x113 + x124;
        const double x170 = x8*x99;
        const double x171 = x8*x144;
        const double x172 = x7*x12;
        const double x173 = 7.40352094557151e-08*x172;
        const double x174 = 1.21040535974e-07*x172;
        const double x175 = 4.51917314044009e-10*x138;
        const double x176 = 2.25958657022004e-10*x90;
        const double x177 = x29*x12;
        const double x178 = 3.70176047278576e-08*x137;
        const double x179 = x8*x12;
        const double x180 = 3.70176047278576e-08*x133;
        const double x181 = 6.05202679869999e-08*x12;
        const double x182 = x31*x181;
        const double x183 = 2.25958657022004e-10*x138;
        const double x184 = x176*x177 - x178*x177 - x179*x180 + x179*x183 - x181*x130 - x29*x182;
        const double x185 = x184 - x129*x174 - x173*x133 + x175*x172;
        const double x186 = x12*x[4];
        const double x187 = 7.47316374932039*x186;
        const double x188 = x14*x187;
        const double x189 = x81*x[4];
        const double x190 = 0.933821137900673*x186;
        const double x191 = x8*x190;
        const double x192 = -x188 - x189 + x191 + x48;
        const double x193 = x24*x192;
        const double x194 = x29*x186;
        const double x195 = x74 - x29*x190 + x93*x[4] + x97*x[4];
        const double x196 = x99*x195;
        const double x197 = x70*x65;
        const double x198 = x38*x[4];
        const double x199 = x167 + x14*x190 + x148*x[4] - x29*x198;
        const double x200 = 9.25440118196439e-09*x79;
        const double x201 = x56*x29;
        const double x202 = 3.70176047278576e-08*x195;
        const double x203 = x[3]*x[4];
        const double x204 = 1.4807041891143e-07*x134;
        const double x205 = x56*x[4];
        const double x206 = 16*x205 + x43*x203;
        const double x207 = 6.05202679869999e-08*x186;
        const double x208 = -x45*x[4] - x78*x186;
        const double x209 = 9.25440118196439e-09*x208;
        const double x210 = 1739.94984196431*x15;
        const double x211 = x65*x144;
        const double x212 = x129*x181;
        const double x213 = 1.4807041891143e-07*x147;
        const double x214 = -x143 - x27 + x29*x187 - x60*x[4];
        const double x215 = 3.70176047278576e-08*x186;
        const double x216 = x20*x[4];
        const double x217 = x19*x[4];
        const double x218 = x117 - x14*x216 - x5*x[4] - x8*x217;
        const double x219 = 2.25958657022004e-10*x8;
        const double x220 = x56*x219;
        const double x221 = 4747.75592843977*x10;
        const double x222 = x219*x186;
        const double x223 = x99*x118;
        const double x224 = x113*x[4];
        const double x225 = -x135 + x188 + x189 + x224;
        const double x226 = 3.70176047278576e-08*x225;
        const double x227 = x8*x56;
        const double x228 = x65*x214;
        const double x229 = x99*x102;
        const double x230 = x85 + x14*x217 + x34*x[4] + x8*x198;
        const double x231 = x24*x218;
        const double x232 = x56*x183;
        const double x233 = 2.42081071948e-07*x131;
        const double x234 = x28*x12;
        const double x235 = 6.05202679869999e-08*x129;
        const double x236 = x203*x106;
        const double x237 = x95*x28;
        const double x238 = x203*x111;
        const double x239 = x7*x19;
        const double x240 = 4747.75592843977*x32;
        const double x241 = 2.42081071948e-07*x33;
        const double x242 = 1739.94984196431*x32;
        const double x243 = x83*x24;
        const double x244 = x67 - x128*x[5] + x29*x216 - x8*x187;
        const double x245 = x8*x244;
        const double x246 = 6.05202679869999e-08*x56;
        const double x247 = x65*x244;
        const double x248 = -x48;
        const double x249 = x135 - x191 - x224 + x248;
        const double x250 = 1.513006699675e-08*x79;
        const double x251 = 29.8926549972816*x205;
        const double x252 = x203*x122;
        const double x253 = x99*x225;
        const double x254 = 6.05202679869999e-08*x31;
        const double x255 = 1.513006699675e-08*x208;
        const double x256 = x24*x136;
        const double x257 = 2.25958657022004e-10*x83;
        const double x258 = 2.25958657022004e-10*x192;
        const double x259 = x1*x106;
        const double x260 = 1.86764227580135*x186;
        const double x261 = x1*x12;
        const double x262 = 29.8926549972816*x261;
        const double x263 = x1*x122;
        const double x264 = 16*x261 + x46 + x1*x43;
        const double x265 = x28*x160;
        const double x266 = 869.974920982154*x12;
        const double x267 = x28*x266;
        const double x268 = 1.4807041891143e-07*x1;
        const double x269 = 4747.75592843977*x1;
        const double x270 = 7.40352094557151e-08*x186;
        const double x271 = x1*x111;
        const double x272 = x28*x17;
        const double x273 = x169 + x272;
        const double x274 = 1.21040535974e-07*x186;
        const double x275 = -x272;
        const double x276 = x105 + x114 + x275;
        const double x277 = -x28*x108;
        const double x278 = x277 + x98;
        const double x279 = x151 - x37*x28;
        const double x280 = 7.40352094557151e-08*x234;
        const double x281 = 1.21040535974e-07*x234;
        const double x282 = -x280*x137 - x31*x281 + 4.51917314044009e-10*x90*x234;
        const double x283 = 2.42081071948e-07*x52;
        const double x284 = x7*x15;
        const double x285 = 6.05202679869999e-08*x44;
        const double x286 = x285*x[5];
        const double x287 = x124*x[5];
        const double x288 = x48*x[3];
        const double x289 = x104*x[5];
        const double x290 = -x287 - x288 + x289;
        const double x291 = 2.25958657022004e-10*x290;
        const double x292 = x20*x[5];
        const double x293 = x19*x[5];
        const double x294 = -x14*x292 + x50*x[3] - x8*x293;
        const double x295 = x24*x294;
        const double x296 = x0*x19;
        const double x297 = -x296 + x50;
        const double x298 = 1.4807041891143e-07*x157;
        const double x299 = x62*x[5];
        const double x300 = -x128*x[4] + x29*x292 + x30*x[3] - x8*x299;
        const double x301 = x8*x300;
        const double x302 = x7*x111;
        const double x303 = x72*x[4];
        const double x304 = x7*x122;
        const double x305 = x12*x[5];
        const double x306 = x219*x305;
        const double x307 = x24*x290;
        const double x308 = x59*x[5];
        const double x309 = x14*x308 + x48*x[4] - x96*x[5];
        const double x310 = x99*x309;
        const double x311 = x65*x300;
        const double x312 = 3.70176047278576e-08*x305;
        const double x313 = x29*x305;
        const double x314 = x0*x38;
        const double x315 = x38*x[5];
        const double x316 = x150*x[5] - x29*x315 + x85*x[4];
        const double x317 = -x315;
        const double x318 = x198*x[3];
        const double x319 = 3.70176047278576e-08*x309;
        const double x320 = x0*x12;
        const double x321 = x95*x[5];
        const double x322 = -x321;
        const double x323 = x75*x[4];
        const double x324 = -x323;
        const double x325 = x0*x95;
        const double x326 = x7*x106;
        const double x327 = 0.933821137900673*x13;
        const double x328 = 9.03834628088017e-10*x86;
        const double x329 = x15*x328;
        const double x330 = 2.25958657022004e-10*x44;
        const double x331 = x330*x[5];
        const double x332 = -x303;
        const double x333 = -x325 + x48;
        const double x334 = x0*x62;
        const double x335 = 4.51917314044009e-10*x86;
        const double x336 = x181*x[5];
        const double x337 = x29*x336;
        const double x338 = 3.70176047278576e-08*x44;
        const double x339 = x338*x[5];
        const double x340 = x50*x[4];
        const double x341 = x113*x[5];
        const double x342 = x288 - x289 + x340 - x341;
        const double x343 = x287 - x340 + x341;
        const double x344 = x99*x343;
        const double x345 = x27*x[3] + x29*x299 - x30*x[4] - x60*x[5];
        const double x346 = 6.05202679869999e-08*x345;
        const double x347 = x65*x345;
        const double x348 = x14*x293 + x8*x315 - x85*x[3];
        const double x349 = 5.64896642555011e-11*x348;
        const double x350 = x68*x[4];
        const double x351 = x28*x15;
        const double x352 = x28*x106;
        const double x353 = x28*x122;
        const double x354 = x28*x111;
        const double x355 = 266.034929842796*x261;
        const double x356 = x8*x305;
        const double x357 = x0*x15;
        const double x358 = x0*x10;
        const double x359 = 3.73528455160269*x358;
        const double x360 = x0*x111;
        const double x361 = x44*x[5];
        const double x362 = x0*x122;
		return {std::vector<double>{-0.0126534059955547, -0.8559692791647, 0.0134776310216689, 0.000810083257082163, -0.0880469813234908, -0.000824225026114151, x185 + x100*(-5897.60501584371*x13 - 1.86764227580135*x56 + x93 + x98 + x2*x92 - x2*x94) - x101*x103 - x101*x119 - x101*x137 + x102*x156 + x118*x166 - x127*(-x107 + x120 + x121 + x123 + x126) + x132*x147 + x134*x132 - x136*x153 + x144*x146 + x152*(-3.73528455160269*x13 + x148 + x151 - 2373.87796421989*x56 - x2*x149 + x32*x35) + x165*(x164 + 869.974920982154*x56 + x2*x159 - x32*x16) + x170*(-x112 + x115 - x121 - x123 + x169 + x81) - x2*x140 - x2*x91 - x25*(1739.94984196431*x13 + x23 + x6 + x15*x16 + x2*x11) + x26*x131 + x33*x26 - x42*(-1064.13971937119*x13 + x34 + x40 - x2*x36 - x35*x10) + x47*x158 + x53*x47 + x66*(14.9463274986408*x56 + x64 + x2*x54 - x2*x55) + x70*x145 - x71*x171 - x71*x31 + 1.85088023639288e-08*x79*x155 - 1.12979328511002e-10*x79*x168 + 3.02601339935e-08*x79*x77 - x83*x84 - x87*x47 + x88*x90 + (x105 + x107 + x110 + x112 + x114 - x115)*x116 - x71*x70*x29 + x8*x88*x136 + x83*x88*x29, -x193 + x196 - x211 - x223 + x228 - x232 + x256 + x100*(-x190 - 5897.60501584371*x205 - x7*x95 + x92*x203 - x94*x203) + x116*(-x217 + x236 - x237 + x238 - x239 + x75) + x152*(-x198 - 3.73528455160269*x205 - x203*x149 + x203*x240 - x7*x38) + x165*(x216 + x68 + x203*x159 - x203*x242 - x62*x28 + x7*x20) + x170*(x217 - x238 + x239 - x251 - x252) + x176*x172 + x176*x186 - x178*x172 - x178*x186 + x197*x[5] + x200*x199 - x201*x202 + x201*x258 - x203*x140 + x203*x204 + x203*x241 + x206*x158 - x207*x171 - 5.64896642555011e-11*x208*x168 + x209*x155 + x213*x203 - 6.05202679869999e-08*x214*x201 - x215*x103 - x215*x119 + x218*x220 + x222*x136 - x227*x226 + x229*x[5] - x231*x[5] + x233*x203 - x234*x180 + x234*x183 - x243*x[5] - x245*x246 + x247*x[5] - x25*(1739.94984196431*x205 + x72 + x11*x203 + x210*x203 - x28*x19) + x250*x249 + x253*x[5] - x254*x186 + x257*x194 - x28*x212 - x42*(-x154 - 1064.13971937119*x205 - x203*x221 - x36*x203 + x38*x28) + x53*x206 + x56*x180 + x56*x235 + x66*(x187 + x69 + x54*x203 - x55*x203 - x59*x28 + x7*x62) - x7*x182 - 6.05202679869999e-08*x70*x194 + x77*x255 - 5.64896642555011e-11*x79*x230 - x87*x206 - x91*x203 - (-x236 + x237 + x251 + x252 - x75)*x127, x184 + x282 - x1*x140 + x1*x233 + x1*x241 - x1*x91 + x100*(-5897.60501584371*x261 + x278 + x93 + x1*x92 - x1*x94) + x152*(x148 - 3.73528455160269*x261 + x279 - x1*x149 + x32*x269) + x165*(x142 + x163 + 14.9463274986408*x186 + x267 + x1*x159 - x1*x242) + x166*x195 - x175*x186 - x192*x153 + 4.51917314044009e-10*x192*x194 + 1.85088023639288e-08*x208*x199 + 3.02601339935e-08*x208*x249 + x214*x146 - x225*x156 - 1.12979328511002e-10*x230*x208 - x244*x145 - x25*(532.069859685593*x186 + x22 + 1739.94984196431*x261 + x6 + x1*x11 + x1*x210) + x264*x158 + x268*x134 + x268*x147 + x270*x133 + x274*x129 - x274*x245 - x42*(-2373.87796421989*x186 - 1064.13971937119*x261 + x34 + x39 - x1*x36 - x10*x269) + x53*x264 + x66*(2948.80250792186*x186 + x265 + x61 + x63 + x1*x54 - x1*x55) + x84*x218 - x87*x264 + (x259 + x260 + x271 + x276)*x116 + (-x262 - x263 - x271 + x273 + x81)*x170 - (x120 + x125 - x259 - x260 + x262 + x263)*x127 - x29*x214*x274 - x29*x270*x195 + 4.51917314044009e-10*x8*x218*x186 - x8*x270*x225, -x307 + x310 + x347 - x0*x212 - x13*x180 + x13*x183 - x13*x235 + x152*(x110 + x317 - x318 - x7*x149 + x7*x240) + x165*(-7.47316374932039*x13 + x292 + x30 - x334 + x216*x[3] + x7*x159 - x7*x242) + x170*(x162 + x293 - x302 + x303 - x304) + x173*x157 - x182*x[5] + x197*x[4] + x200*x316 + x201*x291 - x201*x319 - x201*x346 + x205*x176 - x205*x178 - x205*x254 + x211*x[3] + x223*x[3] - 3.70176047278576e-08*x227*x343 + x229*x[4] - x24*x138 - x243*x[4] - x246*x301 - x25*(x297 - x2*x19 + x7*x11 + x7*x210 + x7*x266) + x250*x342 - x256*x[3] + x257*x313 + x283*x284 + x294*x220 - x295*x[5] + x298*x284 + x305*x176 - x305*x178 + x306*x136 + x311*x[5] - x312*x103 - x312*x119 - x320*x180 + x320*x183 + x331*x168 - x335*x172 - x336*x171 - x339*x155 + x344*x[5] - x42*(x18 + x314 - x85 + x2*x38 - x7*x221 - x7*x36) + x52*x174 + x65*x129 + x66*(x27 + x299 + x350 - x0*x59 - x2*x59 + x7*x54 - x7*x55) - x7*x140 + x7*x204 + x7*x213 + x7*x233 + x7*x241 - x7*x329 - x7*x91 - x70*x337 - x77*x286 - x79*x349 + x99*x133 + (x322 + x324 + x58 + x7*x92 - x7*x94)*x100 - (x161 + x248 + x304 + x325 - x326 + x327)*x127 + (-x293 + x302 + x326 - x327 + x332 + x333)*x116, x295 - x311 - x344 - x0*x182 + x100*(-0.933821137900673*x261 + x333 - x57*x28 + x92*x28 - x94*x28) - x127*(x265 + x322 + x323 - x352 + x353) + x152*(-1186.93898210994*x261 + x277 - x314 + x85 - x28*x149 + x28*x240) + x165*(x141 + 434.987460491077*x261 + x299 - x350 + x0*x20 + x28*x159 - x28*x242) - x193*x[4] + x196*x[4] - x202*x313 - x205*x180 - x207*x301 - x208*x349 + x209*x316 + x212*x[5] - x214*x337 + x218*x306 - x226*x356 + x228*x[4] + x230*x331 - x231*x[3] + x232*x[4] - x234*x335 - x235*x205 - x245*x336 + x247*x[3] - x25*(x267 + x293 + x332 + x28*x11 + x28*x210) + x253*x[3] + x255*x342 + x258*x313 + x261*x176 - x261*x178 - x261*x254 - x28*x140 + x28*x204 + x28*x213 + x28*x233 + x28*x241 - x28*x329 + x280*x157 + x283*x351 - x286*x249 + x291*x194 + x294*x222 + x298*x351 + x305*x180 - x305*x183 - x307*x[5] + x310*x[5] - x319*x194 + x320*x176 - x320*x178 - x339*x199 - x346*x194 + x347*x[5] - x42*(x275 + x317 + x318 - x28*x221 - x36*x28) + x52*x281 + x65*x31 + x66*(7.47316374932039*x261 + x308 + x334 + x67 + x54*x28 - x55*x28 - x69*x[4]) - x90*x24 - x91*x28 + x99*x137 + (x117 - x265 + x296 - x353 - x354 + x355)*x170 + (x297 + x321 + x324 + x352 + x354 - x355)*x116 - x8*x215*x343, x185 + x282 + x0*x213 + x0*x233 + x0*x241 - x0*x91 + x100*(x278 + x0*x92 - x0*x94) - x127*(x126 - x359 + x362) + x152*(x279 - x0*x149 + x0*x240) + x165*(x164 + x267 + 29.8926549972816*x358 - x0*x242) + x170*(x273 - x360 - x362) - x25*(x23 + 1064.13971937119*x358 + x0*x210) + x283*x357 + 4.51917314044009e-10*x290*x313 + 4.51917314044009e-10*x294*x356 + x298*x357 - 7.40352094557151e-08*x309*x313 - 1.21040535974e-07*x313*x345 - x338*x157 - 1.21040535974e-07*x356*x300 - 7.40352094557151e-08*x356*x343 - x357*x328 + 1.4807041891143e-07*x358*x133 - x358*x139 - 7.40352094557151e-08*x361*x316 - 1.21040535974e-07*x361*x342 + 4.51917314044009e-10*x361*x348 - x42*(-4747.75592843977*x358 + x40 - x0*x36) - x52*x285 + x66*(x265 + 5897.60501584371*x358 + x64 - x0*x55) + x86*x330 + (x110 + x276 + x359 + x360)*x116 + x300*x145*x[3] + x309*x156*x[4] + x343*x156*x[3] + x345*x145*x[4] - x84*x290*x[4] - x84*x294*x[3], 1.513006699675e-08, 5.64896642555011e-11, -1.513006699675e-08, 9.25440118196439e-09, -5.64896642555011e-11, -9.25440118196439e-09}, {}, {}, {}, {}, {}};
	}

    std::array<std::vector<double>, 6> gdoptF4(const double *x, const double *u, const double *p, double t)  {
        const double x0 = 1.0*x[3];
        const double x1 = 0.5*x[5];
        const double x2 = pow(x[5], 2);
        const double x3 = pow(x[4], 2);
        const double x4 = pow(x[3], 2);
        const double x5 = 1 + x4;
        const double x6 = x2 + x3 + x5;
        const double x7 = pow(x6, -1);
        const double x8 = 0.00227276774*x7;
        const double x9 = x[5]*x[4];
        const double x10 = x9 - x[3];
        const double x11 = pow(x6, -2);
        const double x12 = 0.00454553548*x11;
        const double x13 = x12*x[3];
        const double x14 = -x8 - x13*x10;
        const double x15 = -x2 - x4;
        const double x16 = x15*x12;
        const double x17 = 0.00454553548*x7;
        const double x18 = -x16*x[3] - x17*x[3];
        const double x19 = 0.01818214192*x11;
        const double x20 = pow(x6, -3);
        const double x21 = 0.01818214192*x20;
        const double x22 = x4*x21;
        const double x23 = -x16 - x17;
        const double x24 = x[3]*x[4];
        const double x25 = x24 - x[5];
        const double x26 = 0.5*x25;
        const double x27 = x24 + x[5];
        const double x28 = x21*x27;
        const double x29 = -x27*x12;
        const double x30 = 0.00909107096*x11;
        const double x31 = x29 - x30*x24;
        const double x32 = 0.5*x5;
        const double x33 = -x12*x10;
        const double x34 = x[3]*x[5];
        const double x35 = 0.5*(x34 + x[4]);
        const double x36 = x8*x[4];
        const double x37 = 0.5*x[3];
        const double x38 = x12*x[4];
        const double x39 = x21*x10;
        const double x40 = -x34*x12;
        const double x41 = 0.00227276774*x11;
        const double x42 = x41*x10;
        const double x43 = x8*x[5];
        const double x44 = x43 - x38*x10;
        const double x45 = 0.00113638387*x7;
        const double x46 = x24*x12;
        const double x47 = 0.00909107096*x25*x20*x15;
        const double x48 = x41*x15;
        const double x49 = x8*x[3];
        const double x50 = x8 - x3*x12;
        const double x51 = x9*x30;
        const double x52 = x33 - x51;
        const double x53 = x12*x[5];
        const double x54 = x36 - x53*x10;
        const double x55 = 0.5*(-x16*x[5] - x17*x[5]);
        const double x56 = x21*x15;
        const double x57 = x2*x21;
		return {std::vector<double>{x0, 0.5*x[4], x1, 1.0*x[0] + 1.0*x14*x[5] + 1.0*x18*x[4] + x26*(x23 + x22*x15 + x4*x19) + x32*(x31 + x4*x28) + x35*(x33 + x22*x10 + x30*x[3]) + x8*x27 + 2.0*x[3]*(x36 - x27*x13), x37, 0.5, -x45 + 0.5*x[1] + x0*(x49 - x38*x27) + x1*x44 - x3*x48 + x32*(x50 + x24*x28 - x4*x12) + x35*(x38 + x40 + x39*x24) + x37*x18 - x42*x[3] + x46*x25 + x47*x24 + 0.000568191935*(1 + 2*x7*x15), x44 - x24*x16 + x3*x47 + x32*(x31 + x3*x28) + x35*(x52 + x3*x39) - x48*x25, -0.5, x37, x49 + 0.5*x[2] + x0*(x8 - x53*x27) + x1*x54 + x32*(-x13 + x34*x28 - x9*x12) + x35*(-x46 + x53 + x34*x39) + x37*x14 + x45*x10 + x48*x[3] + x55*x[4] + (x34*x19 + x56*x34)*x26, x26*(x51 + x9*x56) + x32*(-x38 + x40 + x9*x28) + x35*(x50 - x2*x12 + x9*x39) - x42*x[5] + x44*x37 + x45*x[4] + x48*x[4] + x55*x[3], x43 - x55 + x0*x54 + x26*(x23 + x2*x19 + x57*x15) + x32*(x29 + x2*x28 - x30*x[5]) + x35*(x52 + x57*x10) + x48*x[5]}, {}, {}, {}, {}, {}};
	}

    std::array<std::vector<double>, 6> gdoptF5(const double *x, const double *u, const double *p, double t) {
        const double x0 = 0.5*x[4];
        const double x1 = pow(x[5], 2);
        const double x2 = pow(x[3], 2);
        const double x3 = pow(x[4], 2);
        const double x4 = 1 + x3;
        const double x5 = x1 + x2 + x4;
        const double x6 = pow(x5, -1);
        const double x7 = 0.00227276774*x6;
        const double x8 = x[5]*x[4];
        const double x9 = x8 - x[3];
        const double x10 = pow(x5, -2);
        const double x11 = 0.00454553548*x10;
        const double x12 = x11*x[3];
        const double x13 = -x7 - x9*x12;
        const double x14 = 0.00227276774*x10;
        const double x15 = 0.01818214192*x10;
        const double x16 = -x1 - x2;
        const double x17 = pow(x5, -3);
        const double x18 = 0.01818214192*x17;
        const double x19 = x18*x16;
        const double x20 = 0.00454553548*x6;
        const double x21 = x11*x16;
        const double x22 = -x20 - x21;
        const double x23 = 0.5*x4;
        const double x24 = 0.00113638387*x6;
        const double x25 = -x9*x11;
        const double x26 = x9*x18;
        const double x27 = 0.00909107096*x10;
        const double x28 = 0.5*x9;
        const double x29 = x[3]*x[4];
        const double x30 = x29 + x[5];
        const double x31 = x30*x18;
        const double x32 = -x30*x11;
        const double x33 = x32 - x29*x27;
        const double x34 = 0.5*x30;
        const double x35 = x7*x[4];
        const double x36 = x35 - x30*x12;
        const double x37 = 1.0*x[4];
        const double x38 = 0.5*x[3];
        const double x39 = 0.5*x[5];
        const double x40 = x11*x[4];
        const double x41 = -x12*x[5];
        const double x42 = x4*x16;
        const double x43 = 0.00909107096*x17;
        const double x44 = -x40*x30 + x7*x[3];
        const double x45 = x7 - x3*x11;
        const double x46 = x29*x11;
        const double x47 = x7*x[5] - x9*x40;
        const double x48 = 0.5*x47;
        const double x49 = x8*x27;
        const double x50 = x25 - x49;
        const double x51 = x3*x16;
        const double x52 = x11*x[5];
        const double x53 = x7 - x52*x30;
        const double x54 = x35 - x9*x52;
        const double x55 = 0.5*x54;
        const double x56 = x[3]*x[5];
        const double x57 = x30*x14;
		return {std::vector<double>{x0, -0.5, -0.5*x13 + x24 + x23*(x22 + x2*x15 + x2*x19) + x28*(x25 + x2*x26 + x27*x[3]) + x34*(x33 + x2*x31) + x36*x37 + x9*x14*x[3], x38, x37, x39, -x48 + 0.5*x[0] + x0*x44 + x28*(x40 + x41 + x29*x26) + x30*x24 + x34*(x45 - x2*x11 + x31*x29) + x36*x38 + x39*x13 + x4*x46 + (-x20*x[3] - x21*x[3])*x37 + x42*x43*x29, 1.0*x[1] + x28*(x50 + x3*x26) + x34*(x33 + x3*x31) - x42*x14 + 1.0*x44*x[3] + 1.0*x47*x[5] - x51*x27 + x4*x51*x43 + 0.00113638387*(1 + 2*x6*x16), 0.5, x0, -x55 + x0*x13 + x0*x53 + x24*x[4] + x28*(-x46 + x52 + x56*x26) + x34*(-x12 + x56*x31 - x8*x11) - x57*x[3] + (x56*x15 + x56*x19)*x23, 0.5*x[2] + x23*(x49 + x8*x19) + x24*x[3] + x28*(x45 - x1*x11 + x8*x26) + x34*(-x40 + x41 + x8*x31) + x48*x[4] + x53*x38 + x55*x[5] - x57*x[4] + x9*x24 + (-x20*x[5] - x21*x[5])*x37, x53 + x23*(x22 + x1*x15 + x1*x19) + x28*(x50 + x1*x26) + x34*(x32 + x1*x31 - x27*x[5]) + x54*x37}, {}, {}, {}, {}, {}};
	}

    std::array<std::vector<double>, 6> gdoptF6(const double *x, const double *u, const double *p, double t)  {
        const double x0 = 0.5*x[5];
        const double x1 = pow(x[3], 2);
        const double x2 = pow(x[4], 2);
        const double x3 = pow(x[5], 2);
        const double x4 = 1 + x3;
        const double x5 = x1 + x2 + x4;
        const double x6 = pow(x5, -2);
        const double x7 = 0.01818214192*x6;
        const double x8 = -x1 - x3;
        const double x9 = pow(x5, -3);
        const double x10 = 0.01818214192*x9;
        const double x11 = x1*x10;
        const double x12 = pow(x5, -1);
        const double x13 = 0.00454553548*x12;
        const double x14 = 0.00454553548*x6;
        const double x15 = x8*x14;
        const double x16 = -x13 - x15;
        const double x17 = x[5]*x[4];
        const double x18 = x17 + x[3];
        const double x19 = 0.5*x18;
        const double x20 = 0.00227276774*x12;
        const double x21 = x20*x[4];
        const double x22 = x[3]*x[4];
        const double x23 = x22 + x[5];
        const double x24 = x14*x[3];
        const double x25 = x21 - x24*x23;
        const double x26 = 1.0*x[5];
        const double x27 = x17 - x[3];
        const double x28 = -x27*x14;
        const double x29 = 0.00909107096*x6;
        const double x30 = 0.5*x4;
        const double x31 = -x23*x14;
        const double x32 = x31 - x22*x29;
        const double x33 = x[3]*x[5];
        const double x34 = 0.5*(x33 - x[4]);
        const double x35 = -x13*x[3] - x15*x[3];
        const double x36 = 0.00227276774*x6;
        const double x37 = x8*x36;
        const double x38 = 0.00909107096*x8*x9*x18;
        const double x39 = x14*x[4];
        const double x40 = 0.5*(x20*x[3] - x39*x23);
        const double x41 = x23*x10;
        const double x42 = x20 - x2*x14;
        const double x43 = x36*x23;
        const double x44 = x27*x10;
        const double x45 = -x33*x14;
        const double x46 = x22*x14;
        const double x47 = 0.00113638387*x12;
        const double x48 = x29*x17;
        const double x49 = x28 - x48;
        const double x50 = 0.5*x[3];
        const double x51 = 0.5*x[4];
        const double x52 = x14*x[5];
        const double x53 = x20 - x52*x23;
        const double x54 = x20*x[5];
        const double x55 = x8*x10;
        const double x56 = -x13*x[5] - x15*x[5];
        const double x57 = x3*x10;
		return {std::vector<double>{x0, 0.5, x35 + x19*(x16 + x1*x7 + x8*x11) + x25*x26 + x30*(x28 + x27*x11 + x29*x[3]) + x34*(x32 + x23*x11), -0.5, x0, x0*x35 + x30*(x39 + x45 + x44*x22) + x34*(x42 - x1*x14 + x41*x22) - x37*x[4] + x38*x22 + x40*x[5] + x43*x[3] + x46*x18 - x47*x[4], -x40 - x15*x17 + x2*x38 + x30*(x49 + x2*x44) + x34*(x32 + x2*x41) - x37*x18 + x43*x[4] - x47*x[3], x50, x51, x26, -x54 + 0.5*x[0] + x0*x53 + x19*(x55*x33 + x7*x33) + x26*(-x20 - x24*x27) + x30*(-x46 + x52 + x44*x33) + x34*(-x24 - x14*x17 + x41*x33) - x37*x[5] + x47*x23 + x50*x25 + x51*x35, -0.5*x53 + 0.5*x[1] + x0*x56 + x19*(x48 + x55*x17) - x2*x37 + x26*(x54 - x39*x27) + x30*(x42 - x3*x14 + x44*x17) + x34*(-x39 + x45 + x41*x17) + x40*x[3] + 0.000568191935*(1 + 2*x8*x12), 1.0*x[2] + x19*(x16 + x3*x7 + x8*x57) + x20*x27 + x30*(x49 + x3*x44) + x34*(x31 - x29*x[5] + x57*x23) + 1.0*x53*x[3] + 1.0*x56*x[4] + 2.0*x[5]*(x21 - x52*x27)}, {}, {}, {}, {}, {}};
	}

    std::array<std::vector<double>, 6> gdoptF7(const double *x, const double *u, const double *p, double t) {
		return {std::vector<double>{}, {}, {}, {}, {}, {}};
	}

void attitudeControlGDOPTHessCall(int iterations) {
    double xvals[9] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    double uvals[3] = {10, 11, 12};
    for (int i = 0; i < iterations; i++) {
        gdoptF1(xvals, uvals, NULL, 0);
        gdoptF2(xvals, uvals, NULL, 0);
        gdoptF3(xvals, uvals, NULL, 0);
        gdoptF4(xvals, uvals, NULL, 0);
        gdoptF5(xvals, uvals, NULL, 0);
        gdoptF6(xvals, uvals, NULL, 0);
        gdoptF7(xvals, uvals, NULL, 0);
        gdoptF7(xvals, uvals, NULL, 0);
        gdoptF7(xvals, uvals, NULL, 0);
    }
}

int main() {
    using namespace SymEngine;
    /*
    std::vector<RCP<const Basic>> X =  {symbol("x1"), symbol("x2"), symbol("x3"), symbol("x4"), symbol("x5")};
    RCP<const Basic> f = add(mul(mul(X[0], X[0]), X[1]), div(exp(mul(X[2], sin(X[3]))), X[4]));
    std::vector<double> xvals = {1, 1, 5, 1, 1};
    RCP<const Basic> t = symbol("t");
    std::cout << *f << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto grad = symGrad(f);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "[SYMENGINE]  Time for gradient calculation : " << duration << " microseconds\n";
         
    start = std::chrono::high_resolution_clock::now();
    auto llvmGrad = compileVectorExprLLVM(grad, X, SYM_LLVM_USE_CSE, SYM_LLVM_OPT_LEVEL);
    end = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "[SYMENGINE]  Time for llvm compile gradient: " << duration << " microseconds\n";

    int evals = 10000000;
    std::vector<double> yvals = {0, 0, 0, 0, 0};
    start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < evals; i++) vector
    std::cout << "[SYMENGINE]  Time for call llvm gradient: " << duration << " microseconds\n";
    */
    // new test
    
    std::vector<RCP<const Basic>> x =  {symbol("x"), symbol("y")};
    std::vector<RCP<const Basic>> x0 =  {symbol("x_0")};
    std::vector<RCP<const Basic>> xf =  {symbol("x_f"), symbol("y_f")};
    std::vector<RCP<const Basic>> u =  {symbol("u")};
    std::vector<RCP<const Basic>> p =  {symbol("p")};
    RCP<const Basic> t = symbol("t");
    RCP<const Basic> M = parse("x_f * y_f");
    RCP<const Basic> f1 = parse("p * (x + y * cos(t)) + u");
    RCP<const Basic> f2 = parse("x * exp(u) / (u**2 + 1)");
    RCP<const Basic> l1 = parse("x_f**2 * x_0**2");
    
    double xuvals[3] = {1, 1, 1};
    double xvals[2] = {1, 1};
    double uvals[1] = {1};
    double pvals[2] = {1};
    double tval[1] = {1};

    PhaseVariables variables(x, u, x0, xf, t);
    GlobalParameters parameters(p);
    ModelParameters mparameters;
    PhaseFunctions functions(M, {}, {f1, f2}, {}, {}, {});

    FullSweep fs(variables, parameters, mparameters, functions);
    BoundarySweep bs(variables, parameters, mparameters, functions);
    Linkage lk(0, 1, variables, variables, parameters, mparameters, {l1}, {}, {});
    
    fs.fillInputData(xuvals, pvals, tval);
    fs.callHess();
    fs.callEval();

    auto out = evalDiff2(xvals, uvals, pvals, tval);

    std::string file = "function.elf";
    saveBinaryLLVM(fs.llvmEval, file);
    fs.llvmEval = readBinaryLLVM(file);
    
    fs.callEval();
    volatile int i = 10;
    fs.callEval();

    auto start = std::chrono::high_resolution_clock::now();
    auto fullsweep = attitudeControlInit();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time taken for full NEW Attitude Control Derivative Calculations + Compilation: " << duration.count() << " ms" << std::endl;
    std::cout << "Time taken for full GDOPT Attitude Control Derivative Calculations: " << 106 << " ms" << std::endl;
    std::cout << "Time taken for full GDOPT Attitude Control Compilation: " << 2652 << " ms" << std::endl;

    const int iterations = 1000000;
    start = std::chrono::high_resolution_clock::now();
    attitudeControlHessCall(fullsweep, iterations);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time taken for NEW Attitude Control Hessian Calls: " << duration.count() << " ms" << std::endl;

    start = std::chrono::high_resolution_clock::now();
    attitudeControlGDOPTHessCall(iterations);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time taken for GDOPT Attitude Control Hessian Calls: " << duration.count() << " ms" << std::endl;

    /*
    std::cout << "Gradient\n";
    fs.callGrad();
    for (int i = 0; i < 6; i++) {
        std::cout << fs.gradData[i] << std::endl;
    }
    std::cout << "Hessian\n";

        
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1000; i++)
        fs.callHess();
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "[NEW GDOPT]  Time for call hessian GDOPT : " << duration << " microseconds\n"; */
    return 0;
}
