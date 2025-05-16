#ifndef OPT_OM_GDOP_PROBLEM_H
#define OPT_OM_GDOP_PROBLEM_H

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "evaluations.h"
#include "debug_om.h"

class FullSweep_OM : public FullSweep {
public:
    DATA* data;
    threadData_t* threadData;
    InfoGDOP& info;

    FullSweep_OM(FixedVector<FunctionLFG>&& lfg, std::unique_ptr<AugmentedHessianLFG> aug_hes, std::unique_ptr<AugmentedParameterHessian> aug_pp_hes,
     Collocation& collocation, Mesh& mesh, FixedVector<Bounds>&& g_bounds, DATA* data, threadData_t* threadData, InfoGDOP& info);
                           
    void callback_eval(const F64* xu_nlp, const F64* p) override;
    void callback_jac(const F64* xu_nlp, const F64* p) override;
    void callback_aug_hes(const F64* xu_nlp, const F64* p, const FixedField<F64, 2>& lagrange_factors, const F64* lambda) override;
};

class BoundarySweep_OM : public BoundarySweep {
public:
    DATA* data;
    threadData_t* threadData;
    InfoGDOP& info;

    BoundarySweep_OM(FixedVector<FunctionMR>&& mr, std::unique_ptr<AugmentedHessianMR> aug_hes, Mesh& mesh,
                     FixedVector<Bounds>&& r_bounds, DATA* data, threadData_t* threadData, InfoGDOP& info);
    void callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
    void callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
    void callback_aug_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p, const F64 mayer_factor, const F64* lambda) override;
};

Problem create_gdop(DATA* data, threadData_t* threadData, InfoGDOP& info, Mesh& mesh, Collocation& fLGR);

#endif // OPT_OM_GDOP_PROBLEM_H
