#ifndef OPT_OM_GDOP_PROBLEM
#define OPT_OM_GDOP_PROBLEM

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "evaluations.h"
#include "debug_om.h"

class FullSweep_OM : public FullSweep {
public:
    DATA* data;
    threadData_t* threadData;
    InfoGDOP& info;

    FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds,
                 DATA* data, threadData_t* threadData, InfoGDOP& info);
                           
    void callback_eval(const F64* xu_nlp, const F64* p) override;
    void callback_jac(const F64* xu_nlp, const F64* p) override;
    void callback_hes(const F64* xu_nlp, const F64* p) override;
};

class BoundarySweep_OM : public BoundarySweep {
public:
    DATA* data;
    threadData_t* threadData;
    InfoGDOP& info;
    
    BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                     DATA* data, threadData_t* threadData, InfoGDOP& info);
    void callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
    void callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
    void callback_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
};

Problem create_gdop(DATA* data, threadData_t* threadData, InfoGDOP& info, Mesh& mesh);

#endif // OPT_OM_GDOP_PROBLEM
