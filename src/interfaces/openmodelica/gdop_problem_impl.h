#ifndef OPT_OM_GDOP_PROBLEM_IMPL
#define OPT_OM_GDOP_PROBLEM_IMPL

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

class FullSweepOM : public FullSweep {
public:
    FullSweepOM(FixedVector<FunctionLFG>&& lfg, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds);
    void callback_eval(const double* xu_nlp, const double* p) override;
    void callback_jac(const double* xu_nlp, const double* p) override;
    void callback_hes(const double* xu_nlp, const double* p) override;
};

class BoundarySweepOM : public BoundarySweep {
public:
    BoundarySweepOM(FixedVector<FunctionMR>&& mr, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds);
    void callback_eval(const double* x0_nlp, const double* xf_nlp, const double* p) override;
    void callback_jac(const double* x0_nlp, const double* xf_nlp, const double* p) override;
    void callback_hes(const double* x0_nlp, const double* xf_nlp, const double* p) override;
};

std::shared_ptr<Problem> create_gdop_om(DATA* data, std::shared_ptr<Mesh> mesh);

#endif // OPT_OM_GDOP_PROBLEM_IMPL