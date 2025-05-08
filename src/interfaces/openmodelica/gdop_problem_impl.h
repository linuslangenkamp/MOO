#ifndef OPT_OM_GDOP_PROBLEM_IMPL
#define OPT_OM_GDOP_PROBLEM_IMPL

#include <nlp/instances/gdop/problem.h>

class FullSweepOM : public FullSweep {
public:
    FullSweepOM(FixedVector<FunctionLFG>&& lfg, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds);
    void callbackEval(const double* xu_nlp, const double* p) override;
    void callbackJac(const double* xu_nlp, const double* p) override;
    void callbackHes(const double* xu_nlp, const double* p) override;
};

class BoundarySweepOM : public BoundarySweep {
public:
    BoundarySweepOM(FixedVector<FunctionMR>&& mr, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds);
    void callbackEval(const double* x0_nlp, const double* xf_nlp, const double* p) override;
    void callbackJac(const double* x0_nlp, const double* xf_nlp, const double* p) override;
    void callbackHes(const double* x0_nlp, const double* xf_nlp, const double* p) override;
};

#endif // OPT_OM_GDOP_PROBLEM_IMPL