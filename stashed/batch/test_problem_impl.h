#ifndef OPT_TEST_H
#define OPT_TEST_H

#include "problem.h"

// a first cpp implementation of the GDOP Problem

class FullSweepTestImpl : public FullSweep {
public:
    FullSweepTestImpl(FixedVector<FunctionLFG>&& lfg, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds);

    void callback_eval(const double* xu_nlp, const double* p) override;

    void callbackJac(const double* xu_nlp, const double* p) override;

    void callback_hes(const double* xu_nlp, const double* p) override;
};

class BoundarySweepTestImpl : public BoundarySweep {
public:
    BoundarySweepTestImpl(FixedVector<FunctionMR>&& mr, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds);

    void callback_eval(const double* x0_nlp, const double* xf_nlp, const double* p) override;

    void callbackJac(const double* x0_nlp, const double* xf_nlp, const double* p) override;

    void callback_hes(const double* x0_nlp, const double* xf_nlp, const double* p) override;
};

#endif // OPT_TEST_H