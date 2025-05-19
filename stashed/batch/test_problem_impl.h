#ifndef OPT_TEST_H
#define OPT_TEST_H

#include "problem.h"

// a first cpp implementation of the GDOP Problem

class FullSweepTestImpl : public FullSweep {
public:
    FullSweepTestImpl(FixedVector<FunctionLFG>&& lfg, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& g_bounds);

    void callback_eval(const f64* xu_nlp, const f64* p) override;

    void callbackJac(const f64* xu_nlp, const f64* p) override;

    void callback_hes(const f64* xu_nlp, const f64* p) override;
};

class BoundarySweepTestImpl : public BoundarySweep {
public:
    BoundarySweepTestImpl(FixedVector<FunctionMR>&& mr, std::shared_ptr<Mesh> mesh, FixedVector<Bounds>& r_bounds);

    void callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) override;

    void callbackJac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) override;

    void callback_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p) override;
};

#endif // OPT_TEST_H