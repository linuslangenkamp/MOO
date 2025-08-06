#ifndef MOO_C_PROBLEM_H
#define MOO_C_PROBLEM_H

#include <nlp/instances/gdop/problem.h>
#include <nlp/instances/gdop/gdop.h>

#include <src/interfaces/c/structures.h>

namespace C {

class FullSweep : public GDOP::FullSweep {
public:
    FullSweep(GDOP::FullSweepLayout&& layout_in,
              const GDOP::ProblemConstants& pc,
              c_callbacks_t* c_callbacks)
        : GDOP::FullSweep(std::move(layout_in), pc), c_callbacks(c_callbacks) {};

    void callback_eval(const f64* xu_nlp, const f64* p) override;
    void callback_jac(const f64* xu_nlp, const f64* p) override;
    void callback_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, const f64* lambda) override;

private:
    c_callbacks_t* c_callbacks;
};

class BoundarySweep : public GDOP::BoundarySweep {
public:
    BoundarySweep(GDOP::BoundarySweepLayout&& layout_in,
                  const GDOP::ProblemConstants& pc,
                  c_callbacks_t* c_callbacks)
        : GDOP::BoundarySweep(std::move(layout_in), pc), c_callbacks(c_callbacks) {};

    void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, const f64* lambda) override;

private:
    c_callbacks_t* c_callbacks;
};

GDOP::Problem create_gdop(c_problem_t* c_problem, const Mesh& mesh);

} // namespace C

#endif // MOO_C_PROBLEM_H
