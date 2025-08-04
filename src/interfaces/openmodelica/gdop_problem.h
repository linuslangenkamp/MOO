#ifndef OPT_OM_GDOP_PROBLEM_H
#define OPT_OM_GDOP_PROBLEM_H

#include "simulation_data.h"
#include "simulation/simulation_runtime.h"
#include "simulation/solver/gbode_main.h"
#include "simulation/solver/external_input.h"

#include <nlp/instances/gdop/problem.h>
#include <nlp/instances/gdop/gdop.h>

#include "evaluations.h"
#include "strategies.h"
#include "debug_om.h"

namespace OpenModelica {

class FullSweep_OM : public GDOP::FullSweep {
public:
    InfoGDOP& info;

    FullSweep_OM(GDOP::FullSweepLayout&& lfg_in,
                 const GDOP::ProblemConstants& pc,
                 InfoGDOP& info);

    void callback_eval(const f64* xu_nlp, const f64* p) override;
    void callback_jac(const f64* xu_nlp, const f64* p) override;
    void callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, const f64* lambda) override;
};

class BoundarySweep_OM : public GDOP::BoundarySweep {
public:
    InfoGDOP& info;

    BoundarySweep_OM(GDOP::BoundarySweepLayout&& mr_in,
                     const GDOP::ProblemConstants& pc,
                     InfoGDOP& info);

    void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_aug_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, const f64* lambda) override;
};

GDOP::Problem create_gdop(InfoGDOP& info, Mesh& mesh);

} // namespace OpenModelica

#endif // OPT_OM_GDOP_PROBLEM_H
