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

    FullSweep_OM(GDOP::BlockLFG&& lfg_in,
              std::unique_ptr<AugmentedHessianLFG> aug_hes,
              std::unique_ptr<AugmentedParameterHessian> aug_pp_hes,
              const GDOP::ProblemConstants& pc,
              InfoGDOP& info);

    void callback_eval(const f64* xu_nlp, const f64* p) override;
    void callback_jac(const f64* xu_nlp, const f64* p) override;
    void callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, f64* lambda) override;
};

class BoundarySweep_OM : public GDOP::BoundarySweep {
public:
    InfoGDOP& info;

    BoundarySweep_OM(GDOP::BlockMR&& mr_in,
                  std::unique_ptr<AugmentedHessianMR> aug_hes,
                  const GDOP::ProblemConstants& pc,
                  InfoGDOP& info);

    void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) override;
    void callback_aug_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, f64* lambda) override;
};

GDOP::Problem create_gdop(InfoGDOP& info, Mesh& mesh, Collocation& collocation);

} // namespace OpenModelica

#endif // OPT_OM_GDOP_PROBLEM_H
