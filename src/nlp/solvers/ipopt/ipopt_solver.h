#ifndef OPT_IPOPT_SOLVER_H
#define OPT_IPOPT_SOLVER_H

#include <memory>

#include <nlp/nlp_solver.h>
#include "base/nlp_structs.h"
#include "base/log.h"
#include "ipopt_adapter.h"

namespace IpoptSolver {

class IpoptSolver : public NLP::NLPSolver {
public:
    IpoptSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings);

    virtual ~IpoptSolver() = default;

    Ipopt::SmartPtr<IpoptAdapter> adapter;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app;

    void optimize() override;
    void init_application();
    void set_settings();
};

} // namespace IpoptSolver

#endif // OPT_IPOPT_SOLVER_H
