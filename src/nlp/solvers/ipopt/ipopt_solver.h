#ifndef OPT_IPOPT_SOLVER_H
#define OPT_IPOPT_SOLVER_H

#include <memory>

#include <nlp/nlp_solver.h>
#include "base/nlp_structs.h"
#include "ipopt_adapter.h"

class IpoptSolver : NLPSolver {
public:
    IpoptSolver(NLP& nlp, NLPSolverFlags& solver_flags);

    virtual ~IpoptSolver() = default;

    Ipopt::SmartPtr<IpoptAdapter> adapter;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app;

    void optimize() override;
    void init_IpoptApplication();
};

#endif // OPT_IPOPT_SOLVER_H
