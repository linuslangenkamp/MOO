#ifndef OPT_GDOP_ORCHESTRATOR_H
#define OPT_GDOP_ORCHESTRATOR_H

#include <src/nlp/nlp_solver.h>
#include "gdop.h"

namespace GDOP {

class Orchestrator {
public:
    GDOP& gdop;
    std::unique_ptr<Strategies> strategies;
    NLP::NLPSolver& solver;

    Orchestrator(GDOP& gdop,
                 std::unique_ptr<Strategies> strategies,
                 NLP::NLPSolver& solver)
    : gdop(gdop),
      strategies(std::move(strategies)),
      solver(solver) {}

    virtual ~Orchestrator() = default;

    virtual void optimize() = 0;
};

class MeshRefinementOrchestrator : public Orchestrator {
public:
    MeshRefinementOrchestrator(GDOP& gdop,
                               std::unique_ptr<Strategies> strategies,
                               NLP::NLPSolver& solver)
    : Orchestrator(gdop, std::move(strategies), solver) {}

    void optimize() override;
};

} // namespace GDOP

#endif // OPT_GDOP_ORCHESTRATOR_H
