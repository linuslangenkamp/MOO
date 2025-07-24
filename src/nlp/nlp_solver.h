#ifndef OPT_NLP_SOLVER_H
#define OPT_NLP_SOLVER_H

#include <unordered_map>
#include <memory>

#include "nlp/solvers/nlp_solver_settings.h"
#include "nlp.h"

namespace NLP {

class NLPSolver {
public:
    NLPSolver(NLP& nlp, NLPSolverSettings& solver_settings)
        : nlp(nlp), solver_settings(solver_settings) {}

    virtual ~NLPSolver() = default;

    NLP& nlp;
    NLPSolverSettings& solver_settings;

    virtual void optimize() = 0;
};

}

#endif // OPT_NLP_SOLVER_H
