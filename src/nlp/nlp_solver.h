#ifndef OPT_NLP_SOLVER_H
#define OPT_NLP_SOLVER_H

#include <unordered_map>
#include <memory>

#include "nlp/solvers/nlp_solver_flags.h"
#include "nlp.h"

namespace NLP {

class NLPSolver {
public:
    NLPSolver(NLP& nlp, NLPSolverFlags& solver_flags)
        : nlp(nlp), solver_flags(solver_flags) {}

    virtual ~NLPSolver() = default;

    NLP& nlp;
    NLPSolverFlags& solver_flags;

    virtual void optimize() = 0;
    // virtual void set_default_flags() = 0;
    // virtual void set_flags(NLPSolverFlags& new_flags) = 0;
};

}

#endif // OPT_NLP_SOLVER_H
