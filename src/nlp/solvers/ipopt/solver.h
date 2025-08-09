#ifndef OPT_IPOPT_SOLVER_H
#define OPT_IPOPT_SOLVER_H

#include <memory>

#include <nlp/nlp_solver.h>
#include "base/nlp_structs.h"
#include "base/log.h"


namespace IpoptSolver {

struct IpoptSolverData;

class IpoptSolver : public NLP::NLPSolver {
public:
    IpoptSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings);

    virtual ~IpoptSolver();

    void optimize() override;
    void init_application();
    void set_settings();

private:
   IpoptSolverData* data;
};

} // namespace IpoptSolver

#endif // OPT_IPOPT_SOLVER_H
