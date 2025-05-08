#ifndef OPT_IPOPT_SOLVER_H
#define OPT_IPOPT_SOLVER_H

#include <memory>

#include <nlp/nlp_solver.h>
#include "ipopt_adapter.h"


class IpoptSolver : NLPSolver {
public:
    IpoptSolver(std::shared_ptr<NLP> nlp, std::shared_ptr<std::unordered_map<std::string, std::string>> solver_settings);

    virtual ~IpoptSolver() = default;

    Ipopt::SmartPtr<IpoptAdapter> adapter;
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app;

    void optimize() override;
    void initIpoptApplication();
};

#endif // OPT_IPOPT_SOLVER_H
