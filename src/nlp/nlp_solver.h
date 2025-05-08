#ifndef OPT_NLP_SOLVER_H
#define OPT_NLP_SOLVER_H

#include <unordered_map>
#include <memory>

#include "nlp.h"


class NLPSolver {
public:
    NLPSolver(std::shared_ptr<NLP> nlp, std::shared_ptr<std::unordered_map<std::string, std::string>> solver_settings)
        : nlp(nlp), solver_settings(solver_settings) {}

    virtual ~NLPSolver() = default;

    std::shared_ptr<NLP> nlp;
    std::shared_ptr<std::unordered_map<std::string, std::string>> solver_settings;

    virtual void optimize() = 0;
};

#endif // OPT_NLP_SOLVER_H
