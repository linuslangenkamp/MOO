#ifndef OPT_NLP_SOLVER_H
#define OPT_NLP_SOLVER_H

#include <unordered_map>
#include <memory>

#include "nlp.h"


class NLPSolver {
public:
    NLPSolver(NLP& nlp, std::unordered_map<std::string, std::string>& solver_settings)
        : nlp(nlp), solver_settings(solver_settings) {}

    virtual ~NLPSolver() = default;

    NLP& nlp;
    std::unordered_map<std::string, std::string>& solver_settings;

    virtual void optimize() = 0;

    /* Check if solver_settings contains a specific key-value pair */
    bool check_flag(const std::string& key, const std::string& value) const {
        auto it = solver_settings.find(key);
        return it != solver_settings.end() && it->second == value;
    }
};

#endif // OPT_NLP_SOLVER_H
