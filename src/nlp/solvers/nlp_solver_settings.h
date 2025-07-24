#ifndef NLP_SOLVER_SETTINGS_H
#define NLP_SOLVER_SETTINGS_H

#include <string>
#include <unordered_map>
#include <variant>
#include <optional>
#include <iostream>

#include <src/base/util.h>
#include <src/base/log.h>

namespace NLP {

enum class HessianOption {
    Exact,
    LBFGS,
    CONST
};

enum class LinearSolverOption {
    MUMPS,
    MA27,
    MA57,
    MA77,
    MA86,
    MA97
};

enum class NLPSolverOption {
    Ipopt
};

enum class Option {
/* HessianOption         */    Hessian,
/* f64                   */    Tolerance,
/* int                   */    Iterations,
/* f64                   */    CPUTime,
/* LinearSolverOption    */    LinearSolver,
/* NLPSolverOption       */    NLPSolver,
/* bool                  */    IpoptDerivativeTest,
/* bool                  */    WarmStart,
/* bool                  */    QP,
};

using OptionValue = std::variant<std::string, f64, int, bool, HessianOption, LinearSolverOption, NLPSolverOption>;

class NLPSolverSettings {
public:
    NLPSolverSettings(int argc, char** argv);

    void print() const;

    void set(Option option, const OptionValue& value);
    const OptionValue& get(Option option) const;

    bool option_is_true(Option option) const;
    bool option_matches(Option option, const std::string& str) const;

    template<typename T>
    T get_or_default(Option option) const;

private:
    std::unordered_map<Option, OptionValue> settings;
};

// Option enum to string
std::string to_string(Option option);

// string to Option enum
std::optional<Option> option_from_string(const std::string& name);

extern const std::unordered_map<Option, OptionValue> default_settings;

} // namespace NLP

#endif // NLP_SOLVER_SETTINGS_H
