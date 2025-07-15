#ifndef OPT_NLP_SOLVER_FLAGS_H
#define OPT_NLP_SOLVER_FLAGS_H

#include <unordered_map>
#include <iostream>
#include <string>
#include <iomanip>

#include "base/util.h"
#include "base/log.h"

namespace NLP {

// TODO: this is WIP
/* struct to handle all NLP Solver flags; special solver flags are named via solvername_flagname, e.g. Ipopt_MuInit */
class NLPSolverFlags {
public:
    // constructor: initializes with defaults and parses CLI args
    NLPSolverFlags(int argc, char** argv);

    // print all flags
    void print() const;

    // get value for a flag (returns empty string if not set)
    std::string get(const std::string& flag) const;

    // set flag to a specific value
    void set(const std::string& flag, const std::string& value);

    // check if a flag matches a specific value
    bool check_flag(const std::string& flag, const std::string& value) const;

    // get flag as T, or fallback if not found or invalid
    f64 get_flag_f64_fallback(const std::string& flag, const f64 fallback) const;
    std::string get_flag_string_fallback(const std::string& flag, const std::string& fallback) const;
    std::string get_flag_bool_as_string(const std::string& flag_name, bool default_value) const;
private:
    std::unordered_map<std::string, std::string> solver_settings;
};

} // namespace NLP

#endif // OPT_NLP_SOLVER_FLAGS_H
