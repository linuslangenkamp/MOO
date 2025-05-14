#ifndef OPT_NLP_SOLVER_FLAGS_H
#define OPT_NLP_SOLVER_FLAGS_H

#include <unordered_map>
#include <iostream>
#include <string>
#include <iomanip>

#include "base/util.h"

// TODO: this is WIP
/* struct to handle all NLP Solver flags; special solver flags are named via solvername_flagname, e.g. Ipopt_MuInit */
class NLPSolverFlags {
public:
    // Constructor: initializes with defaults and parses CLI args
    NLPSolverFlags(int argc, char** argv);

    // Print all flags
    void print() const;

    // Get value for a flag (returns empty string if not set)
    std::string get(const std::string& flag) const;

    // Set flag to a specific value
    void set(const std::string& flag, const std::string& value);

    // Check if a flag matches a specific value
    bool check_flag(const std::string& flag, const std::string& value) const;

    // Get flag as T, or fallback if not found or invalid
    F64 get_flag_f64_fallback(const std::string& flag, const F64 fallback) const;
    std::string get_flag_string_fallback(const std::string& flag, const std::string& fallback) const;

private:
    std::unordered_map<std::string, std::string> solver_settings;
};

#endif // OPT_NLP_SOLVER_FLAGS_H
