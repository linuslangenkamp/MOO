#include "nlp_solver_flags.h"

namespace NLP {

NLPSolverFlags::NLPSolverFlags(int argc, char** argv) {
    /* default settings : add when needed */
    solver_settings["Hessian"] = "Exact";
    solver_settings["Tolerance"] = "1e-10";
    solver_settings["Iterations"] = "5000";
    solver_settings["CPUTime"] = "3600";
    solver_settings["LinearSolver"] = "MUMPS";
    solver_settings["NLPSolver"] = "Ipopt";
    solver_settings["IpoptDerivativeTest"] = "false";
    solver_settings["WarmStart"] = "false";
}

void NLPSolverFlags::print() const {
    FixedTableFormat<2> ftf = {
        {25,            15},
        {Align::Center, Align::Center}
    };

    LOG_START_MODULE(ftf, "NLP Solver Flags");

    LOG_HEADER(ftf, "Flag", "Value");
    LOG_DASHES(ftf);

    for (const auto& [flag, value] : solver_settings) {
        LOG_ROW(ftf, flag, value);
    }

    LOG_DASHES_LN(ftf);
}

std::string NLPSolverFlags::get(const std::string& flag) const {
    auto it = solver_settings.find(flag);
    return it != solver_settings.end() ? it->second : "";
}

void NLPSolverFlags::set(const std::string& flag, const std::string& value) {
    solver_settings[flag] = value;
}

bool NLPSolverFlags::check_flag(const std::string& flag, const std::string& value) const {
    auto it = solver_settings.find(flag);
    return it != solver_settings.end() && it->second == value;
}

f64 NLPSolverFlags::get_flag_f64_fallback(const std::string& flag, const f64 fallback) const {
    auto it = solver_settings.find(flag);
    if (it != solver_settings.end()) {
        try {
            return std::stod(it->second);
        } catch (...) {
            LOG_WARNING("Invalid f64 flag: {} with value {} - defaulting to {}.", flag, it->second, fallback);
            return fallback;
        }
    }
    return fallback;
}

std::string NLPSolverFlags::get_flag_string_fallback(const std::string& flag, const std::string& fallback) const {
    auto it = solver_settings.find(flag);
    if (it != solver_settings.end()) {
        return it->second;
    }
    return fallback;
}

} // namespace NLP
