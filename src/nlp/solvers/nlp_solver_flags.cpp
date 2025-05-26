#include "nlp_solver_flags.h"

NLPSolverFlags::NLPSolverFlags(int argc, char** argv) {
    /* Default settings : add when needed */
    solver_settings["Hessian"] = "Exact";
    solver_settings["Tolerance"] = "1e-10";
    solver_settings["CPUTime"] = "3600";
    solver_settings["LinearSolver"] = "MUMPS";
    solver_settings["NLPSolver"] = "Ipopt";

    /* parse from CLI? */
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--hessian" && i + 1 < argc) {
            solver_settings["Hessian"] = argv[++i];
        } else if (arg == "--tolerance" && i + 1 < argc) {
            solver_settings["Tolerance"] = argv[++i];
        }
    }
}

void NLPSolverFlags::print() const {
    size_t max_len = 5;
    for (const auto& [flag, _] : solver_settings) {
        max_len = std::max(max_len, flag.length());
    }

    int col1_width = static_cast<int>(max_len) + 2;
    std::string header = 
        std::string("Flag") + std::string(col1_width - 4, ' ') + "| Value";

    std::cout << header << std::endl;
    std::cout << std::string(header.length(), '-') << std::endl;

    for (const auto& [flag, value] : solver_settings) {
        std::cout << std::left << std::setw(static_cast<int>(max_len) + 2)
                  << flag << ": " << value << std::endl;
    }
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
            std::cout << "Invalid F64 flag: \"" << flag 
                      << "\" with value \"" << it->second 
                      << "\" — defaulting to " << fallback << std::endl;
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
