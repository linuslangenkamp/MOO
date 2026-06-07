// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2026 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include "solver.h"

#include <algorithm>
#include <chrono>
#include <exception>
#include <interfaces/C/Uno_C_API.h>

#include <base/log.h>
#include <base/timing.h>

#include "adapter.h"

namespace UnoSolver {

namespace {

struct LoggerBuffer {
    std::string pending;
    bool logged_identical_bounds = false;
};

bool is_identical_bounds_message(const std::string &line) {
    constexpr const char *prefix = "Variable x";
    return line.compare(0, std::char_traits<char>::length(prefix), prefix) == 0 && line.find(" has identical bounds") != std::string::npos;
}

bool should_log_line(LoggerBuffer &logger_buffer, std::string &line) {
    if (!is_identical_bounds_message(line)) {
        return true;
    }

    if (logger_buffer.logged_identical_bounds) {
        return false;
    }

    Log::warning("At least one variable has the same upper and lower bounds.");
    logger_buffer.logged_identical_bounds = true;

    return false;
}

void flush_logger_buffer(LoggerBuffer &logger_buffer) {
    if (!logger_buffer.pending.empty()) {
        if (should_log_line(logger_buffer, logger_buffer.pending)) {
            Log::info(logger_buffer.pending);
        }
        logger_buffer.pending.clear();
    }
}

uno_int logger_stream_callback(const char *buffer, uno_int length, void *user_data) {
    if (!buffer || length <= 0) {
        return 0;
    }
    auto &logger_buffer = *static_cast<LoggerBuffer *>(user_data);
    logger_buffer.pending.append(buffer, static_cast<size_t>(length));

    size_t pos = 0;
    while (pos < logger_buffer.pending.size()) {
        size_t newline = logger_buffer.pending.find('\n', pos);
        if (newline == std::string::npos) {
            break;
        }

        std::string line = logger_buffer.pending.substr(pos, newline - pos);
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        if (!line.empty() && should_log_line(logger_buffer, line)) {
            Log::info(line);
        }
        pos = newline + 1;
    }
    logger_buffer.pending.erase(0, pos);

    return length;
}

const char *uno_optimization_status_to_string(uno_int status) {
    if (status == UNO_SUCCESS) {
        return "UNO_SUCCESS";
    } else if (status == UNO_ITERATION_LIMIT) {
        return "UNO_ITERATION_LIMIT";
    } else if (status == UNO_TIME_LIMIT) {
        return "UNO_TIME_LIMIT";
    } else if (status == UNO_EVALUATION_ERROR) {
        return "UNO_EVALUATION_ERROR";
    } else if (status == UNO_ALGORITHMIC_ERROR) {
        return "UNO_ALGORITHMIC_ERROR";
    } else if (status == UNO_USER_TERMINATION) {
        return "UNO_USER_TERMINATION";
    } else {
        return "UNO_UNKNOWN_OPTIMIZATION_STATUS";
    }
}

const char *uno_solution_status_to_string(uno_int status) {
    if (status == UNO_NOT_OPTIMAL) {
        return "UNO_NOT_OPTIMAL";
    } else if (status == UNO_FEASIBLE_KKT_POINT) {
        return "UNO_FEASIBLE_KKT_POINT";
    } else if (status == UNO_FEASIBLE_FJ_POINT) {
        return "UNO_FEASIBLE_FJ_POINT";
    } else if (status == UNO_INFEASIBLE_STATIONARY_POINT) {
        return "UNO_INFEASIBLE_STATIONARY_POINT";
    } else if (status == UNO_FEASIBLE_SMALL_STEP) {
        return "UNO_FEASIBLE_SMALL_STEP";
    } else if (status == UNO_INFEASIBLE_SMALL_STEP) {
        return "UNO_INFEASIBLE_SMALL_STEP";
    } else if (status == UNO_UNBOUNDED) {
        return "UNO_UNBOUNDED";
    } else {
        return "UNO_UNKNOWN_SOLUTION_STATUS";
    }
}

} // namespace

struct UnoSolverData {
    UnoAdapter adapter;
    void *solver = nullptr;
    void *model = nullptr;
    LoggerBuffer logger_buffer;
    f64 total_wall_nano = 0.0;

    explicit UnoSolverData(NLP::NLP &nlp)
        : adapter(nlp),
          solver(uno_create_solver()) {}

    ~UnoSolverData() {
        if (model) {
            uno_destroy_model(model);
        }
        if (solver) {
            uno_destroy_solver(solver);
        }
    }
};

UnoSolver::UnoSolver(NLP::NLP &nlp, NLP::NLPSolverSettings &solver_settings)
    : NLPSolver(nlp, solver_settings),
      udata(new UnoSolverData(nlp)) {
    if (!udata->solver) {
        Log::error("[Uno Interface] Failed to create Uno solver.");
        abort();
    }
}

UnoSolver::~UnoSolver() {
    delete udata;
}

void UnoSolver::optimize() {
    ScopedTimer<UnoTimingNode, UnoSolver *> timer("UnoSolver::optimize", this);
    const auto start = std::chrono::high_resolution_clock::now();

    if (udata->model) {
        uno_destroy_model(udata->model);
        udata->model = nullptr;
    }

    try {
        set_settings();
        udata->model = udata->adapter.create_model();
        if (!udata->model) {
            abort();
        }

        udata->adapter.reset_timing();
        uno_optimize(udata->solver, udata->model);
        flush_logger_buffer(udata->logger_buffer);

        log_status();
        udata->adapter.finalize_solution(udata->solver, get_return_code());
    } catch (const std::exception &e) {
        flush_logger_buffer(udata->logger_buffer);
        Log::error("[Uno Interface] Unhandled exception: {}", e.what());
        std::abort();
    } catch (...) {
        flush_logger_buffer(udata->logger_buffer);
        Log::error("[Uno Interface] Unhandled non-standard exception.");
        std::abort();
    }

    udata->total_wall_nano = std::chrono::duration<f64, std::nano>(std::chrono::high_resolution_clock::now() - start).count();
}

void UnoSolver::set_settings() {
    uno_set_logger_stream_callback(&logger_stream_callback, &udata->logger_buffer);

    const std::string preset = solver_settings.get_or_default<std::string>(NLP::Option::UnoPreset);
    if (!uno_set_solver_preset(udata->solver, preset.c_str())) {
        Log::warning("[Uno Interface] Failed to set Uno preset '{}'.", preset);
    }
    uno_set_solver_string_option(udata->solver, "logger", "INFO");

    uno_set_solver_integer_option(udata->solver, "max_iterations", solver_settings.get_or_default<int>(NLP::Option::Iterations));
    uno_set_solver_double_option(udata->solver, "time_limit", solver_settings.get_or_default<f64>(NLP::Option::CPUTime));

    const f64 tol = solver_settings.get_or_default<f64>(NLP::Option::Tolerance);
    uno_set_solver_double_option(udata->solver, "primal_tolerance", tol);
    uno_set_solver_double_option(udata->solver, "dual_tolerance", tol);
    uno_set_solver_double_option(udata->solver, "loose_primal_tolerance", tol * 1e3);
    uno_set_solver_double_option(udata->solver, "loose_dual_tolerance", tol * 1e3);

    switch (solver_settings.get_or_default<NLP::HessianOption>(NLP::Option::Hessian)) {
        case NLP::HessianOption::Exact:
        case NLP::HessianOption::Const:
            uno_set_solver_string_option(udata->solver, "hessian_model", "exact");
            break;
        case NLP::HessianOption::LBFGS:
            uno_set_solver_string_option(udata->solver, "hessian_model", "LBFGS");
            break;
    }

    uno_set_solver_bool_option(udata->solver, "print_solution", false);
}

NLP::ReturnCode UnoSolver::get_return_code() const {
    void *solver = udata->solver;

    const uno_int optimization_status = uno_get_optimization_status(solver);
    const uno_int solution_status = uno_get_solution_status(solver);

    if (optimization_status == UNO_ITERATION_LIMIT) {
        return NLP::ReturnCode::ITERATION_LIMIT_EXCEEDED;
    } else if (optimization_status == UNO_TIME_LIMIT) {
        return NLP::ReturnCode::TIME_LIMIT_EXCEEDED;
    } else if (optimization_status == UNO_EVALUATION_ERROR || optimization_status == UNO_ALGORITHMIC_ERROR || optimization_status == UNO_USER_TERMINATION) {
        return NLP::ReturnCode::GENERIC_FAILURE;
    } else if (optimization_status != UNO_SUCCESS) {
        return NLP::ReturnCode::GENERIC_FAILURE;
    }

    if (solution_status == UNO_FEASIBLE_KKT_POINT) {
        return NLP::ReturnCode::OPTIMAL;
    } else if (solution_status == UNO_FEASIBLE_FJ_POINT) {
        return NLP::ReturnCode::ACCEPTABLE;
    } else if (solution_status == UNO_FEASIBLE_SMALL_STEP) {
        return NLP::ReturnCode::STEP_TOO_SMALL;
    } else if (solution_status == UNO_INFEASIBLE_STATIONARY_POINT || solution_status == UNO_INFEASIBLE_SMALL_STEP) {
        return NLP::ReturnCode::INFEASIBLE;
    } else if (solution_status == UNO_UNBOUNDED) {
        return NLP::ReturnCode::DIVERGENCE;
    } else if (solution_status == UNO_NOT_OPTIMAL) {
        return NLP::ReturnCode::UNKNOWN_SUCCESS;
    } else {
        return NLP::ReturnCode::UNKNOWN_SUCCESS;
    }
}

void UnoSolver::log_status() const {
    void *solver = udata->solver;

    const uno_int optimization_status = uno_get_optimization_status(solver);
    const uno_int solution_status = uno_get_solution_status(solver);

    const char *optimization_status_name = uno_optimization_status_to_string(optimization_status);
    const char *solution_status_name = uno_solution_status_to_string(solution_status);

    if (optimization_status == 0) {
        Log::success("[Uno Interface] optimization status: {}, solution status: {}.", optimization_status_name, solution_status_name);
    } else {
        Log::error("[Uno Interface] optimization status: {}, solution status: {}.", optimization_status_name, solution_status_name);
    }
}

int UnoSolver::get_iterations() const {
    return static_cast<int>(uno_get_number_iterations(udata->solver));
}

f64 UnoSolver::get_total_time() const {
    return udata->total_wall_nano;
}

f64 UnoSolver::get_solver_time() const {
    return std::max(0.0, get_total_time() - get_callback_time());
}

f64 UnoSolver::get_callback_time() const {
    return udata->adapter.get_callback_timing().total_nano();
}

UnoTimingNode::UnoTimingNode(std::string n, TimingNode *p, UnoSolver *uno_solver)
    : TimingNode(n, p),
      uno_solver(uno_solver) {}

void UnoTimingNode::finalize() {
    const auto &timing = uno_solver->udata->adapter.get_callback_timing();
    const f64 total_nano = uno_solver->udata->total_wall_nano;
    const f64 callback_nano = timing.total_nano();
    const f64 internal_nano = std::max(0.0, total_nano - callback_nano);

    VirtualTimer timer_wall_total("Uno Total", total_nano);
    {
        { VirtualTimer timer_wall_uno_self("Uno Internal", internal_nano); }
        {
            VirtualTimer timer_wall_func("Uno Callbacks", callback_nano);
            { VirtualTimer<CountedTimingNode, int> timer_wall_f("Objective", timing.objective_nano, static_cast<int>(timing.objective_count)); }
            { VirtualTimer<CountedTimingNode, int> timer_wall_grad_f("Objective Gradient", timing.objective_gradient_nano, static_cast<int>(timing.objective_gradient_count)); }
            { VirtualTimer<CountedTimingNode, int> timer_wall_g_eq("Equality Constraints", timing.constraints_nano, static_cast<int>(timing.constraints_count)); }
            { VirtualTimer<CountedTimingNode, int> timer_wall_jac_g_eq("Jacobian Equality Constraints", timing.jacobian_nano, static_cast<int>(timing.jacobian_count)); }
            { VirtualTimer<CountedTimingNode, int> timer_wall_lag_hes("Lagrangian Hessian", timing.hessian_nano, static_cast<int>(timing.hessian_count)); }
        }
    }
}

} // namespace UnoSolver
