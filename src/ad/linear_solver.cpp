// SPDX-License-Identifier: LGPL-3.0-or-later
#include "linear_solver.h"

#include "detail/matrix_ops.h"

#include <stdexcept>
#include <string>
#include <vector>

extern "C" {
void dgesv_(const int *n,
            const int *nrhs,
            double *a,
            const int *lda,
            int *ipiv,
            double *b,
            const int *ldb,
            int *info);
}

namespace ad {
namespace {

double matrix_value(const LinearSolveProblem &problem, int row, int col) {
    return problem.matrix[detail::matrix_flat_index(row, col, problem.size, problem.size, problem.layout)];
}

void validate_problem(const LinearSolveProblem &problem) {
    if (problem.size < 0) {
        throw std::runtime_error("linear solve size must be non-negative");
    }
    if (problem.size > 0 && problem.matrix == nullptr) {
        throw std::runtime_error("linear solve matrix pointer is null");
    }
    if (problem.size > 0 && problem.rhs == nullptr) {
        throw std::runtime_error("linear solve rhs pointer is null");
    }
    if (problem.size > 0 && problem.solution == nullptr) {
        throw std::runtime_error("linear solve solution pointer is null");
    }
}

void lapack_solve_checked(const LinearSolveProblem &problem) {
    validate_problem(problem);
    if (problem.size == 0) {
        return;
    }

    const int n = problem.size;
    const int nrhs = 1;
    const int lda = n;
    const int ldb = n;
    int info = 0;
    std::vector<double> a(static_cast<std::size_t>(n * n), 0.0);
    std::vector<int> ipiv(static_cast<std::size_t>(n), 0);

    for (int col = 0; col < n; ++col) {
        for (int row = 0; row < n; ++row) {
            a[static_cast<std::size_t>(row + col * n)] = problem.transpose
                ? matrix_value(problem, col, row)
                : matrix_value(problem, row, col);
        }
    }
    for (int i = 0; i < n; ++i) {
        problem.solution[i] = problem.rhs[i];
    }

    dgesv_(&n, &nrhs, a.data(), &lda, ipiv.data(), problem.solution, &ldb, &info);
    if (info < 0) {
        throw std::runtime_error("LAPACK dgesv argument " + std::to_string(-info) + " had an illegal value");
    }
    if (info > 0) {
        throw std::runtime_error("LAPACK dgesv found a singular matrix");
    }
}

class LapackLinearSolver final : public LinearSolver {
public:
    const char *name() const override {
        return "lapack_lu";
    }

    void solve(const LinearSolveProblem &problem) const override {
        lapack_solve_checked(problem);
    }
};

} // namespace

const LinearSolver &lapack_linear_solver() {
    static const LapackLinearSolver solver;
    return solver;
}

const LinearSolver &linear_solver(LinearSolverKind kind) {
    switch (kind) {
        case LinearSolverKind::LapackLU:
            return lapack_linear_solver();
    }
    throw std::runtime_error("unsupported linear solver kind");
}

} // namespace ad

extern "C" int moo_ad_lapack_linear_solve(const MOOAdLinearSolveProblem *problem) {
    if (problem == nullptr) {
        return -1000;
    }
    try {
        ad::LinearSolveProblem cpp_problem;
        cpp_problem.size = problem->size;
        cpp_problem.layout = problem->layout == 0 ? ad::MatrixLayout::ColumnMajor : ad::MatrixLayout::RowMajor;
        cpp_problem.transpose = problem->transpose != 0;
        cpp_problem.matrix = problem->matrix;
        cpp_problem.rhs = problem->rhs;
        cpp_problem.solution = problem->solution;
        ad::lapack_linear_solver().solve(cpp_problem);
        return 0;
    } catch (const std::runtime_error &) {
        return -1001;
    } catch (...) {
        return -1002;
    }
}
