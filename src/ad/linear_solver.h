// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_LINEAR_SOLVER_H
#define MOO_AD_LINEAR_SOLVER_H

#include "matrix.h"

namespace ad {

struct LinearSolveProblem {
    int size = 0;
    MatrixLayout layout = MatrixLayout::ColumnMajor;
    bool transpose = false;
    const double *matrix = nullptr;
    const double *rhs = nullptr;
    double *solution = nullptr;
};

class LinearSolver {
public:
    virtual ~LinearSolver() = default;
    virtual const char *name() const = 0;
    virtual void solve(const LinearSolveProblem &problem) const = 0;
};

const LinearSolver &lapack_linear_solver();
const LinearSolver &linear_solver(LinearSolverKind kind);

} // namespace ad

extern "C" {

struct MOOAdLinearSolveProblem {
    int size;
    int layout;
    int transpose;
    const double *matrix;
    const double *rhs;
    double *solution;
};

int moo_ad_lapack_linear_solve(const MOOAdLinearSolveProblem *problem);

} // extern "C"

#endif // MOO_AD_LINEAR_SOLVER_H
