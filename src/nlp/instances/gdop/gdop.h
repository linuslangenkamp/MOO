#ifndef OPT_GDOP_H
#define OPT_GDOP_H

#include <cassert>

#include <base/block_sparsity.h>
#include <base/collocation.h>
#include <base/fixed_vector.h>
#include <base/linalg.h>
#include <base/nlp_state.h>
#include <base/nlp_structs.h>
#include <base/mesh.h>
#include <base/util.h>
#include <nlp/nlp.h>

#include "problem.h"


class GDOP : public NLP {
public:
    GDOP(std::shared_ptr<Problem> problem, std::unique_ptr<Collocation> collocation, std::shared_ptr<Mesh> mesh, std::shared_ptr<Trajectory> guess)
        : NLP(),
          mesh(mesh),
          problem(problem),
          collocation(std::move(collocation)),
          guess(guess) {
        init();
    }

    // structures
    std::shared_ptr<Mesh> mesh;                // grid / mesh
    std::shared_ptr<Problem> problem;          // continuous GDOP
    std::unique_ptr<Collocation> collocation;  // collocation data
    std::shared_ptr<Trajectory> guess;         // initial guess / trajectory, will be interpolated accordingly
    NLP_State evaluation_state;                // simple state to check which actions are / have to be performed for an iteration

    // constant NLP derivative matrix part of the jacobian
    FixedVector<double> const_der_jac;

    // offsets
    int off_x;         // offset #xVars
    int off_u;         // offset #uVars
    int off_p;         // offset #pVars
    int off_xu;        // number of vars for one collocation grid point
    int off_last_xu;   // last collocation grid point *x_nm, u_nm
    int off_xu_total;  // first parameter variable index
    int off_fg_total;  // first boundary constraint index

    // note off_acc_xu[0][0] = off_x for time t = first collocation node, since there are no controls at time t=0
    FixedField<int, 2>   off_acc_xu;  // offset to NLP_X first index of (x, u)(t_ij), i.e. NLP_X[off_acc_xu[i][j]] = x[i][j], u[i][j]
    FixedField<int, 2>   off_acc_fg;  // offset to NLP_G first index of (f, g)(t_ij), i.e. NLP_G[off_acc_fg[i][j]] = f[i][j], g[i][j]
    FixedVector<int> off_acc_jac_fg;  // offset to NLP_JAC_G first index of nabla (f, g)(t_ij)

    // hessian sparsity helpers, O(1/2 * (x + u)² + p * (p + x + u)) memory, but no need for hashmaps, these are still fairly cheap
    // for further info see hessian layout at the bottom
    BlockSparsity hes_a = BlockSparsity::createLowerTriangular(problem->x_size, BlockType::Exact);
    BlockSparsity hes_b = BlockSparsity::createLowerTriangular(problem->x_size + problem->u_size, BlockType::Offset);
    BlockSparsity hes_c = BlockSparsity::createSquare(problem->x_size, BlockType::Exact);
    BlockSparsity hes_d = BlockSparsity::createLowerTriangular(problem->x_size + problem->u_size, BlockType::Exact);
    BlockSparsity hes_e = BlockSparsity::createRectangular(problem->p_size, problem->x_size, BlockType::Exact);
    BlockSparsity hes_f = BlockSparsity::createRectangular(problem->p_size, problem->x_size + problem->u_size, BlockType::RowOffset);
    BlockSparsity hes_g = BlockSparsity::createRectangular(problem->p_size, problem->x_size + problem->u_size, BlockType::Exact);
    BlockSparsity hes_h = BlockSparsity::createLowerTriangular(problem->p_size, BlockType::Exact);

    // init nlp and sparsity
    void init();
    void initSizesOffsets();
    void initBuffers();
    void initBounds();
    void initStartingPoint();
    void initJacobian();
    void initJacobianNonzeros();
    void initJacobianSparsityPattern();
    void initHessian();

    // hessian updates
    void updateHessianLFG(FixedVector<double>& values, const HessianLFG& hes, const int i, const int j, const BlockSparsity* ptr_map_xu_xu,
                          const BlockSparsity* ptr_map_p_xu, const double factor);
    void updateHessianMR(FixedVector<double>& values, const HessianMR& hes, const double factor);

    // inline methods to jump (i, j) callback buffer blocks
    int jac_offset(int i, int j);
    int hes_offset(int i, int j);

    // get callback data
    void callback_evaluation();
    void callback_jacobian();
    void callback_hessian();

    // nlp solver calls
    void check_new_x(const double* nlp_solver_x, bool new_x);
    void check_new_lambda(const double* nlp_solver_lambda, const bool new_lambda);
    void check_new_sigma(const double obj_factor);
    void eval_f_internal();
    void eval_g_internal();
    void eval_grad_f_internal();
    void eval_jac_g_internal();
    void eval_hes_internal();

    // virtuals in NLP
    void eval_f(const double* nlp_solver_x, bool new_x);
    void eval_g(const double* nlp_solver_x, bool new_x);
    void eval_grad_f(const double* nlp_solver_x, bool new_x);
    void eval_jac_g(const double* nlp_solver_x, bool new_x);
    void eval_hes(const double* nlp_solver_x, const double* nlp_solver_lambda, double sigma, bool new_x, bool new_lambda);
};

#endif  // OPT_GDOP_H

/*
Hessian Sparsity Layout (lower triangle):
    L: lower triangular matrix with diagonal
    X: square / rectangular matrix
    Note that blocks [[L, 0], [X, L]] are also triangular

     | x00 | x01 u01 | x02 u02 | x** u** | xnm1 unm1 | xnm unm | p |
-------------------------------------------------------------------|
 x00 |  L  |         |         |         |           |         |   |
-------------------------------------------------------------------|
 x01 |     |  L   0  |         |         |           |         |   |
 u01 |     |  X   L  |         |         |           |         |   |
 ------------------------------------------------------------------|
 x02 |     |         |  L   0  |         |           |         |   |
 u02 |     |         |  X   L  |         |           |         |   |
 ------------------------------------------------------------------|
 x** |     |         |         |  L   0  |           |         |   |
 u** |     |         |         |  X   L  |           |         |   |
 ------------------------------------------------------------------|
 xnm1|     |         |         |         |   L   0   |         |   |
 unm1|     |         |         |         |   X   L   |         |   |
-------------------------------------------------------------------|
 xnm |  X  |         |         |         |           |  L   0  |   |
 unm |     |         |         |         |           |  X   L  |   |
-------------------------------------------------------------------|
  p  |  X  |  X   X  |  X   X  |  X   X  |   X   X   |  X   X  | L |
-------------------------------------------------------------------*

Block Sparsity Patterns: A - H
where A=triang(x) B=triang(x + u), C=sq(x), D=triang(x + u), E=rect(p, x),
      F=rect(p, x + u), G=rect(p, x + u), H=triang(p, p)

     | x00 | x01 u01 | x02 u02 | x** u** | xnm1 unm1 | xnm unm | p |
-------------------------------------------------------------------|
 x00 |  A  |         |         |         |           |         |   |
-------------------------------------------------------------------|
 x01 |     |    B    |         |         |           |         |   |
 u01 |     |         |         |         |           |         |   |
 ------------------------------------------------------------------|
 x02 |     |         |    B    |         |           |         |   |
 u02 |     |         |         |         |           |         |   |
 ------------------------------------------------------------------|
 x** |     |         |         |    B    |           |         |   |
 u** |     |         |         |         |           |         |   |
 ------------------------------------------------------------------|
 xnm1|     |         |         |         |     B     |         |   |
 unm1|     |         |         |         |           |         |   |
-------------------------------------------------------------------|
 xnm |  C  |         |         |         |           |    D    |   |
 unm |     |         |         |         |           |         |   |
-------------------------------------------------------------------|
  p  |  E  |    F    |    F    |    F    |     F     |    G    | H |
-------------------------------------------------------------------*
*/