#ifndef OPT_GDOP_H
#define OPT_GDOP_H

#include <cassert>

#include <base/block_sparsity.h>
#include <base/collocation.h>
#include <base/fixed_vector.h>
#include <base/linalg.h>
#include <base/nlp_structs.h>
#include <base/mesh.h>
#include <base/util.h>
#include <nlp/nlp.h>

#include "problem.h"
#include "gdop_strategies.h"

namespace GDOP {

// TODO: think about this. Should we also delay the strategies? Maybe we want to use different strategies from time to time?!
class GDOP : public NLP {
public:
    GDOP(Problem& problem, Collocation& collocation, Mesh& mesh, std::unique_ptr<Strategies> strategies)
        : NLP(),
          mesh(mesh),
          problem(problem),
          collocation(collocation) {

        if (strategies != nullptr) {
            this->strategies = std::move(strategies);
        }
        else {
            this->strategies = std::make_unique<Strategies>(Strategies::default_strategies());
        }

        init();
    }

    // structures
    Mesh& mesh;                 // grid / mesh
    Problem& problem;           // continuous GDOP
    Collocation& collocation;   // collocation data
    NLP_State evaluation_state; // simple state to check which callbacks are performed for an iteration

    // strategies for initialization, simulation, mesh refinement etc.
    std::unique_ptr<Strategies> strategies;

    std::unique_ptr<Trajectory> optimal_solution = std::make_unique<Trajectory>(); // gets filled in finalize_solution()

    // constant NLP derivative matrix part of the jacobian
    FixedVector<f64> const_der_jac;

    // offsets
    int off_x;        // offset #xVars
    int off_u;        // offset #uVars
    int off_p;        // offset #pVars
    int off_xu;       // number of vars for one collocation grid point
    int off_last_xu;  // last collocation grid point *x_nm, u_nm
    int off_xu_total; // first parameter variable index
    int off_fg_total; // first boundary constraint index

    // note off_acc_xu[0][0] = off_x for time t = first collocation node, since there are no controls at time t=0
    FixedField<int, 2>   off_acc_xu; // offset to NLP_X first index of (x, u)(t_ij), i.e. NLP_X[off_acc_xu[i][j]] = x[i][j], u[i][j]
    FixedField<int, 2>   off_acc_fg; // offset to NLP_G first index of (f, g)(t_ij), i.e. NLP_G[off_acc_fg[i][j]] = f[i][j], g[i][j]
    FixedVector<int> off_acc_jac_fg; // offset to NLP_JAC_G first index of nabla (f, g)(t_ij)

    /* scaling factors for lagrange terms in augmented Hessian callback */
    FixedField<f64, 2> lagrange_obj_factors; // = sigma_f * collocation.b[mesh.intervals[i]][mesh.nodes[j]] mesh.delta_t[i]

    // hessian sparsity helpers, O(1/2 * (x + u)² + p * (p + x + u)) memory, but no need for hashmaps, these are still fairly cheap
    // for further info see hessian layout at the bottom
    BlockSparsity hes_a = BlockSparsity::create_lower_triangular(problem.x_size, BlockType::Exact);
    BlockSparsity hes_b = BlockSparsity::create_lower_triangular(problem.x_size + problem.u_size, BlockType::Offset);
    BlockSparsity hes_c = BlockSparsity::create_rectangular(problem.x_size + problem.u_size, problem.x_size, BlockType::Exact);
    BlockSparsity hes_d = BlockSparsity::create_lower_triangular(problem.x_size + problem.u_size, BlockType::Exact);
    BlockSparsity hes_e = BlockSparsity::create_rectangular(problem.p_size, problem.x_size, BlockType::Exact);
    BlockSparsity hes_f = BlockSparsity::create_rectangular(problem.p_size, problem.x_size + problem.u_size, BlockType::RowOffset);
    BlockSparsity hes_g = BlockSparsity::create_rectangular(problem.p_size, problem.x_size + problem.u_size, BlockType::Exact);
    BlockSparsity hes_h = BlockSparsity::create_lower_triangular(problem.p_size, BlockType::Exact);

    // init nlp and sparsity
    void init();
    void init_sizes_offsets();
    void init_buffers();
    void init_bounds();
    void init_starting_point();
    void init_jacobian();
    void init_jacobian_nonzeros();
    void init_jacobian_sparsity_pattern();
    void init_hessian();

    /* mutiply lambda (dual) with mesh factors => callbacks (except Lagrange) can use exact multipliers */
    void update_curr_lambda_obj_factors();

    /* augmented hessian updates */
    void update_augmented_hessian_lfg(const AugmentedHessianLFG& hes, const int i, const int j,
                                      const BlockSparsity* ptr_map_xu_xu, const BlockSparsity* ptr_map_p_xu);
    void update_augmented_parameter_hessian_lfg(const AugmentedParameterHessian& aug_hes); // sum of all weighted Hessian(Lfg)_pp
    void update_augmented_hessian_mr(const AugmentedHessianMR& hes);

    // get callback data
    void callback_evaluation();
    void callback_jacobian();
    void callback_hessian();

    // inline methods for getting and providing current variable / dual addresses in callback
    // x0 => x(t0), xu => xu(t_01), xuf => xu(t_f), p => p, lamb_fg => fg(t_01), lamb_r => r
    inline f64* get_curr_x_x0()    { return off_x        != 0 ? curr_x.raw()          : nullptr; }
    inline f64* get_curr_x_xu()    { return off_xu       != 0 ? &curr_x[off_x]        : nullptr; }
    inline f64* get_curr_x_xuf()   { return off_xu       != 0 ? &curr_x[off_last_xu]  : nullptr; }
    inline f64* get_curr_x_p()     { return off_p        != 0 ? &curr_x[off_xu_total] : nullptr; }
    inline f64* get_curr_lamb_fg() { return off_fg_total != 0 ? curr_lambda.raw()     : nullptr; }
    inline f64* get_curr_lamb_r()  { return problem.boundary->r_size != 0 ? &curr_lambda[off_fg_total] : nullptr; }

    // nlp solver calls
    void check_new_x(bool new_x);
    void check_new_lambda(const bool new_lambda);
    void eval_f_internal();
    void eval_g_internal();
    void eval_grad_f_internal();
    void eval_jac_g_internal();
    void eval_hes_internal();

    // virtuals in NLP
    void eval_f(bool new_x);
    void eval_g(bool new_x);
    void eval_grad_f(bool new_x);
    void eval_jac_g(bool new_x);
    void eval_hes(bool new_x, bool new_lambda);
    void finalize_solution();
};

} // namespace GDOP

#endif // OPT_GDOP_H

/*
Hessian Sparsity Layout (lower triangle):
    L: lower triangular matrix with diagonal
    X: square / rectangular matrix
    Note that blocks [[L, 0], [X, L]] are also triangular

                                            {n,m-1}
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
 unm |  X  |         |         |         |           |  X   L  |   |
-------------------------------------------------------------------|
  p  |  X  |  X   X  |  X   X  |  X   X  |   X   X   |  X   X  | L |
-------------------------------------------------------------------*

Block Sparsity Patterns: A - H
where A=triang(x) B=triang(x + u), C=sq(x), D=triang(x + u), E=rect(p, x),
      F=rect(p, x + u), G=rect(p, x + u), H=triang(p, p)
                                            {n,m-1}
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