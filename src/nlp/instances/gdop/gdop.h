#ifndef OPT_GDOP_H
#define OPT_GDOP_H

#include <cassert>

#include <base/block_sparsity.h>
#include <base/collocation.h>
#include <base/trajectory.h>
#include <base/fixed_vector.h>
#include <base/linalg.h>
#include <base/nlp_structs.h>
#include <base/mesh.h>
#include <base/util.h>
#include <nlp/nlp.h>

#include "problem.h"
#include "gdop_strategies.h"

namespace GDOP {

class GDOP : public NLP::NLP {
public:
  GDOP(Problem& problem,
        Collocation& collocation,
        Mesh& mesh)
      : NLP::NLP(),
        mesh(mesh),
        problem(problem),
        collocation(collocation) {}
  // TODO: what can be made private, also make layout (args) nice like in Ipopt

  // structures
  Mesh& mesh;                 // grid / mesh
  Problem& problem;           // continuous GDOP
  Collocation& collocation;   // collocation data
  NLP_State evaluation_state; // simple state to check which callbacks are performed for an iteration

  // initial guess
  std::unique_ptr<PrimalDualTrajectory> initial_guess; // set by initialization strategy

  // optimal solution
  std::unique_ptr<PrimalDualTrajectory> optimal_solution;

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
  BlockSparsity hes_a, hes_b, hes_c, hes_d, hes_e, hes_f, hes_g, hes_h;
  OrderedIndexSet A, B, C, D, E, F, G, H; // rename

  void update(Mesh&& new_mesh);

  void set_initial_guess(std::unique_ptr<PrimalDualTrajectory> initial_trajectory);

  // virtuals in NLP
  void get_sizes(
    int& number_vars,
    int& number_constraints) override;

  void get_bounds(
    FixedVector<f64>& x_lb,
    FixedVector<f64>& x_ub,
    FixedVector<f64>& g_lb,
    FixedVector<f64>& g_ub) override;

  void get_initial_guess(
    bool init_x,
    FixedVector<f64>& x_init,
    bool init_lambda,
    FixedVector<f64>& lambda_init,
    bool init_z,
    FixedVector<f64>& z_lb_init,
    FixedVector<f64>& z_ub_init) override;

  void GDOP::get_nnz(
    int& nnz_jac,
    int& nnz_hes) override;

  void GDOP::get_jac_sparsity(
      FixedVector<int> i_row_jac,
      FixedVector<int> j_col_jac);

  void GDOP::get_hes_sparsity(
      FixedVector<int> i_row_hes,
      FixedVector<int> j_col_hes);

  void eval_f(
    bool new_x,
    const FixedVector<f64>& curr_x,
    f64& curr_obj) override;

  void eval_g(
    bool new_x,
    const FixedVector<f64>& curr_x,
    FixedVector<f64>& curr_g) override;

  void eval_grad_f(
    bool new_x,
    const FixedVector<f64>& curr_x,
    FixedVector<f64>& curr_grad_f) override;

  void eval_jac_g(
    bool new_x,
    const FixedVector<f64>& curr_x,
    const FixedVector<int>& i_row_jac,
    const FixedVector<int>& j_col_jac,
    FixedVector<f64>& curr_jac) override;

  void eval_hes(
    bool new_x,
    const FixedVector<f64>& curr_x,
    bool new_lambda,
    FixedVector<f64>& curr_lambda,
    f64& curr_obj_factor,
    const FixedVector<int>& i_row_hes,
    const FixedVector<int>& j_col_hes,
    FixedVector<f64>& curr_hes) override;

  void finalize_solution() override;


private:
  // inline methods for getting and providing current variable / dual addresses in callback
  // x0 => x(t0), xu => xu(t_01), xuf => xu(t_f), p => p, lamb_fg => fg(t_01), lamb_r => r
  inline const f64* get_curr_x_x0(const FixedVector<f64>& curr_x)         { return off_x        != 0 ? curr_x.raw()          : nullptr; }
  inline const f64* get_curr_x_xu(const FixedVector<f64>& curr_x)         { return off_xu       != 0 ? &curr_x[off_x]        : nullptr; }
  inline const f64* get_curr_x_xuf(const FixedVector<f64>& curr_x)        { return off_xu       != 0 ? &curr_x[off_last_xu]  : nullptr; }
  inline const f64* get_curr_x_p(const FixedVector<f64>& curr_x)          { return off_p        != 0 ? &curr_x[off_xu_total] : nullptr; }
  inline f64* get_curr_lamb_fg(FixedVector<f64>& curr_lambda) { return off_fg_total != 0 ? curr_lambda.raw()     : nullptr; }
  inline f64* get_curr_lamb_r(FixedVector<f64>& curr_lambda)  { return problem.boundary->r_size != 0 ? &curr_lambda[off_fg_total] : nullptr; }

  // helpers for initialize offsets 
  void create_acc_offset_xu(int off_x, int off_xu);
  void create_acc_offset_fg(int off_fg);

  // init nlp and sparsity
  void init_jacobian_nonzeros(int& nnz_jac);
  void init_hessian_nonzeros(int& nnz_hes);

  /* mutiply lambda (dual) with mesh factors => callbacks (except Lagrange) can use exact multipliers */
  void update_curr_lambda_obj_factors(FixedVector<f64>& curr_lambda, f64 curr_sigma_f);

  /* augmented hessian updates */
  void update_augmented_hessian_lfg(const AugmentedHessianLFG& hes, const int i, const int j,
                                    const BlockSparsity* ptr_map_xu_xu, const BlockSparsity* ptr_map_p_xu, FixedVector<f64>& curr_hes);
  void update_augmented_parameter_hessian_lfg(const AugmentedParameterHessian& aug_hes, FixedVector<f64>& curr_hes); // sum of all weighted Hessian(Lfg)_pp
  void update_augmented_hessian_mr(const AugmentedHessianMR& hes, FixedVector<f64>& curr_hes);

  // get callback data
  void callback_evaluation(const FixedVector<f64>& curr_x);
  void callback_jacobian(const FixedVector<f64>& curr_x);
  void callback_hessian(const FixedVector<f64> x, FixedVector<f64>& curr_lambda, f64 curr_sigma_f);

  // internal evaluations
  void check_new_x(bool new_x);
  void check_new_lambda(bool new_lambda);
  void eval_f_internal(f64& curr_obj);
  void eval_g_internal(const FixedVector<f64>& curr_x, FixedVector<f64>& curr_g);
  void eval_grad_f_internal(FixedVector<f64>& curr_grad);
  void eval_jac_g_internal(FixedVector<f64>& curr_jac);
  void eval_hes_internal(FixedVector<f64>& curr_hes);


  void flatten_trajectory_to_layout(const Trajectory& Trajectory, FixedVector<f64>& flat_buffer);
  void transform_duals_costates(FixedVector<f64>& lambda, bool to_costate);
  void transform_duals_costates_bounds(FixedVector<f64>& zeta, bool to_costate);

  std::unique_ptr<Trajectory> finalize_optimal_primals();
  std::unique_ptr<CostateTrajectory> finalize_optimal_costates();
  std::pair<std::unique_ptr<Trajectory>, std::unique_ptr<Trajectory>> finalize_optimal_bound_duals();
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