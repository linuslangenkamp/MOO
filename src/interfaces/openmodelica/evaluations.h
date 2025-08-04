#ifndef OPT_OM_EVALUATIONS_H
#define OPT_OM_EVALUATIONS_H

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>
#include <base/block_sparsity.h>

#include "debug_om.h"
#include "info_gdop.h"
#include "sim_runtime_ext.h"

namespace OpenModelica {

// TODO: try to not call functionDAE that often. Maybe create a workspace buffer where we memcpy the realVars after evaluation
//       main advantage: no need to solve NLS several times!!
// => maybe allocate an additional buffer for the realVars + parameters and memcpy them inside the evaluation again

/* init evaluations */
void init_eval(InfoGDOP& info, GDOP::FullSweepLayout& layout_lfg, GDOP::BoundarySweepLayout& layout_mr);
void init_eval_lfg(InfoGDOP& info, GDOP::FullSweepLayout& layout_lfg);
void init_eval_mr(InfoGDOP& info, GDOP::BoundarySweepLayout& layout_mr);

/* init Jacobians  */
void init_jac(InfoGDOP& info, GDOP::FullSweepLayout& layout_lfg, GDOP::BoundarySweepLayout& layout_mr);
void init_jac_lfg(InfoGDOP& info, GDOP::FullSweepLayout& layout_lfg);
void init_jac_mr(InfoGDOP& info, GDOP::BoundarySweepLayout& layout_mr);

/* init Hessians */
void init_hes(InfoGDOP& info, GDOP::FullSweepLayout& layout_lfg, GDOP::BoundarySweepLayout& layout_mr);
void init_hes_lfg(InfoGDOP& info, GDOP::FullSweepLayout& layout_lfg);
void init_hes_mr(InfoGDOP& info, GDOP::BoundarySweepLayout& layout_mr);

/* set values in OM realVars array / time value */
void set_parameters(InfoGDOP& info, const f64* p);
void set_states(InfoGDOP& info, const f64* x_ij);
void set_inputs(InfoGDOP& info, const f64* u_ij);
void set_states_inputs(InfoGDOP& info, const f64* xu_ij);
void set_time(InfoGDOP& info, const f64 t_ij);

inline void eval_current_point(InfoGDOP& info) {
    info.data->callback->functionDAE(info.data, info.threadData);
}

/* write previous evaluation to buffer */
void eval_lfg_write(InfoGDOP& info, f64* eval_lfg_buffer);
void eval_mr_write(InfoGDOP& info, f64* eval_mr_buffer);

/* call evalJacobian and write to buffer in *CSC* form; just passes the current buffer with offset to OM Jacobian */
inline void jac_eval_write_as_csc(InfoGDOP& info, JACOBIAN* jacobian, f64* eval_jac_buffer) {
    assert(jacobian != NULL && jacobian->sparsePattern != NULL);
    mooEvalJacobian(info.data, info.threadData, jacobian, NULL, eval_jac_buffer);
}

/* eval full jacobian (full_buffer) but only fill eval_jac_buffer with elements of first, *moved* row in Exchange's COO structure
 * (see construction of CscToCoo structures for further info)
 * clearly order doesnt matter for one row; CSC == COO order for this row, but the entries in original CSC are
 * stored in CscToCoo.coo_to_csc(nz). */
void jac_eval_write_first_row_as_csc(InfoGDOP& info, JACOBIAN* jacobian, f64* full_buffer,
                                     f64* eval_jac_buffer, CscToCoo& exc);

} // namespace OpenModelica

#endif // OPT_OM_EVALUATIONS_H
