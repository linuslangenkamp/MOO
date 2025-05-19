#ifndef OPT_OM_EVALUATIONS_H
#define OPT_OM_EVALUATIONS_H

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "debug_om.h"
#include "info_gdop.h"
#include "sim_runtime_ext.h"

/* init evaluations */
void init_eval(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr);
void init_eval_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg);
void init_eval_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr);

/* init Jacobians  */
void init_jac(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr);
void init_jac_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg);
void init_jac_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr);

/* set values in OM realVars array / time value */
void set_parameters(DATA* data, threadData_t* threadData, InfoGDOP& info, const f64* p);
void set_states(DATA* data, threadData_t* threadData, InfoGDOP& info, const f64* x_ij);
void set_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const f64* u_ij);
void set_states_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const f64* xu_ij);
void set_time(DATA* data, threadData_t* threadData, InfoGDOP& info, const f64 t_ij);

/* TODO: is there a better combination e.g. functionODE + someother? */
/* evaluation at current point */
inline void eval_current_point(DATA* data, threadData_t* threadData, InfoGDOP& info) {
    data->callback->functionDAE(data, threadData);
}

/* write previous evaluation to buffer */
void eval_lfg_write(DATA* data, threadData_t* threadData, InfoGDOP& info, f64* eval_lfg_buffer);
void eval_mr_write(DATA* data, threadData_t* threadData, InfoGDOP& info, f64* eval_mr_buffer);

/* call evalJacobian and write to buffer in *CSC* form; just passes the current buffer with offset to OM Jacobian */
inline void jac_eval_write_as_csc(DATA* data, threadData_t* threadData, InfoGDOP& info, JACOBIAN* jacobian, f64* eval_jac_buffer) {
    assert(jacobian != NULL && jacobian->sparsePattern != NULL);
    __evalJacobian(data, threadData, jacobian, NULL, eval_jac_buffer);
}

/* eval full jacobian (full_buffer) but only fill eval_jac_buffer with elements of first, *moved* row in Exchange's COO structure
 * (see construction of Exchange_COO_CSC structures for further info)
 * clearly order doesnt matter for one row; CSC == COO order for this row, but the entries in original CSC are
 * stored in Exchange_COO_CSC.coo_to_csc(nz). */
void jac_eval_write_first_row_as_csc(DATA* data, threadData_t* threadData, InfoGDOP& info, JACOBIAN* jacobian, f64* full_buffer,
                                     f64* eval_jac_buffer, Exchange_COO_CSC& exc);

#endif // OPT_OM_EVALUATIONS_H
