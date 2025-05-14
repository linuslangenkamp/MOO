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
void set_parameters(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* p);
void set_states(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* x_ij);
void set_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* u_ij);
void set_states_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* xu_ij);
void set_time(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64 t_ij);

/* TODO: is there a better combination e.g. functionODE + someother? */
/* evaluation at current point */
inline void eval_current_point(DATA* data, threadData_t* threadData, InfoGDOP& info) {
    data->callback->functionDAE(data, threadData);
}

/* write previous evaluation to buffer */
void eval_lfg_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_lfg_buffer);
void eval_mr_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_mr_buffer);

/* eval jacobian and write to buffer in *CSC* form; just passes the current buffer with offset to OM Jacobian */
inline void jac_eval_write_csc_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, JACOBIAN* jacobian, F64* eval_jac_buffer) {
    assert(jacobian != NULL && jacobian->sparsePattern != NULL);
    new_evalJacobian(data, threadData, jacobian, NULL, eval_jac_buffer, new_JACOBIAN_OUTPUT_FORMAT::new_JAC_OUTPUT_CSC);
}

void jac_CD_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_jac_mr_buffer);

#endif // OPT_OM_EVALUATIONS_H
