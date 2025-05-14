#include "evaluations.h"

/* just enumerate them from 0 ... #lfg - 1 abd 0 ... #mr - 1, the correct placement will be handled in eval */
void init_eval(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr) {
    init_eval_lfg(data, threadData, info, lfg);
    init_eval_mr(data, threadData, info, mr);
}

void init_eval_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg) {
    for (int i = 0; i < lfg.int_size(); i++) {
        lfg[i].buf_index = i;
    }
}

void init_eval_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr) {
    for (int i = 0; i < mr.int_size(); i++) {
        mr[i].buf_index = i;
    }
}

/* since its no problem and conversion from COO <-> CSC has been carried out, we just use the CSC ordering in the Jacobians
 * thus hopefully only using a memcpy to the OPT buffer for Lfg and 2 memcpy's for C(Mayer CSC) and D(r CSC) */
void init_jac(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr) {
    init_jac_lfg(data, threadData, info, lfg);
    init_jac_mr(data, threadData, info, mr);
}

void init_jac_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg) {
    /* full B Jacobian */
    for (int nz = 0; nz < info.exc_jac->B_coo.nnz; nz++) {
        int row = info.exc_jac->B_coo.row[nz];
        int col = info.exc_jac->B_coo.col[nz];
        int csc_buffer_entry_B = info.exc_jac->B_coo.coo_to_csc(nz); // jac_buffer == OpenModelica B CSC buffer!
        if (col < info.x_size) {
            lfg[row].jac.dx.push_back(JacobianSparsity{col, csc_buffer_entry_B});
        }
        else if (col < info.xu_size) {
            lfg[row].jac.du.push_back(JacobianSparsity{col - info.x_size, csc_buffer_entry_B});
        }
        else {
            lfg[row].jac.dp.push_back(JacobianSparsity{col - info.xu_size, csc_buffer_entry_B});
        }
    }
}

void init_jac_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr) {
    /* M (first row) in C(COO) Jacobian */
    int nz_C = 0;
    if (info.mayer_exists) {
        while (info.exc_jac->C_coo.row[nz_C] == 0) {
            int col = info.exc_jac->C_coo.col[nz_C];
            int csc_buffer_entry_C = info.exc_jac->C_coo.coo_to_csc(nz_C); // jac_buffer == OpenModelica C(row=0) CSC buffer!

            /* for now only final constraints! no parameters, no dx0, esp. no du! */
            if (col < info.x_size) {
                mr[0].jac.dxf.push_back(JacobianSparsity{col, csc_buffer_entry_C});
            }

            nz_C++;
        }
    }

    /* r in D Jacobian */
    int r_start = (int)info.mayer_exists;
    for (int nz_D = 0; nz_D < info.exc_jac->D_coo.nnz; nz_D++) {
        int row = info.exc_jac->D_coo.row[nz_D];
        int col = info.exc_jac->D_coo.col[nz_D];
        int csc_buffer_entry_D = info.exc_jac->D_coo.coo_to_csc(nz_D); // jac_buffer == OpenModelica D CSC buffer!
        if (col < info.x_size) {
            /* add the Mayer offset, since the values F64* is [M, r] */
            mr[r_start + row].jac.dxf.push_back(JacobianSparsity{col, info.exc_jac->D_coo.nnz_offset + csc_buffer_entry_D});
        }
    }
}

/* TODO: add me */
void set_parameters(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* p) {
    return;
}

void set_states(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* x_ij) {
    for (int x = 0; x < info.x_size; x++) {
        data->localData[0]->realVars[info.index_x_real_vars + x] = (modelica_real) x_ij[x];
    }
}

void set_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* u_ij) {
    for (int u = 0; u < info.u_size; u++) {
        data->localData[0]->realVars[info.u_indices_real_vars[u]] = (modelica_real) u_ij[info.x_size + u];
    }
}

void set_states_inputs(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64* xu_ij) {
    for (int x = 0; x < info.x_size; x++) {
        data->localData[0]->realVars[info.index_x_real_vars + x] = (modelica_real) xu_ij[x];
    }

    for (int u = 0; u < info.u_size; u++) {
        data->localData[0]->realVars[info.u_indices_real_vars[u]] = (modelica_real) xu_ij[info.x_size + u];
    }
}

void set_time(DATA* data, threadData_t* threadData, InfoGDOP& info, const F64 t_ij) {
    data->localData[0]->timeValue = (modelica_real) t_ij;
}

void eval_lfg_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_lfg_buffer) {
    int nz = 0;
    /* L */
    if (info.lagrange_exists) {
        eval_lfg_buffer[nz++] = data->localData[0]->realVars[info.index_lagrange_real_vars];
    }
    /* f */
    for (int der_x = 0; der_x < info.f_size; der_x++) {
        eval_lfg_buffer[nz++] = data->localData[0]->realVars[info.index_der_x_real_vars + der_x];
    }
    /* g */
    for (int g = 0; g < info.g_size; g++) {
        eval_lfg_buffer[nz++] = data->localData[0]->realVars[info.index_g_real_vars + g];
    }
}

void eval_mr_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, F64* eval_mr_buffer) {
    int nz = 0;
    /* M */
    if (info.mayer_exists) {
        eval_mr_buffer[nz++] = data->localData[0]->realVars[info.index_mayer_real_vars];
    }
    /* r */
    for (int r = 0; r < info.r_size; r++) {
        eval_mr_buffer[nz++] = data->localData[0]->realVars[info.index_r_real_vars + r];
    }
}

void jac_CD_write_to_buffer(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<F64>& eval_jac_mr_buffer) {

}
