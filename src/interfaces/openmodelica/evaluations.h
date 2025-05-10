#ifndef OPT_OM_EVALUATIONS
#define OPT_OM_EVALUATIONS

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "helper.h"

struct ExchangeJacobians {
    JACOBIAN* A;
    JACOBIAN* B;
    JACOBIAN* C;
    JACOBIAN* D;

    bool A_exists;
    bool B_exists;
    bool C_exists;
    bool D_exists;

    Exchange_COO_CSC A_coo;
    Exchange_COO_CSC B_coo;
    Exchange_COO_CSC C_coo;
    Exchange_COO_CSC D_coo;

    ExchangeJacobians(DATA* data, threadData_t* threadData, InfoGDOP& info);

    void init_jac(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg, FixedVector<FunctionMR>& mr);
    void init_jac_lfg(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionLFG>& lfg);
    void init_jac_mr(DATA* data, threadData_t* threadData, InfoGDOP& info, FixedVector<FunctionMR>& mr);

};

#endif // OPT_OM_EVALUATIONS
