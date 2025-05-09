#ifndef OPT_OM_GDOP_PROBLEM_IMPL
#define OPT_OM_GDOP_PROBLEM_IMPL

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "print_data_structures_om.h"

class FullSweep_OM : public FullSweep {
public:
    FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds, 
                           bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size);
    void callback_eval(const F64* xu_nlp, const F64* p) override;
    void callback_jac(const F64* xu_nlp, const F64* p) override;
    void callback_hes(const F64* xu_nlp, const F64* p) override;
};

class BoundarySweep_OM : public BoundarySweep {
public:
    BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                     bool has_mayer, int r_size, int x_size, int p_size);
    void callback_eval(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
    void callback_jac(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
    void callback_hes(const F64* x0_nlp, const F64* xf_nlp, const F64* p) override;
};

struct Jacobians_OM {
    JACOBIAN* A;
    JACOBIAN* B;
    JACOBIAN* C;
    JACOBIAN* D;

    bool A_exists;
    bool B_exists;
    bool C_exists;
    bool D_exists;

    Exchange_COO_CSC B_coo;
    Exchange_COO_CSC C_coo;
    Exchange_COO_CSC D_coo;

    std::unique_ptr<RowExchange_COO_CSC> C_mayer_coo;

    Jacobians_OM(DATA* data, threadData_t* threadData, bool mayer_exists, bool lagrange_exists, int x_size);

    /* for the GDOP, we get the following Jacobian structure (Note ^T means its a column, so xu = (x, u) columns and the functions are rows)
     * A = (f_xu)^T
     * B = (f_xu, L_xu, g_xu)^T
     * C = (f_xu, L_xu, M_xu, g_xu)^T
     * D = (r_xu)^T, but why not (M_xu, r_xu)?? so much easier!
    */
};

Problem* create_gdop_om(DATA* data, threadData_t* threadData, Mesh& mesh);

#endif // OPT_OM_GDOP_PROBLEM_IMPL
