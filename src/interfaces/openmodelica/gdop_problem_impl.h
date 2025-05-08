#ifndef OPT_OM_GDOP_PROBLEM_IMPL
#define OPT_OM_GDOP_PROBLEM_IMPL

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

class FullSweep_OM : public FullSweep {
public:
    FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds, 
                           bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size);
    void callback_eval(const f64* xu_nlp, const f64* p) override;
    void callback_jac(const f64* xu_nlp, const f64* p) override;
    void callback_hes(const f64* xu_nlp, const f64* p) override;
};

class BoundarySweep_OM : public BoundarySweep {
public:
    BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                     bool has_mayer, int r_size, int x_size, int p_size);
    void callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) override;
    void callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) override;
    void callback_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p) override;
};

Problem* create_gdop_om(DATA* data, Mesh& mesh);

#endif // OPT_OM_GDOP_PROBLEM_IMPL