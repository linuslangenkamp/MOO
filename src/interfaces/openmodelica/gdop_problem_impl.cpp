#include "gdop_problem_impl.h"

// Dummy implementation of FullSweep_OM
FullSweep_OM::FullSweep_OM(FixedVector<FunctionLFG>&& lfg, Mesh& mesh, FixedVector<Bounds>&& g_bounds, 
                           bool has_lagrange, int f_size, int g_size, int x_size, int u_size, int p_size)
    : FullSweep(std::move(lfg), mesh, std::move(g_bounds), has_lagrange, f_size, g_size, x_size, u_size, p_size) {
    // Dummy constructor (does nothing)
}
void FullSweep_OM::callback_eval(const f64* xu_nlp, const f64* p) {
    // Dummy evaluation (does nothing)
}

void FullSweep_OM::callback_jac(const f64* xu_nlp, const f64* p) {
    // Dummy Jacobian (does nothing)
}

void FullSweep_OM::callback_hes(const f64* xu_nlp, const f64* p) {
    // Dummy Hessian (does nothing)
}

// Dummy implementation of BoundarySweep_OM
BoundarySweep_OM::BoundarySweep_OM(FixedVector<FunctionMR>&& mr, Mesh& mesh, FixedVector<Bounds>&& r_bounds,
                                   bool has_mayer, int r_size, int x_size, int p_size)
    : BoundarySweep(std::move(mr), mesh, std::move(r_bounds), has_mayer, r_size, x_size, p_size) {
    // Dummy constructor (does nothing)
}

void BoundarySweep_OM::callback_eval(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    // Dummy evaluation (does nothing)
}

void BoundarySweep_OM::callback_jac(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    // Dummy Jacobian (does nothing)
}

void BoundarySweep_OM::callback_hes(const f64* x0_nlp, const f64* xf_nlp, const f64* p) {
    // Dummy Hessian (does nothing)
}

Problem* create_gdop_om(DATA* data, Mesh& mesh) {
    // sizes
    int size_x = data->modelData->nStates;
    int size_u = data->modelData->nInputVars;
    int size_p = 0; // TODO: add this feature

    // TODO: figure this out we get some derivative? ptrs i guess?
    short der_index_m = -1;
    short der_indices_l[2] = {-1, -1};
    f64* address_m;
    f64* address_l;
    // this is really ugly IMO, fix this when ready for master!
    bool mayer = (data->callback->mayer(data, &address_m, &der_index_m) >= 0);
    bool lagrange = (data->callback->lagrange(data, &address_l, &der_indices_l[0], &der_indices_l[1]) >= 0);

    int size_g = data->modelData->nOptimizeConstraints;
    int size_r = data->modelData->nOptimizeFinalConstraints; // TODO: add *generic boundary* constraints later also at t=t0

    /* create functions and bounds */
    FixedVector<FunctionMR> mr((int)mayer + size_r);
    FixedVector<FunctionLFG> lfg((int)lagrange + size_x + size_g); // size_f == size_x

    /* create sparsity patterns */


    /* variable bounds */
    FixedVector<Bounds> x_bounds(size_x);
    FixedVector<Bounds> u_bounds(size_u);
    FixedVector<Bounds> p_bounds(size_p);

    for (int x = 0; x < size_x; x++) {
        x_bounds[x].lb = data->modelData->realVarsData[x].attribute.min;
        x_bounds[x].ub = data->modelData->realVarsData[x].attribute.max;
    }

    /* new generated function getInputVarIndices, just fills a index list */
    auto u_data_indices = std::make_unique<int>(size_u);
    data->callback->getInputVarIndices(data, u_data_indices.get());
    for (int u = 0; u < size_u; u++) {
        int u_index =  *(u_data_indices.get() + u);
        u_bounds[u].lb = data->modelData->realVarsData[u_index].attribute.min;
        u_bounds[u].ub = data->modelData->realVarsData[u_index].attribute.max;
    }

    /* constraint bounds */
    FixedVector<Bounds> g_bounds(size_g);
    FixedVector<Bounds> r_bounds(size_r);

    int first_index_g = data->modelData->nVariablesReal - (size_g + size_r);
    int first_index_r = data->modelData->nVariablesReal - size_r;
    printf("%d\n", first_index_g);
    for (int g = 0; g < size_g; g++) {
        g_bounds[g].lb = data->modelData->realVarsData[first_index_g + g].attribute.min;
        g_bounds[g].ub = data->modelData->realVarsData[first_index_g + g].attribute.max;
        
    }    
    printf("%f, %f\n", g_bounds[0].lb, g_bounds[0].ub);
    for (int r = 0; r < size_r; r++) {
        r_bounds[r].lb = data->modelData->realVarsData[first_index_r + r].attribute.min;
        r_bounds[r].ub = data->modelData->realVarsData[first_index_r + r].attribute.max;
    }

    /* for now we ignore xf fixed (need some steps in Backend to detect)
       and also ignore x0 non fixed, since too complicated
       => assume x(t_0) = x0 fixed, x(t_f) free to r constraint */
    FixedVector<std::optional<f64>> x0_fixed(size_x);
    FixedVector<std::optional<f64>> xf_fixed(size_x);

    
    
    /* set bounds and initial values */

    /*
    std::unique_ptr<FullSweep> fs(new FullSweep_OM(std::move(lfg), mesh, g_bounds));
    std::unique_ptr<BoundarySweep> bs(new BoundarySweep_OM(std::move(mr), mesh, r_bounds));
    std::make_shared<Problem>(std::move(fs), std::move(bs), 
                                     std::move(x_bounds), std::move(u_bounds), std::move(p_bounds),
                                     std::move(x0_fixed), std::move(xf_fixed))
    */

    /*
        for (size_t idx = 0; idx < 100; idx++)
        printf("%s\n", (data->modelData->realVarsData[idx].info.name));
    */

    return NULL;
}


