#include "nlp.h"

void NLP::initSizesOffsets() {
    off_x = problem->xSize;
    off_u = problem->uSize;
    off_p = problem->pSize;
    off_xu = off_x + off_u;
    off_acc_xu = mesh->createAccOffsetXU(off_x, off_xu);
    off_xu_total = off_acc_xu.back().back() + off_xu;
    number_vars = off_xu_total + problem->pSize;
    number_constraints;
}

void NLP::init() {
    // init all the structures for optimization
    // n, m, nnz_jac, nnz_hess
    // (create scaling), create NLP bounds
    // sparsity pattern
    initSizesOffsets();


    // calculate n, m, nnz_jac, nnz_hes
    // create Jacobian and Hessian Sparsity patterns on the way

    // #constraints
    //m = sz(problem->A) + sz(problem->R) + (sz(problem->F) + sz(problem->G)) * rk.steps * mesh.intervals;


}
