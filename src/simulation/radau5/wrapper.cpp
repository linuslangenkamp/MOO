#include "wrapper.h"

// TODO: make this more simple: just provide tolerances time ControlTrajectory etc. => Obtain simulated Trajectory (optionally at predefined points)

void radau5_solver(
    int* n,
    void (*fcn)(int*, double*, double*, double*),
    double* x,
    double* y,
    double* xend,
    double* h,
    double* rtol,
    double* atol,
    int* itol,
    void (*jac)(int*, double*, double*, int*, int*, double*, double*),
    int* ijac,
    int* mljac,
    int* mujac,
    void (*mas)(int*, int*, double*),
    int* imas,
    int* mlmas,
    int* mumas,
    void (*solout)(int*, double*, double*, double*, double*),
    int* iout,
    double* work,
    int* lwork,
    int* iwork,
    int* liwork,
    double* rpar,
    int* ipar,
    int* idid)
{
    radau5_(
        n, fcn, x, y, xend, h, rtol, atol, itol,
        jac, ijac, mljac, mumas, mas, imas, mlmas, mumas,
        solout, iout, work, lwork, iwork, liwork, rpar, ipar, idid
    );
}
