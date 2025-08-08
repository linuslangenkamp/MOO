#ifndef RADAU5_WRAPPER_H
#define RADAU5_WRAPPER_H

extern "C" {
    void radau5_(
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
        int* idid
    );

    void radau_(
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
        int* idid
    );
}

void radau_solver(
    int* n,           // N: size of the problem
    void (*fcn)(int*, double*, double*, double*), // FCN: function for dy/dx
    double* x,        // X: initial value of the independent variable
    double* y,        // Y: solution array
    double* xend,     // XEND: final value of the independent variable
    double* h,        // H: initial step size
    double* rtol,     // RTOL: relative tolerance
    double* atol,     // ATOL: absolute tolerance
    int* itol,        // ITOL: tolerance type
    void (*jac)(int*, double*, double*, int*, int*, double*, double*), // JAC: function for the Jacobian
    int* ijac,        // IJAC: jacobian type
    int* mljac,       // MLJAC: lower bandwidth of the Jacobian
    int* mujac,       // MUJAC: upper bandwidth of the Jacobian
    void (*mas)(int*, int*, double*), // MAS: function for the mass matrix
    int* imas,        // IMAS: mass matrix type
    int* mlmas,       // MLMAS: lower bandwidth of the mass matrix
    int* mumas,       // MUMAS: upper bandwidth of the mass matrix
    void (*solout)(int*, double*, double*, double*, double*), // SOLOUT: output function
    int* iout,        // IOUT: output type
    double* work,     // WORK: double precision work array
    int* lwork,       // LWORK: size of WORK
    int* iwork,       // IWORK: integer work array
    int* liwork,      // LIWORK: size of IWORK
    double* rpar,     // RPAR: user-defined double precision parameters
    int* ipar,        // IPAR: user-defined integer parameters
    int* idid         // IDID: success indicator
);

#endif // RADAU5_WRAPPER_H
