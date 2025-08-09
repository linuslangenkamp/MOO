#include <iomanip>
#include <cmath>

#include <simulation/radau/wrapper.h>
#include <base/log.h>

// ODE function fcn: dy/dx = f(x,y)
void fcn(int* n, double* x, double* y, double* dydx) {
    dydx[0] = -y[0];
}

void jac(int* n, double* x, double* y, int* ml, int* mu, double* pd, double* pdata) {
    pd[0] = -1.0; // df/dy = -1
}

void mas(int* n, int* m, double* data) {
    // no mass matrix for this problem
}

void solout(int* nr, double* x_old, double* x, double* y, double* rpar) {
    // CONTRA(I,S,CONT,LRC) 
    FixedTableFormat<4> table_format = {{12, 22, 22, 12}, {Align::Center, Align::Center, Align::Center, Align::Center}};

    double expected = std::exp(-(*x));
    double diff = expected - y[0];

    LOG_ROW(table_format,
        fmt::format("{:+.8f}", *x),
        fmt::format("{:+.15e}", y[0]),
        fmt::format("{:+.15e}", expected),
        fmt::format("{:+.8e}", diff)
    );
}

int radau_wrapper_test() {
    FixedTableFormat<4> table_format = {{12, 22, 22, 16}, {Align::Center, Align::Center, Align::Center, Align::Center}};

    LOG_START_MODULE(table_format, "Test: RADAU interface");
    LOG_HEADER(table_format, "Time", "y(t)", "Expected y(t)", "Error");
    LOG_DASHES(table_format);

    int n = 1;            // system size
    double x = 0.0;       // initial x
    double y[1] = {1.0};  // initial condition y(0) = 1
    double xend = 1.0;    // final x
    double h = 0.1;      // initial step size
    double rtol = 1e-6;   // relative tolerance
    double atol = 1e-6;   // absolute tolerance
    int itol = 0;         // tolerance type (0 = scalar)
    int ijac = 1;         // jacobian supplied (1 = yes)
    int mljac = 0;        // lower band width (not used here)
    int mujac = 0;        // upper band width (not used here)
    int imas = 0;         // mass matrix flag
    int mlmas = 0;
    int mumas = 0;
    int iout = 1;         // output function calls
    int lwork = 100;      // size of work array
    int liwork = 100;     // size of iwork array
    int idid = 0;         // status flag (output)
    double work[100] = {0};
    int iwork[100] = {0};
    double rpar[10] = {0};
    int ipar[10] = {0};

    radau_solver(
        &n, fcn, &x, y, &xend, &h, &rtol, &atol, &itol,
        jac, &ijac, &mljac, &mujac, mas, &imas, &mlmas, &mumas,
        solout, &iout, work, &lwork, iwork, &liwork, rpar, ipar, &idid
    );

    if (idid < 0) {
        LOG_ERROR("RADAU failed with error code: {}", idid);
        return idid;
    }

    LOG_DASHES(table_format);
    LOG_SUCCESS("RADAU finished successfully with return code: {}", idid);
    LOG("");

    return 0;
}
