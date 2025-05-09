/* Linearization */
#include "include_test_model.h"
#if defined(__cplusplus)
extern "C" {
#endif
const char *include_test_linear_model_frame()
{
  return "model linearized_model \"include_test\"\n"
  "  parameter Integer n = 2 \"number of states\";\n"
  "  parameter Integer m = 1 \"number of inputs\";\n"
  "  parameter Integer p = 3 \"number of outputs\";\n"
  "\n"
  "  parameter Real x0[n] = %s;\n"
  "  parameter Real u0[m] = %s;\n"
  "\n"
  "  parameter Real A[n, n] =\n\t[%s];\n\n"
  "  parameter Real B[n, m] =\n\t[%s];\n\n"
  "  parameter Real C[p, n] =\n\t[%s];\n\n"
  "  parameter Real D[p, m] =\n\t[%s];\n\n"
  "\n"
  "  Real x[n](start=x0);\n"
  "  input Real u[m](start=u0);\n"
  "  output Real y[p];\n"
  "\n"
  "  Real 'x_x1' = x[1];\n"
  "  Real 'x_x2' = x[2];\n"
  "  Real 'u_u' = u[1];\n"
  "  Real 'y_$OMC$objectMayerTerm' = y[1];\n"
  "  Real 'y_$con$CONSTR' = y[2];\n"
  "  Real 'y_cost_m' = y[3];\n"
  "equation\n"
  "  der(x) = A * x + B * u;\n"
  "  y = C * x + D * u;\n"
  "end linearized_model;\n";
}
const char *include_test_linear_model_datarecovery_frame()
{
  return "model linearized_model \"include_test\"\n"
  "  parameter Integer n = 2 \"number of states\";\n"
  "  parameter Integer m = 1 \"number of inputs\";\n"
  "  parameter Integer p = 3 \"number of outputs\";\n"
  "  parameter Integer nz = 8 \"data recovery variables\";\n"
  "\n"
  "  parameter Real x0[n] = %s;\n"
  "  parameter Real u0[m] = %s;\n"
  "  parameter Real z0[nz] = %s;\n"
  "\n"
  "  parameter Real A[n, n] =\n\t[%s];\n\n"
  "  parameter Real B[n, m] =\n\t[%s];\n\n"
  "  parameter Real C[p, n] =\n\t[%s];\n\n"
  "  parameter Real D[p, m] =\n\t[%s];\n\n"
  "  parameter Real Cz[nz, n] =\n\t[%s];\n\n"
  "  parameter Real Dz[nz, m] =\n\t[%s];\n\n"
  "\n"
  "  Real x[n](start=x0);\n"
  "  input Real u[m](start=u0);\n"
  "  output Real y[p];\n"
  "  output Real z[nz];\n"
  "\n"
  "  Real 'x_x1' = x[1];\n"
  "  Real 'x_x2' = x[2];\n"
  "  Real 'u_u' = u[1];\n"
  "  Real 'y_$OMC$objectMayerTerm' = y[1];\n"
  "  Real 'y_$con$CONSTR' = y[2];\n"
  "  Real 'y_cost_m' = y[3];\n"
  "  Real 'z_$OMC$objectMayerTerm' = z[1];\n"
  "  Real 'z_cost_m' = z[2];\n"
  "  Real 'z_k1' = z[3];\n"
  "  Real 'z_k2' = z[4];\n"
  "  Real 'z_k3' = z[5];\n"
  "  Real 'z_k4' = z[6];\n"
  "  Real 'z_k5' = z[7];\n"
  "  Real 'z_u' = z[8];\n"
  "equation\n"
  "  der(x) = A * x + B * u;\n"
  "  y = C * x + D * u;\n"
  "  z = Cz * x + Dz * u;\n"
  "end linearized_model;\n";
}
#if defined(__cplusplus)
}
#endif

