#include "sim_runtime_ext.h"

void __evalJacobian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac)
{
  size_t color, column, nz_csc;
  const SPARSE_PATTERN* sparse = jacobian->sparsePattern;

  memset(jac, 0.0, (jacobian->sparsePattern->numberOfNonZeros) * sizeof(modelica_real));

  /* evaluate constant equations of Jacobian */
  if (jacobian->constantEqns != NULL) {
    jacobian->constantEqns(data, threadData, jacobian, parentJacobian);
  }

  /* evaluate Jacobian */
  for (color = 0; color < sparse->maxColors; color++) {
    /* activate seed variable for the corresponding color */
    for (column = 0; column < jacobian->sizeCols; column++) /* TODO: maybe refactor colors as int** of dim (#colors, #size_color_j) of col indices? */
      if (sparse->colorCols[column] - 1 == color)
        jacobian->seedVars[column] = 1.0;

    /* evaluate Jacobian column */
    jacobian->evalColumn(data, threadData, jacobian, parentJacobian);

    for (column = 0; column < jacobian->sizeCols; column++) { 
      if (sparse->colorCols[column] - 1 == color) { /* TODO: maybe refactor colors as int** of dim (#colors, #size_color_j) of col indices? => just loop here over col indices */
        for (nz_csc = sparse->leadindex[column]; nz_csc < sparse->leadindex[column+1]; nz_csc++) {
          /* CSC sparse output buffer */
          jac[nz_csc] = jacobian->resultVars[sparse->index[nz_csc]]; /* solverData->xScaling[column]; */
        }
        /* de-activate seed variable for the corresponding color */
        jacobian->seedVars[column] = 0.0;
      }
    }
  }
}
