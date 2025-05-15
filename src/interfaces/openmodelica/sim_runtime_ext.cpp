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
/*
void __evalHessian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac) {
  for (int color = 0; color < num_colors; ++color) {
    // Find all columns in this color class
    int* seed_cols = get_columns_with_color(color);
    int num_seeds = count(seed_cols);

    // For each seed direction e_i (i in seed_cols)
    for (int s = 0; s < num_seeds; ++s) {
      int i = seed_cols[s];

      // x + h e_i and x - h e_i
      perturb_x(x_plus, x, i, +h);
      perturb_x(x_minus, x, i, -h);

      // Evaluate Jacobians
      __evalJacobian(data, threadData, jacobian, NULL, J_plus);
      __evalJacobian(data, threadData, jacobian, NULL, J_minus);

      // For each col j with J(r, j) != 0 (can restrict using Jacobian sparsity pattern)
      for (int j = 0; j < n; ++j) {
        if (is_nonzero_jacobian_entry(i, j)) {
          for (int r = 0; r < m; ++r) {
            double d2f = (J_plus[r * n + j] - J_minus[r * n + j]) / (2 * h);
            store_hessian_entry(r, i, j, d2f); // Symmetric storage
          }
        }
      }
    }
  }
}
*/
