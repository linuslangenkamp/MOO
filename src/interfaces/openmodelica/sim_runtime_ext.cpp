#include "sim_runtime_ext.h"

void new_evalJacobian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac, new_JACOBIAN_OUTPUT_FORMAT bufferFormat)
{
  int i,j,k,l,ii = 0;
  const SPARSE_PATTERN* sp = jacobian->sparsePattern;

  if (bufferFormat == new_JACOBIAN_OUTPUT_FORMAT::new_JAC_OUTPUT_DENSE) {
    memset(jac, 0.0, (jacobian->sizeRows) * (jacobian->sizeCols) * sizeof(modelica_real));
  } else if (bufferFormat == new_JACOBIAN_OUTPUT_FORMAT::new_JAC_OUTPUT_CSC) {
    memset(jac, 0.0, (jacobian->sparsePattern->numberOfNonZeros) * sizeof(modelica_real));
  }

  /* evaluate constant equations of Jacobian */
  if (jacobian->constantEqns != NULL) {
    jacobian->constantEqns(data, threadData, jacobian, parentJacobian);
  }

  /* evaluate Jacobian */
  for (i = 0; i < sp->maxColors; i++) {
    /* activate seed variable for the corresponding color */
    for (j = 0; j < jacobian->sizeCols; j++) /* TODO: maybe refactor colors as int** of dim (#colors, #size_color_j) of col indices? */
      if (sp->colorCols[j]-1 == i)
        jacobian->seedVars[j] = 1.0;

    /* evaluate Jacobian column */
    jacobian->evalColumn(data, threadData, jacobian, parentJacobian);

    for (j = 0; j < jacobian->sizeCols; j++) { 
      if (sp->colorCols[j]-1 == i) { /* TODO: maybe refactor colors as int** of dim (#colors, #size_color_j) of col indices? => just loop here over col indices */
        for (ii = sp->leadindex[j]; ii < sp->leadindex[j+1]; ii++) {
          l = sp->index[ii];
          if (bufferFormat == new_JACOBIAN_OUTPUT_FORMAT::new_JAC_OUTPUT_DENSE) {
            /* dense output buffer */
            
            k = j*jacobian->sizeRows + l;
            jac[k] = jacobian->resultVars[l];  /* solverData->xScaling[j]; */
          }
          else if (bufferFormat == new_JACOBIAN_OUTPUT_FORMAT::new_JAC_OUTPUT_CSC) {
            /* CSC sparse output buffer */
            //printf("CSC: %d, %d = %f\n", j, l, jacobian->resultVars[l]);
            jac[ii] = jacobian->resultVars[l]; /* solverData->xScaling[j]; */
          }
        }
        /* de-activate seed variable for the corresponding color */
        jacobian->seedVars[j] = 0.0;
      }
    }
  }
}
