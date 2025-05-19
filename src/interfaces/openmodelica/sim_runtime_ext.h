#ifndef OPT_OM_EXTENSIONS_H
#define OPT_OM_EXTENSIONS_H

/** 
 * has all the missing structures and functions that should be implemented in OpenModelica SimulationRuntime
 * we collect them for now
 */
#include <stdlib.h>
#include <unordered_map>
#include <map>
#include <utility>

#include "simulation_data.h"

#include <nlp/instances/gdop/problem.h>

#include "debug_om.h"

/* Maps a (color1, color2) pair to all (i,j) variable pairs sharing these colors.
 * For each (i,j), stores the list of function rows f where both ∂f/∂xi and ∂f/∂xj are nonzero (overestimate).
 * Also stores the flat COO index for (i,j). */
typedef struct {
  int** rowIndices;  // rowIndices[k][]: function rows for k-th variable pair
  int* rowSizes;     // number of rows for each pair
  int* lnnzIndices;  // mapping from variable pair to Hessian COO index
  int size;          // number of variable pairs in this color group
} HessianEntry;

/* Holds the compressed Hessian structure derived from a Jacobian.
 * COO format row/col lists lower-triangular nonzeros (∂²G/∂xi∂xj).
 * Variable pairs are grouped by (color1, color2) into HessianEntry blocks. */
typedef struct {
  HessianEntry** entries;  // [c1][c2] -> variable pairs for color pair
  int* col;                // COO column indices (j)
  int* row;                // COO row indices (i)
  int size;                // number of variables (Hessian is size × size)
  int lnnz;                // number of lower triangular nonzeros
  int numColors;           // number of seed vector colors
  JACOBIAN* jac;           // input Jacobian with sparsity + coloring
} HESSIAN_PATTERN;

static inline int entryIndexFromColors(int c1, int c2, int numColors) { return c1 * numColors + c2; }

HESSIAN_PATTERN* generateHessianPattern(JACOBIAN* jac);
void printHessianPattern(const HESSIAN_PATTERN* hes_pattern);
void freeHessianPattern(HESSIAN_PATTERN* hes_pattern);

/* simple extension to evalJacobian of SimulationRuntime */
void __evalJacobian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac);

/* numerical Hessian using foward differences on OpenModelica Jacobian */
void __evalForwardDifferencesHessian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian,
                                     modelica_real h, modelica_real* jac, modelica_real* hes);

/* numerical Hessian using extrapolation on foward differences and OpenModelica Jacobian */
void __evalNumericalHessianExtrapolation(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian,
                                         modelica_real h0, int steps, modelica_real* jac, modelica_real** hes);


#endif // OPT_OM_EXTENSIONS_H
