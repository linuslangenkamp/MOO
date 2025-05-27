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

#include "debug_om.h"

typedef struct {
  int i;  // first variable index
  int j;  // second variable index
} VarPair;

/* Maps a (color1, color2) pair to all (i, j) variable pairs sharing these colors.
 * For each (i, j), stores the list of function rows f where both ∂f/∂xi and ∂f/∂xj are nonzero (overestimate).
 * Also stores the flat COO nz index for (i, j). */
typedef struct {
  // actual variable pairs, i.e. varPair[k] == (v1, v2); can also be accessed via HESSIAN->(row, col)[lnnzIndices[k]]
  VarPair* varPairs;        // is the variable pair contributing to the functions contributingRows[k]
  int** contributingRows;   // contributingRows[k] = functions affecting the varPair[k]
  int* numContributingRows; // number of rows for each pair
  int* lnnzIndices;         // mapping from variable pair to Hessian COO index
  int size;                 // number of variable pairs in this color group
} ColorPair;

/* Holds the compressed Hessian structure derived from a Jacobian.
 * COO format row/col lists lower-triangular nonzeros (∂²G/∂xi∂xj).
 * Variable pairs are grouped by (color1, color2) inside ColorPair blocks. */
typedef struct {
  /* this is an array of ptrs to ColorPair, is NULL if (c1, c2) is not contained */
  ColorPair** colorPairs;        // __getColorPairIndex(c1, c2) with c1 >= c2 -> variable pairs for color pair
  int* row;                      // flat COO row indices (i)
  int* col;                      // flat COO column indices (j)
  int size;                      // number of variables (Hessian is size × size)
  int numFuncs;                  // number of functions in the augmented Hessian
  int lnnz;                      // number of lower triangular nonzeros
  int** colsForColor;            // colsForColor[c] is an array of column indices in color c
  int* colorSizes;               // colorSizes[c] is the number of columns in colorCols[c]
  int numColors;                 // number of seed vector colors
  JACOBIAN* jac;                 // input Jacobian with sparsity + coloring
  int** cscJacIndexFromRowColor; // mapping of J[function / row][color] -> index in flat Jacobian CSC buffer
  modelica_real* ws_oldX;        // workspace array to remember old x values and seed vector for JVPs
  modelica_real** ws_baseJac;    // workspace stores all rows x colors of the base Jacobian J(x)
} HESSIAN_PATTERN;

/* always use this if accessing HESSIAN_PATTERN.colorPairs
 * returns the index of a colorPair (c1, c2) in the HESSIAN_PATTERN.colorPairs */
static inline int __getColorPairIndex(int c1, int c2) {
  if (c1 >= c2) return c1 * (c1 + 1) / 2 + c2;
  else return c2 * (c2 + 1) / 2 + c1;
}

static inline void __setSeedVector(int size, const int* cols, modelica_real value, modelica_real* seeds) {
  for (int i = 0; i < size; i++) { seeds[cols[i]] = value; }
}

HESSIAN_PATTERN* __generateHessianPattern(JACOBIAN* jac);

void __printHessianPattern(const HESSIAN_PATTERN* hes_pattern);

void __freeHessianPattern(HESSIAN_PATTERN* hes_pattern);

/* simple extension to evalJacobian of SimulationRuntime */
void __evalJacobian(DATA* data, threadData_t* threadData, JACOBIAN* jacobian, JACOBIAN* parentJacobian, modelica_real* jac);

void __evalHessianForwardDifferences(DATA* data, threadData_t* threadData, HESSIAN_PATTERN* hes_pattern, modelica_real h,
                                     modelica_real* lambda, modelica_real* jac_csc, modelica_real* hes);

void __forwardDiffHessianWrapper(void* args, modelica_real h, modelica_real* result);

// ===== EXTRAPOLATION =====

/* generic computation function of the form "result := f(args, h0)" */
typedef void (*Computation_fn_ptr)(void* args, modelica_real h0, modelica_real* result);

/**
 * @brief Workspace for Richardson extrapolation.
 * Stores intermediate results and metadata for extrapolation steps.
 */
typedef struct {
  modelica_real** ws_results;
  int resultSize;
  int maxSteps;
} ExtrapolationData;

/* augmented Hessian structure for richardson extrapolation scheme */
typedef struct {
  DATA* data;
  threadData_t* threadData;
  HESSIAN_PATTERN* hes_pattern;
  modelica_real* lambda;
  modelica_real* jac_csc;
} HessianFiniteDiffArgs;

ExtrapolationData* __initExtrapolationData(int resultSize, int maxSteps);

void __freeExtrapolationData(ExtrapolationData* extrData);

void __richardsonExtrapolation(ExtrapolationData* extrData, Computation_fn_ptr fn, void* args, modelica_real h0,
                              int steps, modelica_real stepDivisor, int methodOrder, modelica_real* result);

#endif // OPT_OM_EXTENSIONS_H
