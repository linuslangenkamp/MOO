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

void __extractJacobianCSCColumnDense(int col, JACOBIAN* jacobian, modelica_real* jac, modelica_real* denseCol) {
  SPARSE_PATTERN* sp = jacobian->sparsePattern;
  for (unsigned int idx = sp->leadindex[col]; idx < sp->leadindex[col + 1]; ++idx) {
    unsigned int row = sp->index[idx];
    denseCol[row] = jac[idx];
  }
}

void __resetJacobianCSCColumnDense(int col, JACOBIAN* jacobian, modelica_real* denseCol) {
  SPARSE_PATTERN* sp = jacobian->sparsePattern;
  for (unsigned int idx = sp->leadindex[col]; idx < sp->leadindex[col + 1]; ++idx) {
    unsigned int row = sp->index[idx];
    denseCol[row] = 0;
  }
}

/**
 * @brief Constructs a compressed Hessian sparsity pattern using Jacobian coloring.
 *
 * Given a sparse Jacobian J(x) of F: R^n → R^m, this builds the structure of the Hessian
 * of a scalar adjoint G(x) = Σ λ[i]·F[i](x), based on co-occurrence of variables in J(x).
 *
 * Variable pairs (i,j) are collected if they appear together in any function row f.
 * These define the nonzero Hessian structure (i.e. where ∂²G/∂xi∂xj != 0).
 *
 * A second coloring is induced: if variable x_i ∈ color c₁ and x_j ∈ color c₂, then (i,j) is assigned to the color pair (c₁, c₂).
 * This allows evaluating Hessian entries H[i,j] via directional finite differences: perturb x along seed vector s₁ (color c₁),
 * and apply a Jacobian-vector product with seed vector s₂ (color c₂).
 *
 * For each color pair (c₁, c₂), only a subset of Hessian entries is affected. Each entry H[i,j] receives contributions
 * only from those function rows r where both variables x_i and x_j appear (i.e. where ∂f/∂x_i and ∂f/∂x_j are nonzero).
 * The directional second derivative is given using:
 *
 *     H[i,j] = ∑_r λ[r] · ((J(x + h · s_1) - J(x)) · s_2 / h)[r]
 *
 * where λ ∈ ℝᵐ is the adjoint vector. This corresponds to evaluating the contraction λᵗ · ∇²F(x) · v without forming the full Hessian.
 *
 * The structure enables efficient scheduling of these evaluations via pre-colored seed vectors, reducing the number of required sweeps.
 
 * The result is a `HESSIAN_PATTERN` pure C struct that includes:
 * - COO row/col index arrays for Hessian nonzeros (lower triangle).
 * - A lookup from (color₁, color₂) to `HessianEntry`, listing variable pairs and contributing rows.
 *
 * @param jac [in]  Pointer to a `JACOBIAN` struct (sparsity pattern and coloring).
 * @return    [out] Pointer to newly allocated `HESSIAN_PATTERN` struct, or NULL on error.
 *
 */
HESSIAN_PATTERN* __generateHessianPattern(JACOBIAN* jac) {
  if (jac == nullptr or jac->sparsePattern == nullptr) { return nullptr; }

  int n_vars = jac->sizeCols;
  int n_funcs = jac->sizeRows;
  SPARSE_PATTERN* sp = jac->sparsePattern;
  int numColors = sp->maxColors;

  // 1. build adjacency list: which variables affect which functions
  std::vector<std::vector<int>> adj(n_funcs);
  for (int col = 0; col < n_vars; col++) {
    for (int nz = (int)sp->leadindex[col]; nz < (int)sp->leadindex[col + 1]; nz++) {
      int row = sp->index[nz];
      adj[row].push_back(col);
    }
  }

  // 2. build M[v1, v2] = list of function rows where both variables appear
  std::map<std::pair<int, int>, std::vector<int>> M;
  for (int f = 0; f < n_funcs; f++) {
    const auto& vars = adj[f];
    for (size_t i = 0; i < vars.size(); i++) {
      for (size_t j = 0; j <= i; j++) {
        int v1 = vars[i];
        int v2 = vars[j];
        if (v1 < v2) std::swap(v1, v2);
        M[{v1, v2}].push_back(f);
      }
    }
  }

  // 3. assign flat indices (lower nnz) directly from sorted M keys
  std::map<std::pair<int, int>, int> cooMap;
  int lnnz = 0;
  for (const auto& [pair, _] : M) {
    cooMap[pair] = lnnz++;
  }

  // 4. build color groups :: TODO: implement this in OpenModelica for the JACOBIAN
  std::vector<std::vector<int>> colorCols(numColors);
  for (int col = 0; col < n_vars; ++col) {
    int c = sp->colorCols[col];
    if (c > 0) {
      colorCols[c - 1].push_back(col);
    }
  }

  // 5. allocate pattern
  HESSIAN_PATTERN* hes_pattern = (HESSIAN_PATTERN*)malloc(sizeof(HESSIAN_PATTERN));
  hes_pattern->entries = (HessianEntry**)calloc(numColors * (numColors + 1) / 2, sizeof(HessianEntry*));
  hes_pattern->row = (int*)malloc(lnnz * sizeof(int));
  hes_pattern->col = (int*)malloc(lnnz * sizeof(int));
  hes_pattern->colsForColor = (int**)malloc(numColors * sizeof(int*));
  hes_pattern->colorSizes = (int*)malloc(numColors * sizeof(int));
  hes_pattern->numColors = numColors;
  hes_pattern->numFuncs = n_funcs;
  hes_pattern->size = n_vars;
  hes_pattern->lnnz = lnnz;
  hes_pattern->jac = jac;

  // 6. remember columns in each color
  for (int i = 0; i < numColors; ++i) {
    int size = colorCols[i].size();
    hes_pattern->colorSizes[i] = size;
    hes_pattern->colsForColor[i] = (int*)malloc(size * sizeof(int));
    memcpy(hes_pattern->colsForColor[i], colorCols[i].data(), size * sizeof(int));
  }

  // 7. fill the coordinate format sparsity
  for (const auto& coo : cooMap) {
    int var_row = coo.first.first;
    int var_col = coo.first.second;
    int nz = coo.second;

    hes_pattern->row[nz] = var_row;
    hes_pattern->col[nz] = var_col;
  }

  HessianEntry* entry;

  // 8. fill HESSIAN_PATTERN.entries[c1][c2] -> HessianEntry
  for (int c1 = 0; c1 < numColors; c1++) {
    for (int c2 = 0; c2 <= c1; c2++) {
      std::vector<std::vector<int>> rowsVec;
      std::vector<int> nnzIndicesVec;

      for (int i1 : colorCols[c1]) {
        for (int i2 : colorCols[c2]) {
          // copy and swap if needed
          int v1 = i1;
          int v2 = i2;
          if (v1 < v2){
            std::swap(v1, v2);
          }

          auto v_pair = std::make_pair(v1, v2);
          auto it = M.find(v_pair);
          if (it == M.end()) continue;

          auto cooIt = cooMap.find(v_pair);
          if (cooIt == cooMap.end()) continue;

          rowsVec.push_back(it->second);          // function rows
          nnzIndicesVec.push_back(cooIt->second); // flat Hessian index, nz index
        }
      }

      // create and allocate HessianEntry
      int entryCount = rowsVec.size();
      if (entryCount == 0) {
        entry = nullptr;
      }
      else {
        entry = (HessianEntry*)malloc(sizeof(HessianEntry));
        entry->size = entryCount;
        entry->rowIndices = (int**)malloc(entryCount * sizeof(int*));
        entry->rowSizes = (int*)malloc(entryCount * sizeof(int));
        entry->lnnzIndices = (int*)malloc(entryCount * sizeof(int));

        for (int i = 0; i < entryCount; i++) {
          int sz = rowsVec[i].size();
          entry->rowIndices[i] = (int*)malloc(sz * sizeof(int));
          memcpy(entry->rowIndices[i], rowsVec[i].data(), sz * sizeof(int));
          entry->rowSizes[i] = sz;
          entry->lnnzIndices[i] = nnzIndicesVec[i];
        }
      }

      hes_pattern->entries[__entryIndexFromColors(c1, c2)] = entry;
    }
  }

  return hes_pattern;
}

/**
 * @brief Compute Hessian-vector product λᵗH(x) using forward finite differences of the Jacobian.
 *
 * Approximates the entries of the Hessian matrix H(x) using first-order directional derivatives.
 * The method uses seed vector coloring for efficient evaluation and supports sparse Hessian structure.
 * Assumes the current point x has all controls and states set in `data->localData[0]->realVars`.
 *
 * @param[in]  data         Runtime simulation data structure.
 * @param[in]  threadData   Thread-local data (unused internally).
 * @param[in]  hes_pattern  Precomputed sparsity and coloring pattern for Hessian and Jacobian.
 * @param[in]  h            Perturbation step size (for now without nominal incorporation).
 * @param[in]  lambda       Adjoint vector (size = number of functions).
 * @param[out] hes          Output sparse Hessian values (COO format of hes_pattern, length = hes_pattern->nnz).
 */
void __evalHessianForwardDifferences(DATA* data, threadData_t* threadData, HESSIAN_PATTERN* hes_pattern, modelica_real h,
                                     modelica_real* lambda, modelica_real* hes) {
  /* 0. retrieve data */
  JACOBIAN* jacobian = hes_pattern->jac;
  modelica_real* seeds = jacobian->seedVars;
  modelica_real* jvp = jacobian->resultVars;

  /* TODO: Attention: for now we assume all inputs are control variables (to optimize); update this when needed! => iterate over all controls */
  int nStates = data->modelData->nStates;
  int uOffset = data->modelData->nVariablesReal - data->modelData->nInputVars - data->modelData->nStates;

  /* 1. evaluate all JVPs J(x) * s_{c} of the current point x */
  modelica_real** baseJacCols = (modelica_real**)malloc(hes_pattern->numColors * sizeof(modelica_real*));

  for (int c = 0; c < hes_pattern->numColors; c++) {
    baseJacCols[c] = (modelica_real*)calloc(hes_pattern->numFuncs, sizeof(modelica_real));
    __setSeedVector(hes_pattern->colorSizes[c], hes_pattern->colsForColor[c], 1, seeds);
    jacobian->evalColumn(data, threadData, jacobian, NULL);

    for (int colIndex = 0; colIndex < hes_pattern->colorSizes[c]; colIndex++) {
      int col = hes_pattern->colsForColor[c][colIndex];
      for (unsigned int nz = jacobian->sparsePattern->leadindex[col]; nz < jacobian->sparsePattern->leadindex[col + 1]; nz++) {
        int row = jacobian->sparsePattern->index[nz];
        baseJacCols[c][row] = jacobian->resultVars[row];
      }
    }

    __setSeedVector(hes_pattern->colorSizes[c], hes_pattern->colsForColor[c], 0, seeds);
  }

  /* allocate temp array to remember old x values and seed vector for JVPs */
  modelica_real* ws_old_x = (modelica_real*)malloc(data->modelData->nVariablesReal * sizeof(modelica_real));

  /* 2. loop over all colors c1 */
  for (int c1 = 0; c1 < hes_pattern->numColors; c1++) {
    /* 3. define seed vector s_{c_1} with all cols in c_1 active (implicitly) */
    /* 4. peturbate current x_{c_1} := x + h * s_{c_1} */
    for (int columnIndex = 0; columnIndex < hes_pattern->colorSizes[c1]; columnIndex++) {
      int col = hes_pattern->colsForColor[c1][columnIndex];
      int realVarsIndex = (col < nStates ? col : uOffset + col);

      /* remember the current realVars (to be perturbated) and perturbate */
      ws_old_x[columnIndex] = data->localData[0]->realVars[realVarsIndex];
      data->localData[0]->realVars[realVarsIndex] += h; /* TODO: incorporate nominals here for perturbation * nom */
    }

    /* 5. loop over all colors c2 with index less or equal to c_1 */
    for (int c2 = 0; c2 <= c1; c2++) {
      /* 6. define seed vector s_{c_2} with all cols in c_2 active */
      __setSeedVector(hes_pattern->colorSizes[c2], hes_pattern->colsForColor[c2], 1, seeds);

      /* 7. evaluate JVP J(x_{c_1}) * s_{c_2}: writes column to jvp = jacobian->resultVars */
      jacobian->evalColumn(data, threadData, jacobian, NULL);

      /* 8. retrieve augmented Hessian approximation */
      HessianEntry* colorPair = hes_pattern->entries[__entryIndexFromColors(c1, c2)];
      if (colorPair != nullptr) {
        for (int varPair = 0; varPair < colorPair->size; varPair++) {
          /* nz index in flattened Hessian array (COO format) */
          int nz = colorPair->lnnzIndices[varPair];

          /* rows (functions) where both ∂f/∂xi and ∂f/∂xj are nonzero */
          int* contributingRows = colorPair->rowIndices[varPair];
          int numberContributingRows = colorPair->rowSizes[varPair];

          /* second derivative eval at nz index */
          modelica_real der = 0;

          /* 10. Approximate directional second derivative:
           *     (1/h) ∑_{f ∈ rows} λ[f] · (J(x + h·s_{c₁})[s_{c₂}][f] - J(x)[s_{c₂}][f])
           *     where:
           *       - f / fnRow indexes function rows where both ∂f/∂xᵢ and ∂f/∂xⱼ are nonzero */
          for (int fIdx = 0; fIdx < numberContributingRows; fIdx++) {
            int fnRow = contributingRows[fIdx];
            der += lambda[fnRow] * (jvp[fnRow] - baseJacCols[c2][fnRow]);
          }

          /* store and divide by step size, TODO: make h depend on variable nominal
           * divide by nominal of hes_pattern->rows or cols variable at nz index */
          hes[nz] = der / h;
        }
      }

      /* 11. reset s_{c_2} */
      __setSeedVector(hes_pattern->colorSizes[c2], hes_pattern->colsForColor[c2], 0, seeds);
    }

    /* 12. reset perturbated x */ 
    for (int columnIndex = 0; columnIndex < hes_pattern->colorSizes[c1]; columnIndex++) {
      int col = hes_pattern->colsForColor[c1][columnIndex];
      int realVarsIndex = (col < nStates ? col : uOffset + col);
      data->localData[0]->realVars[realVarsIndex] = ws_old_x[columnIndex];
    }
  }

  for (int c = 0; c < hes_pattern->numColors; c++) {
    free(baseJacCols[c]);
  }
  free(baseJacCols);

  free(ws_old_x);
}

void __printHessianPattern(const HESSIAN_PATTERN* hes_pattern) {
  if (!hes_pattern) {
    printf("Hessian pattern is NULL.\n");
    return;
  }

  printf("\n=== AUGMENTED HESSIAN SPARISTY INFO ===\n");
  printf("Matrix size: %d x %d\n", hes_pattern->size, hes_pattern->size);
  printf("Lower triangle NNZ: %d\n", hes_pattern->lnnz);
  printf("Number of colors: %d\n", hes_pattern->numColors);
  printf("\nBase Jacobian Colors:\n");
  for (int c = 0; c < hes_pattern->numColors; ++c) {
      printf("    Color %d (size %d): ", c, hes_pattern->colorSizes[c]);
      for (int j = 0; j < hes_pattern->colorSizes[c]; ++j) {
          printf("%d ", hes_pattern->colsForColor[c][j]);
      }
      printf("\n");
  }
  printf("\nCoordinate Format (COO, lower triangle):\n");
  printf(" lnnz | Row | Col\n");
  printf("------------------\n");
  for (int i = 0; i < hes_pattern->lnnz; ++i) {
    printf(" %4d | %3d | %3d\n", i, hes_pattern->row[i], hes_pattern->col[i]);
  }

  printf("\nColor Pair Entries:\n");
  for (int c1 = 0; c1 < hes_pattern->numColors; c1++) {
    for (int c2 = 0; c2 <= c1; c2++) {  // symmetric lower triangle
      int idx = __entryIndexFromColors(c1, c2);
      HessianEntry* entry = hes_pattern->entries[idx];
      if (!entry) continue;

      printf("  Color pair (%d, %d): %d variable pairs\n", c1, c2, entry->size);
      for (int i = 0; i < entry->size; i++) {
        int nnzIdx = entry->lnnzIndices[i];
        printf("    VarPair: (%d, %d), nnz_index = %d, Functions = [",
               hes_pattern->row[nnzIdx], hes_pattern->col[nnzIdx], nnzIdx);

        for (int j = 0; j < entry->rowSizes[i]; j++) {
          printf("%d%s", entry->rowIndices[i][j], (j + 1 < entry->rowSizes[i]) ? ", " : "");
        }

        printf("]\n");
      }
      printf("----------------------------------------------------------------------------\n");  // separator between color pairs
    }
  }

  // sparsity
  int n = hes_pattern->size;
  printf("\n=== HESSIAN SPARSITY PLOT (λᵗ·∇²F) ===\n    ");
  for (int j = 0; j < n; ++j)
    printf("%d", j);
  printf("\n");

  // allocate and build symmetric sparsity map
  char** sparsity = (char**)calloc(n, sizeof(char*));
  for (int i = 0; i < n; ++i)
    sparsity[i] = (char*)calloc(n, sizeof(char));

  for (int k = 0; k < hes_pattern->lnnz; ++k) {
    int i = hes_pattern->row[k];
    int j = hes_pattern->col[k];
    sparsity[i][j] = 1;
    sparsity[j][i] = 1;  // symmetric for display
  }

  for (int i = 0; i < n; ++i) {
    printf("%2d: ", i);
    for (int j = 0; j < n; ++j)
      printf("%c", sparsity[i][j] ? '*' : ' ');
    printf("\n");
    free(sparsity[i]);
  }
  free(sparsity);
  printf("=====================================\n");
}

void __freeHessianPattern(HESSIAN_PATTERN* hes_pattern) {
  if (!hes_pattern) return;

  int numColorPairs = hes_pattern->numColors * (hes_pattern->numColors + 1) / 2;
  for (int i = 0; i < numColorPairs; i++) {
    HessianEntry* entry = hes_pattern->entries[i];
    if (!entry) continue;

    for (int j = 0; j < entry->size; j++) {
      free(entry->rowIndices[j]);
    }

    free(entry->rowIndices);
    free(entry->rowSizes);
    free(entry->lnnzIndices);
    free(entry);
  }

  free(hes_pattern->entries);
  free(hes_pattern->row);
  free(hes_pattern->col);

  if (hes_pattern->colsForColor) {
    for (int i = 0; i < hes_pattern->numColors; i++) {
      free(hes_pattern->colsForColor[i]);
    }
    free(hes_pattern->colsForColor);
  }
  free(hes_pattern->colorSizes);

  free(hes_pattern);
}


// ====== EXTRAPOLATION ======

int __richardsonExtrapolation(Computation_fn_ptr fn, void* args, modelica_real h0,
                              int steps, modelica_real stepDivisor, int methodOrder,
                              int resultSize, modelica_real* result) {
  /* simple wrapper fallback */
  if (steps <= 1) {
    fn(args, h0, result);
    return 0;
  }

  /* compute all stages for extrapolation */
  modelica_real** stageResults = (modelica_real**)malloc(steps * sizeof(modelica_real*));
  for (int i = 0; i < steps; i++) {
    stageResults[i] = (modelica_real*)malloc(resultSize * sizeof(modelica_real));
    modelica_real h = h0 / pow(stepDivisor, i);
    fn(args, h, stageResults[i]);
  }

  /* perform extrapolation: cancel taylor terms */
  for (int j = 0; j < resultSize; j++) {
    for (int i = 1; i < steps; i++) {
      modelica_real factor = pow(stepDivisor, methodOrder * i);
      stageResults[i][j] = (factor * stageResults[i][j] - stageResults[i - 1][j]) / (factor - 1);
    }
    result[j] = stageResults[steps - 1][j];
  }

  /* free stages */
  for (int i = 0; i < steps; i++) {
    free(stageResults[i]);
  }
  free(stageResults);

  return 0;
}

/* wrapper for __evalHessianForwardDifferences */
void __forwardDiffHessianWrapper(void* args, modelica_real h, modelica_real* result) {
  HessianFiniteDiffArgs* hessianArgs = (HessianFiniteDiffArgs*)args;
  __evalHessianForwardDifferences(hessianArgs->data, hessianArgs->threadData, hessianArgs->hes_pattern, h, hessianArgs->lambda, result);
}
