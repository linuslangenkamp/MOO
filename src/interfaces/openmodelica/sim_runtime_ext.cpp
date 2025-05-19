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
 * only from those function rows r where both variables x_i and x_j appear (i.e., where ∂f/∂x_i and ∂f/∂x_j are nonzero).
 * The directional second derivative is accumulated using:
 *
 *     H[i,j] += ∑_r λ[r] · (J(x + h·s₁)[r] - J(x)[r]) / h
 *
 * where λ ∈ ℝᵐ is the adjoint vector. This corresponds to evaluating the contraction λᵗ · ∇²F(x) · v without forming the full Hessian.
 *
 * The structure enables efficient scheduling of these evaluations via pre-colored seed vectors, reducing the number of required sweeps.
 
 * The result is a `HESSIAN_PATTERN` C struct that includes:
 * - COO row/col index arrays for Hessian nonzeros (lower triangle).
 * - A lookup from (color₁, color₂) to `HessianEntry`, listing variable pairs and contributing rows.
 *
 * @param jac [in]  Pointer to a `JACOBIAN` struct (sparsity pattern and coloring).
 * @return    [out] Pointer to newly allocated `HESSIAN_PATTERN` struct, or NULL on error.
 *
 * This is a pure C data structure. Use `extern "C"` guards to include in C++ code.
 */
HESSIAN_PATTERN* generateHessianPattern(JACOBIAN* jac) {
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

  // 2. build M[v1, v2] = list of function rows where both appear
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
  hes_pattern->entries = (HessianEntry**)calloc(numColors * numColors, sizeof(HessianEntry*));
  hes_pattern->row = (int*)malloc(lnnz * sizeof(int));
  hes_pattern->col = (int*)malloc(lnnz * sizeof(int));
  hes_pattern->numColors = numColors;
  hes_pattern->lnnz = lnnz;
  hes_pattern->size = n_vars;
  hes_pattern->jac = jac;

  // 6. fill the coordinate format sparsity
  for (const auto& coo : cooMap) {
    int var_row = coo.first.first;
    int var_col = coo.first.second;
    int nz = coo.second;

    hes_pattern->row[nz] = var_row;
    hes_pattern->col[nz] = var_col;
  }

  // 7. fill HESSIAN_PATTERN.entries[c1][c2] -> HessianEntry
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
      if (entryCount == 0) continue;

      HessianEntry* entry = (HessianEntry*)malloc(sizeof(HessianEntry));
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

      hes_pattern->entries[entryIndexFromColors(c1, c2, numColors)] = entry;
    }
  }

  return hes_pattern;
}

void printHessianPattern(const HESSIAN_PATTERN* hes_pattern) {
  if (!hes_pattern) {
    printf("Hessian pattern is NULL.\n");
    return;
  }

  printf("\nAUGMENTED HESSIAN PATTERN\n");
  printf("Matrix size: %d x %d\n", hes_pattern->size, hes_pattern->size);
  printf("Number of colors: %d\n", hes_pattern->numColors);
  printf("Lower triangle NNZs: %d\n", hes_pattern->lnnz);
  printf("\nCoordinate Format (COO):\n");
  printf(" lnnz | Row | Col\n");
  printf("------------------\n");
  for (int i = 0; i < hes_pattern->lnnz; ++i) {
    printf(" %4d | %3d | %3d\n", i, hes_pattern->row[i], hes_pattern->col[i]);
  }

  printf("\nColor Pair Entries:\n");
  for (int c1 = 0; c1 < hes_pattern->numColors; c1++) {
    for (int c2 = 0; c2 <= c1; c2++) {  // symmetric lower triangle
      int idx = entryIndexFromColors(c1, c2, hes_pattern->numColors);
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
      printf("------------------------------------------------------------\n");  // separator between color pairs
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

void freeHessianPattern(HESSIAN_PATTERN* hes_pattern) {
  if (!hes_pattern) return;

  int numColorPairs = hes_pattern->numColors * hes_pattern->numColors;
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
  free(hes_pattern);
}
