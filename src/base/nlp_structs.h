#ifndef OPT_NLP_STRUCTS_H
#define OPT_NLP_STRUCTS_H

#include <vector>
#include <algorithm>

#include "fixed_vector.h"
#include "util.h"

struct Bounds {
    F64 lb = MINUS_INFINITY;
    F64 ub = PLUS_INFINITY;
};

struct JacobianSparsity {
    int col;
    int buf_index;
};

struct HessianSparsity {
    int row;
    int col;
    int buf_index;
};

// LFG generic global function f(x, u, p, t)
// used for Lagrange term (L), dynamic (F), path (G) in GDOP

struct JacobianLFG {
    // coordinate format jacobian for LFGH functions
    std::vector<JacobianSparsity> dx;
    std::vector<JacobianSparsity> du;
    std::vector<JacobianSparsity> dp;

    inline int nnz() const {
        return dx.size() + du.size() + dp.size();
    }
};

struct HessianLFG {
    // coordinate format hessian for LFG functions
    std::vector<HessianSparsity> dx_dx;
    std::vector<HessianSparsity> du_dx;
    std::vector<HessianSparsity> du_du;
    std::vector<HessianSparsity> dp_dx;
    std::vector<HessianSparsity> dp_du;
    std::vector<HessianSparsity> dp_dp;

    inline int nnz() const {
        return dx_dx.size() + du_dx.size() + du_du.size() + dp_dx.size() + dp_du.size() + dp_dp.size();
    }
};

struct FunctionLFG {
    int buf_index;
    JacobianLFG jac;
    HessianLFG hes;
};

// MRsemi-generic boundary function r(x(t0), x(tf), p)
// used for Mayer term (M), boundary constraints (R) in GDOP

struct JacobianMR {
    // coordinate format jacobian for MR functions
    std::vector<JacobianSparsity> dx0;
    std::vector<JacobianSparsity> dxf;
    std::vector<JacobianSparsity> dp;

    inline int nnz() const {
        return dx0.size() + dxf.size() + dp.size();
    }
};

struct HessianMR {
    // coordinate format hessian for MR functions
    std::vector<HessianSparsity> dx0_dx0;
    std::vector<HessianSparsity> dxf_dx0;
    std::vector<HessianSparsity> dxf_dxf;
    std::vector<HessianSparsity> dp_dx0;
    std::vector<HessianSparsity> dp_dxf;
    std::vector<HessianSparsity> dp_dp;

    inline int nnz() const {
        return dx0_dx0.size() + dxf_dx0.size() + dxf_dxf.size() + dp_dx0.size() + dp_dxf.size() + dp_dp.size();
    }
};

struct FunctionMR {
    int buf_index;
    JacobianMR jac;
    HessianMR hes;
};

/* exchange form from CSC -> COO and back */
struct Exchange_COO_CSC {
public:
    FixedVector<int> row;
    FixedVector<int> col;

    FixedVector<int> __coo_to_csc;
    FixedVector<int> __csc_to_coo;

    int nnz = 0;
    int nnz_offset = 0;

    Exchange_COO_CSC() = default;

    Exchange_COO_CSC(int nnz) 
        : row(FixedVector<int>(nnz)),
          col(FixedVector<int>(nnz)),
          __coo_to_csc(FixedVector<int>(nnz)),
          __csc_to_coo(FixedVector<int>(nnz)),
          nnz(nnz),
          nnz_offset(0)
    {}

    /* local access (for standard CSC blocks) */
    inline int csc_to_coo(int index) {
        return __csc_to_coo[index];
    }

    inline int coo_to_csc(int local_index) {
        return __coo_to_csc[local_index];
    }

    /* embedded/global access (for offset blocks, e.g., [*, B](COO)) */
    inline int csc_to_coo_get_global(int index) {
        return __csc_to_coo[index] + nnz_offset;
    }

    inline int coo_to_csc_from_global(int global_index) {
        return __coo_to_csc[global_index - nnz_offset];
    }

    int row_nnz(int row_index) {
        int count = 0;
        for (int nz = 0; nz < row.int_size(); nz++) {
            if (row[nz] == row_index) {
                count++;
            }
            else if (row[nz] > row_index) {
                return count;
            }
        }
        return count;
    }

    /**
     * @brief Converts a matrix from CSC format to COO format with optional row reordering and global offset mapping.
     *
     * This struct performs a transformation from the CSC (Compressed Sparse Column) format
     * to the COO (Coordinate) format and stores permutation mappings between the two.
     *
     * Optionally, a specific row in the CSC format can be moved to index 0 before conversion.
     * This is useful for reordering matrices when a particular row must appear first in the
     * COO representation (e.g. Lagrange term in Lfg(COO) / fLg(CSC) sortings).
     *
     * @param lead_col Index of the first non-zero of each column (size = number_cols + 1).
     * @param row_csc  Row indices for non-zero elements in CSC format (size = nnz).
     * @param number_cols Number of columns in the matrix.
     * @param nnz Total number of non-zero elements in this CSC block.
     * @param move_to_first_row
     *       `-1` → No permutation, standard behavior (preserve row order).
     *       `r`  → Treat original row `r` as new row 0. All rows `< r` are shifted up by +1.
     * @param nnz_offset Offset of this CSC block in a global COO matrix.
     *        If a global COO is composed as [*, B], where this block is B, then nnz_offset = nnz(*)
     *        allows correct indexing through `*_embedded()` methods.
     *
     * @note The resulting COO matrix is sorted by row, then by column.
     *       The mappings `coo_to_csc` and `csc_to_coo` store **local (block-wise)** indices.
     *       Use the `*_global()` methods to map to/from global COO indices when embedded in a larger sparsity structure.
     */
    static Exchange_COO_CSC from_csc(const int* lead_col, const int* row_csc, int number_cols, int nnz, int move_to_first_row = -1, int nnz_offset = 0) {
        return Exchange_COO_CSC(lead_col, row_csc, number_cols, nnz, move_to_first_row, nnz_offset);
    }

private:
    Exchange_COO_CSC(const int* lead_col, const int* row_csc, int number_cols, int nnz, int move_to_first_row = -1, int nnz_offset = 0)
        : row(nnz), col(nnz), __coo_to_csc(nnz), __csc_to_coo(nnz), nnz(nnz), nnz_offset(nnz_offset) {
        int nz = 0;
        for (int curr_col = 0; curr_col < number_cols; ++curr_col) {
            for (int i = lead_col[curr_col]; i < lead_col[curr_col + 1]; i++) {
                int curr_row = row_csc[i];
                if (move_to_first_row >= 0) {
                    if      (curr_row == move_to_first_row) curr_row = 0;
                    else if (curr_row <  move_to_first_row) curr_row++;
                }
                row[nz]          = curr_row;
                col[nz]          = curr_col;
                __coo_to_csc[nz] = nz;
                nz++;
            }
        }

        std::sort(__coo_to_csc.begin(), __coo_to_csc.end(), [&](int a, int b) {
            return (row[a] != row[b]) ? (row[a] < row[b]) : (col[a] < col[b]);}
        );

        FixedVector<int> sorted_row(nnz);
        FixedVector<int> sorted_col(nnz);
        for (int i = 0; i < nnz; i++) {
            sorted_row[i] = row[__coo_to_csc[i]];
            sorted_col[i] = col[__coo_to_csc[i]];
            __csc_to_coo[__coo_to_csc[i]] = i;
        }

        row = std::move(sorted_row);
        col = std::move(sorted_col);
    }
};

#endif // OPT_NLP_STRUCTS_H
