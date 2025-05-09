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
    F64* value;
};

struct HessianSparsity {
    int row;
    int col;
    F64* value;
};

// LFG - generic global function f(x, u, p, t)
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
    F64* eval;
    JacobianLFG jac;
    HessianLFG hes;
};

// MR - semi-generic boundary function r(x(t0), x(tf), p)
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
    F64* eval;
    JacobianMR jac;
    HessianMR hes;
};

/* exchange form from CSC -> COO and back */
struct Exchange_COO_CSC {
    FixedVector<int> row;
    FixedVector<int> col;

    // exchange mappings
    FixedVector<int> csc_to_coo;
    FixedVector<int> coo_to_csc;

    int nnz;

    Exchange_COO_CSC(int nnz) 
        : row(FixedVector<int>(nnz)),
          col(FixedVector<int>(nnz)),
          csc_to_coo(FixedVector<int>(nnz)),
          coo_to_csc(FixedVector<int>(nnz)),
          nnz(nnz)
    {}

    // static method for better readibility
    static Exchange_COO_CSC from_csc(const int* lead_col, const int* row_csc, int number_cols, int nnz) {
        return Exchange_COO_CSC(lead_col, row_csc, number_cols, nnz);
    }

private:
    /* simple sorting based CSC -> COO constructor with csc_to_coo and coo_to_csc */
    Exchange_COO_CSC(const int* lead_col, const int* row_csc, int number_cols, int nnz) :
                     row(nnz), col(nnz), csc_to_coo(nnz), coo_to_csc(nnz), nnz(nnz) {
        int nz = 0;
        for (int j = 0; j < number_cols; j++) {
            for (int i = lead_col[j]; i < lead_col[j + 1]; i++) {
                row[nz] = row_csc[i];
                col[nz] = j;
                csc_to_coo[nz] = nz;
                nz++;
            }
        }

        std::sort(csc_to_coo.begin(), csc_to_coo.end(), [&](int a, int b) {
            if (row[a] != row[b]) return row[a] < row[b];
            return col[a] < col[b];
        });

        FixedVector<int> sorted_row(nnz);
        FixedVector<int> sorted_col(nnz);
        for (int i = 0; i < nnz; i++) {
            sorted_row[i] = row[csc_to_coo[i]];
            sorted_col[i] = col[csc_to_coo[i]];
        }

        row = std::move(sorted_row);
        col = std::move(sorted_col);
        
        // permutation inverse
        for (int i = 0; i < nnz; ++i) {
            coo_to_csc[csc_to_coo[i]] = i;
        }
    }
};

/* one row of COO_CSC */
struct RowExchange_COO_CSC {
    FixedVector<int> col;

    // exchange mappings
    FixedVector<int> csc_to_coo;
    FixedVector<int> coo_to_csc;

    int nnz;

    // static method for better readibility
    static RowExchange_COO_CSC extract_row(Exchange_COO_CSC& exchange, int row_index) {
        return RowExchange_COO_CSC(exchange, row_index);
    }

private:
    RowExchange_COO_CSC(Exchange_COO_CSC& exchange, int row_index) {
        printf("int: %d", row_index);
        nnz = 0;
        int nz = 0;
        for (; nz < exchange.row.int_size(); nz++) {
            if (exchange.row[nz] == row_index) {
                nnz++;
            } 
            else if (exchange.row[nz] > row_index) {
                break;
            }
        }
        int start = nz - nnz;

        col        = (FixedVector<int>(nnz));
        csc_to_coo = (FixedVector<int>(nnz));
        coo_to_csc = (FixedVector<int>(nnz));

        for (int i = 0; i < nnz; i++) {
            col[i]        = exchange.col[start + i];
            csc_to_coo[i] = exchange.csc_to_coo[start + i];
            coo_to_csc[i] = exchange.coo_to_csc[start + i];
        }
    }
};

#endif // OPT_NLP_STRUCTS_H
