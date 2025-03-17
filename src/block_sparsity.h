#include "util.h"

struct BlockSparsity {
    std::unique_ptr<std::unique_ptr<int[]>[]> block;
    int nnz;

    /* creates a dense diagonal, lower, diagonal block structure (DLD)
       0 | 1 2    0 | 1 2
       -------    -------
       1 | x 0    1 | D 0
       2 | x x    2 | L D

       size_1 = size of D_11, size_2 = size of D_22
    */
    static BlockSparsity createDLD(const int size_1, const int size_2) {
        BlockSparsity b;
        int size_total = size_1 + size_2;
        b.block = std::make_unique<std::unique_ptr<int[]>[]>(size_total);
        for (int i = 0; i < size_1; i++) {
            b.block[i] = std::make_unique<int[]>(size_1);
        }
        for (int i = size_1; i < size_total; i++) {
            b.block[i] = std::make_unique<int[]>(size_total);
        }
        return b;
    }

    /* creates a dense rectangular block structure Rows x Cols
       0 | 2
       -----
       1 | x
    */
    static BlockSparsity createRectangular(const int rows, const int cols) {
        BlockSparsity b;
        b.block = std::make_unique<std::unique_ptr<int[]>[]>(rows);
        for (int i = 0; i < rows; i++) {
            b.block[i] = std::make_unique<int[]>(cols);
        }
        return b;
    }

    /* creates a dense square block structure Size x Size
       0 | 1
       -----
       1 | x
    */
    inline static BlockSparsity createSquare(const int size) {
        return createRectangular(size, size);
    }

    void insert(const int i, const int j, const int index) {
        block[i][j] = index;
        nnz++;
    }
};