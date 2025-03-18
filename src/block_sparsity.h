#include "util.h"
#include "problem.h"

enum class BlockType {
    Exact,
    Offset,
    RowOffset,
};

struct BlockSparsity {
    BlockType type;
    std::unique_ptr<std::unique_ptr<int[]>[]> block;
    int nnz;

    // block F only
    std::unique_ptr<int[]> row_offset_prev;
    std::unique_ptr<int[]> row_size;

    // block B only
    int off_prev;

    /* creates a dense lower triangular (with diagonal) block structure
       0 | 1 2 . s
       -----------
       1 | x      
       2 | x x    
       . | x x x  
       s | x x x x
    */
    static BlockSparsity createLowerTriangular(const int size, const BlockType block_type) {
        BlockSparsity b;
        b.block = std::make_unique<std::unique_ptr<int[]>[]>(size);
        for (int i = 0; i < size; i++) {
            b.block[i] = std::make_unique<int[]>(i+1);
        }
        switch (block_type) {
            case BlockType::Offset:
                b.off_prev = 0;
                break;
            case BlockType::RowOffset:
                b.row_offset_prev = std::make_unique<int[]>(size);
                b.row_size = std::make_unique<int[]>(size);
                break;
            default:
                break;
        }
        return b;
    }

    /* creates a dense rectangular block structure Rows x Cols
       0 | 1 2 . C
       -----------
       1 | x x x x
       . | x x x x
       R | x x x x
    */
    static BlockSparsity createRectangular(const int rows, const int cols, const BlockType block_type) {
        BlockSparsity b;
        b.block = std::make_unique<std::unique_ptr<int[]>[]>(rows);
        for (int i = 0; i < rows; i++) {
            b.block[i] = std::make_unique<int[]>(cols);
        }
        switch (block_type) {
            case BlockType::Offset:
                b.off_prev = 0;
                break;
            case BlockType::RowOffset:
                b.row_offset_prev = std::make_unique<int[]>(rows);
                b.row_size = std::make_unique<int[]>(rows);
                break;
            default:
                break;
        }
        return b;
    }

    /* creates a dense square block structure Size x Size
       0 | 1 . C
       -----------
       1 | x x x
       . | x x x
       R | x x x
    */
    inline static BlockSparsity createSquare(const int size, const BlockType block_type) {
        return createRectangular(size, size, block_type);
    }

    inline void insert(const int row, const int col, const int index) {
        block[row][col] = index;
        nnz++;
    }

    inline int get(const int row, const int col) {
        switch (type) {
            case BlockType::Exact:
                return block[row][col];
            default:
               throw std::runtime_error("Unknown BlockType in BlockSparsity::get().");
        }
    }


    inline int get(const int row, const int col, const int block_count) {
        switch (type) {
            // for A -> B: row, col, offset_prev - #nnz at end of blocktype (e.g. |A|), block_count (e.g. (i,j) = (0, 2) => 2)
            case BlockType::Offset:
                return off_prev + nnz * block_count + block[row][col];

            // for E -> F
            case BlockType::RowOffset:
                return row_offset_prev[row] + row_size[row] * block_count + block[row][col];
                
            default:
               return get(row, col);
        }
    }
};

struct OrderedIndexSet {
    struct Compare {
        bool operator()(const std::pair<int, int>& a, const std::pair<int, int>& b) const {
            if (a.first != b.first) {
                return a.first < b.first;
            } else {
                return a.second < b.second;
            }
        }
    };

    std::set<std::pair<int, int>, Compare> set;

    void insertSparsity(std::vector<HessianSparsity>& hes, int row_off, int col_off) {
        for (auto& coo : hes) {
            set.insert({coo.index1 + row_off, coo.index2 + col_off});
        }
    }

    inline size_t size() const {
        return set.size();
    }
};