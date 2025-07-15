#include "mesh.h"

/** 
 * @brief Create a mesh with fixed (h, p) for final time tf
 *
 * @param intervals  Number of Intervals
 * @param tf         Final Time
 * @param stages     Number of Nodes for each Interval
 * @return Mesh      Mesh
 */
Mesh Mesh::create_equidistant_fixed_stages(f64 tf, int intervals, int stages, Collocation& collocation) {
    FixedVector<f64> grid(intervals + 1);
    FixedVector<f64> delta_t(intervals);
    FixedVector<int> nodes(intervals);
    FixedField<int, 2> acc_nodes(intervals, stages);
    FixedField<f64, 2> t(intervals, stages);

    f64 h = tf / intervals;
    for (int i = 0; i < intervals; i++) {
        grid[i] = i * h;
    }
    grid[intervals] = tf;

    for (int i = 0; i < intervals; i++) {
        delta_t[i] = h;
        nodes[i] = stages;
        for (int j = 0; j < stages; j++) {
            acc_nodes[i][j] = stages * i + j;
            t[i][j] = grid[i] + delta_t[i] * collocation.c[stages][j];
        }
    }
    int node_count = stages * intervals;
    return {intervals, tf, std::move(grid), std::move(delta_t), std::move(t), std::move(nodes), std::move(acc_nodes), node_count};
}

FixedField<int, 2> Mesh::create_acc_offset_xu(int off_x, int off_xu) {
    FixedField<int, 2> off_acc_xu(intervals);
    int off = off_x;
    for (int i = 0; i < intervals; i++) {
        off_acc_xu[i] = FixedVector<int>(nodes[i]);
        for (int j = 0; j < nodes[i]; j++) {
            off_acc_xu[i][j] = off;
            off += off_xu;
        }
    }
    return off_acc_xu;
}

FixedField<int, 2> Mesh::create_acc_offset_fg(int off_fg) {
    FixedField<int, 2> acc_fg = acc_nodes;
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < nodes[i] ; j++) {
            acc_fg[i][j] *= off_fg;
        }
    }
    return acc_fg;
}
