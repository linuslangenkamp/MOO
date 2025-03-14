#include "mesh.h"

/** 
 * @brief Create a mesh with fixed (h, p) for final time tf
 *
 * @param intervals  Number of Intervals
 * @param tf         Final Time
 * @param steps      Number of Nodes for each Interval
 * @return Mesh      Mesh
 */
Mesh Mesh::createEquidistantMeshFixedDegree(int intervals, double tf, int steps) {
    std::vector<double> grid(intervals + 1);
    std::vector<double> delta_t(intervals);
    std::vector<int> nodes(intervals);
    std::vector<std::vector<int>> acc_nodes(intervals, std::vector<int>(steps, 0));

    double h = tf / intervals;
    for (int i = 0; i < intervals; i++) {
        grid[i] = i * h;
    }
    grid[intervals - 1] = tf;

    for (int i = 0; i < intervals; i++) {
        delta_t[i] = h;
        nodes[i] = steps;
        for (int j = 0; j < steps; j++) {
            acc_nodes[i][j] = steps * i + j;
        }
    }
    int node_count = steps * intervals;
    return {intervals, tf, std::move(grid), std::move(delta_t), std::move(nodes), std::move(acc_nodes), node_count};
}

std::vector<std::vector<int>> Mesh::createAccOffsetXU(int off_x, int off_xu) {
    std::vector<std::vector<int>> off_acc_xu(intervals);
    int off = off_x;
    for (int i = 0; i < intervals; i++) {
        off_acc_xu[i].reserve(nodes[i]);
        for (int j = 0; j < nodes[i]; j++) {
            off_acc_xu[i].emplace_back(off);
            off += off_xu;
        }
    }
    return off_acc_xu;
}
