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

void Mesh::update(std::unique_ptr<MeshUpdate> mesh_update, Collocation& collocation) {
    grid = std::move(mesh_update->new_grid);
    nodes = std::move(mesh_update->new_nodes_per_interval);
    intervals = nodes.size();

    // update delta_t
    FixedVector<f64> new_delta_t(intervals);
    for (int i = 0; i < intervals; i++) {
        new_delta_t[i] = grid[i + 1] - grid[i];
    }

    // update node_count
    node_count = std::accumulate(nodes.begin(), nodes.end(), 0);

    // update t and acc_nodes: allocate per row based on actual number of nodes
    FixedField<f64, 2> new_t(intervals);
    FixedField<int, 2> new_acc_nodes(intervals);

    int global_index = 0;
    for (int i = 0; i < intervals; i++) {
        int p = nodes[i];
        f64 h = new_delta_t[i];

        new_t[i] = FixedVector<f64>(p);
        new_acc_nodes[i] = FixedVector<int>(p);

        for (int j = 0; j < p; j++) {
            new_t[i][j] = grid[i] + h * collocation.c[p][j];
            new_acc_nodes[i][j] = global_index++;
        }
    }

    // move into structure
    delta_t = std::move(new_delta_t);
    t = std::move(new_t);
    acc_nodes = std::move(new_acc_nodes);
}
