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
    return Mesh{intervals, tf, std::move(grid), std::move(delta_t), std::move(t), std::move(nodes), std::move(acc_nodes), node_count};
}

Mesh::Mesh(std::unique_ptr<MeshUpdate> mesh_update, Collocation& collocation)
    : grid(std::move(mesh_update->new_grid)),
      nodes(std::move(mesh_update->new_nodes_per_interval))
{
    intervals = nodes.int_size();
    tf = grid.back();

    // compute delta_t
    delta_t = FixedVector<f64>(intervals);
    for (int i = 0; i < intervals; i++) {
        delta_t[i] = grid[i + 1] - grid[i];
    }

    // compute node_count
    node_count = std::accumulate(nodes.begin(), nodes.end(), 0);

    // allocate t and acc_nodes
    t = FixedField<f64, 2>(intervals);
    acc_nodes = FixedField<int, 2>(intervals);

    int global_index = 0;
    for (int i = 0; i < intervals; i++) {
        int p = nodes[i];
        f64 h = delta_t[i];

        t[i] = FixedVector<f64>(p);
        acc_nodes[i] = FixedVector<int>(p);

        for (int j = 0; j < p; j++) {
            t[i][j] = grid[i] + h * collocation.c[p][j];
            acc_nodes[i][j] = global_index++;
        }
    }
}

void Mesh::move_from(Mesh&& other) {
    intervals = other.intervals;
    tf = other.tf;
    grid = std::move(other.grid);
    delta_t = std::move(other.delta_t);
    t = std::move(other.t);
    nodes = std::move(other.nodes);
    acc_nodes = std::move(other.acc_nodes);
    node_count = other.node_count;
}

std::vector<f64> Mesh::get_flat_t(bool with_zero) const {
    std::vector<f64> flat_t;

    flat_t.reserve(node_count + (with_zero ? 1 : 0));

    if (with_zero) {
        flat_t.push_back(0.0);
    }

    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < nodes[i]; j++) {
            flat_t.push_back(t[i][j]);
        }
    }

    return flat_t;
}
