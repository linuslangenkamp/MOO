#ifndef OPT_MESH_H
#define OPT_MESH_H

#include "fixed_vector.h"
#include "collocation.h"
#include "util.h"


enum class InterpolationMethod {
    LINEAR = 0
};

struct Mesh {
    Mesh(int intervals, f64 tf, FixedVector<f64>&& grid, FixedVector<f64>&& delta_t, FixedField<f64, 2>&& t,
         FixedVector<int>&& nodes, FixedField<int, 2>&& acc_nodes, int node_count)
        : intervals(intervals), node_count(node_count), tf(tf), grid(std::move(grid)), delta_t(std::move(delta_t)),
          t(std::move(t)), nodes(std::move(nodes)), acc_nodes(std::move(acc_nodes)) {
    }

    int intervals;                // number of intervals
    int node_count;               // number of collocation nodes (sum over all intervals)
    f64 tf;                       // final time
    FixedVector<f64>   grid;      // grid base points
    FixedVector<f64>   delta_t;   // step size h for each interval
    FixedField<f64, 2> t;         // mesh points t_{i, j} = grid[i] + delta_t[i] * c[i][j] (c: collocation nodes for interval i)
    FixedVector<int>   nodes;     // number of collocation nodes p for each interval
    FixedField<int, 2> acc_nodes; // number of nodes to the left of index (i, j)


    static Mesh create_equidistant_fixed_stages(f64 tf, int intervals, int p, Collocation& collocation); 
    FixedField<int, 2> create_acc_offset_xu(int off_x, int off_xu);
    FixedField<int, 2> create_acc_offset_fg(int off_fg);
};

// given some data trajectories t, x(t), u(t), p -> interpolate w.r.t. mesh and collocation scheme -> new fitting guess
struct Trajectory {
    // x[i][j] = x_i at time t[j]
    std::vector<f64> t;
    std::vector<std::vector<f64>> x;
    std::vector<std::vector<f64>> u;
    std::vector<f64> p;
    InterpolationMethod interpolation;

    Trajectory() = default;

    // TODO: move me?
    Trajectory(std::vector<f64> t, std::vector<std::vector<f64>> x, std::vector<std::vector<f64>> u,
               std::vector<f64> p, InterpolationMethod interpolation = InterpolationMethod::LINEAR)
        : t(t), x(x), u(u), p(p), interpolation(interpolation) {
    }

    Trajectory interpolate(Mesh& mesh, Collocation& collocation);
    Trajectory linear_interpolation(Mesh& mesh, Collocation& collocation);

    void print();
};

#endif  // OPT_MESH_H
