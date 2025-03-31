#ifndef OPT_MESH_H
#define OPT_MESH_H

#include "collocation.h"
#include "fixed_vector.h"
#include "util.h"


enum class InterpolationMethod {
    LINEAR = 0
};

struct Mesh {
    Mesh(int intervals, double tf, FixedVector<double>&& grid, FixedVector<double>&& delta_t, FixedField<double, 2>&& t,
         FixedVector<int>&& nodes, FixedField<int, 2>&& acc_nodes, int node_count)
        : intervals(intervals), node_count(node_count), tf(tf), grid(std::move(grid)), delta_t(std::move(delta_t)),
          t(std::move(t)), nodes(std::move(nodes)), acc_nodes(std::move(acc_nodes)) {
    }

    int intervals;                   // number of intervals
    int node_count;                  // number of collocation nodes (sum over all intervals)
    double tf;                       // final time
    FixedVector<double>   grid;      // grid base points
    FixedVector<double>   delta_t;   // step size h for each interval
    FixedField<double, 2> t;         // mesh points t_{i, j} = grid[i] + delta_t[i] * c[i][j] (c: collocation nodes for interval i)
    FixedVector<int>      nodes;     // number of collocation nodes p for each interval
    FixedField<int, 2>    acc_nodes; // number of nodes to the left of index (i, j)


    static Mesh createEquidistantMeshFixedDegree(int intervals, double tf, int p, Collocation& collocation); 
    FixedField<int, 2> createAccOffsetXU(int off_x, int off_xu);
    FixedField<int, 2> createAccOffsetFG(int off_fg);
};

// given some data trajectories t, x(t), u(t), p -> interpolate w.r.t. mesh and collocation scheme -> new fitting guess
struct Trajectory {
    std::vector<double> t;
    std::vector<std::vector<double>> x;
    std::vector<std::vector<double>> u;
    std::vector<double> p;
    InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    Trajectory interpolate(Mesh& mesh, Collocation& collocation);
    Trajectory linearInterpolation(Mesh& mesh, Collocation& collocation);
};

#endif  // OPT_MESH_H
