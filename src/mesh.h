#ifndef OPT_MESH_H
#define OPT_MESH_H

#include "collocation.h"
#include "base/fixed_vector.h"


enum class InterpolationMethod {
    LINEAR = 0
};

struct Mesh {
    Mesh(int intervals, double tf, FixedVector<double> grid, FixedVector<double> delta_t, FixedVector<int> nodes, FixedField<int, 2> acc_nodes, int node_count)
        : intervals(intervals), tf(tf), grid(std::move(grid)), delta_t(std::move(delta_t)), nodes(std::move(nodes)),
          acc_nodes(std::move(acc_nodes)), node_count(node_count) {
    }

    int intervals;                 // number of intervals
    int node_count;                // number of collocation nodes (sum over all intervals)
    double tf;                     // final time
    FixedVector<double> grid;      // grid base points
    FixedVector<double> delta_t;   // step size h for each interval
    FixedVector<int>    nodes;     // number of collocation nodes p for each interval
    FixedField<int, 2> acc_nodes;  // number of nodes to the left of index (i, j)

    static Mesh createEquidistantMeshFixedDegree(int intervals, double tf, int p);
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
