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
    grid[intervals] = tf;

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

std::vector<std::vector<int>> Mesh::createAccOffsetFG(int off_fg) {
    std::vector<std::vector<int>> acc_fg = acc_nodes;
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < nodes[i] ; j++) {
            acc_fg[i][j] *= off_fg;
        }
    }
    return acc_fg;
}

Trajectory Trajectory::interpolate(Mesh& mesh, Collocation& collocation) {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return linearInterpolation(mesh, collocation);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

Trajectory Trajectory::linearInterpolation(Mesh& mesh, Collocation& collocation) {
    Trajectory new_guess;

    std::vector<double> new_t = {0};
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            double new_time = mesh.grid[i] + mesh.delta_t[i] * collocation.c[mesh.nodes[i]][j];
            new_t.push_back(new_time);
        }
    }

    new_guess.t = new_t;
    new_guess.x.resize(x.size());
    new_guess.u.resize(u.size());

    // interpolate x(t)
    for (int k = 0; k < x.size(); k++) {
        new_guess.x[k].resize(new_t.size());
        for (int i = 0; i < new_t.size(); i++) {
            double t_new = new_t[i];
            auto it = std::lower_bound(t.begin(), t.end(), t_new);
            if (it == t.begin()) {
                new_guess.x[k][i] = x[k][0];
            }
            else if (it == t.end()) {
                new_guess.x[k][i] = x[k].back();
            }
            else {
                int idx = std::distance(t.begin(), it);
                double t1 = t[idx - 1];
                double t2 = t[idx];
                double x1 = x[k][idx - 1];
                double x2 = x[k][idx];
                new_guess.x[k][i] = x1 + (t_new - t1) * (x2 - x1) / (t2 - t1);
            }
        }
    }

    // interpolate u(t)
    for (int k = 0; k < u.size(); k++) {
        new_guess.u[k].resize(new_t.size());
        for (int i = 0; i < new_t.size(); i++) {
            double t_new = new_t[i];
            auto it = std::lower_bound(t.begin(), t.end(), t_new);
            if (it == t.begin()) {
                new_guess.u[k][i] = u[k][0];
            }
            else if (it == t.end()) {
                new_guess.u[k][i] = u[k].back();
            }
            else {
                int idx = std::distance(t.begin(), it);
                double t1 = t[idx - 1];
                double t2 = t[idx];
                double u1 = u[k][idx - 1];
                double u2 = u[k][idx];
                new_guess.u[k][i] = u1 + (t_new - t1) * (u2 - u1) / (t2 - t1);
            }
        }
    }

    // static parameters
    new_guess.p = p;

    return new_guess;
}
