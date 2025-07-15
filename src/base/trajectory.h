#ifndef OPT_TRAJECTORY_H
#define OPT_TRAJECTORY_H

#include "mesh.h"
#include "linalg.h"
#include "collocation.h"

enum class InterpolationMethod {
    LINEAR = 0
};

struct ControlTrajectory {
    std::vector<f64> t;                       // time grid, monotonic increasing
    std::vector<std::vector<f64>> u;          // u[k][j] = value of k-th control at t[j]
    InterpolationMethod interpolation;

    // for repeated interpolation, cache last index
    mutable size_t last_index = 0;

    void interpolate_at(f64 t_query, f64* interpolation_values) const;
    void interpolate_at_linear(f64 t_query, f64* interpolation_values) const;
};

// given some data trajectories t, x(t), u(t), p -> interpolate w.r.t. mesh and collocation scheme -> new fitting guess
struct Trajectory {
    // x[i][j] = x_i(t_j) = x_i at time t[j] 
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

    // create new trajectories based on mesh & collocation
    Trajectory interpolate_onto_mesh(Mesh& mesh, Collocation& collocation);
    Trajectory interpolate_onto_mesh_linear(Mesh& mesh, Collocation& collocation);

    // extract + copy information from the trajectory
    ControlTrajectory copy_extract_controls() const;
    FixedVector<f64> extract_initial_states() const; 

    // compare with other trajectory
    FixedVector<f64> state_errors(const Trajectory& other, Linalg::Norm norm) const;
    FixedVector<f64> state_errors_inf_norm(const Trajectory& other) const;

    // dumps
    void print();
    int to_csv(const std::string& filename) const;
};

#endif // OPT_TRAJECTORY_H
