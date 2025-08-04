#ifndef OPT_TRAJECTORY_H
#define OPT_TRAJECTORY_H

#include "log.h"
#include "mesh.h"
#include "linalg.h"
#include "fLGR.h"

enum class InterpolationMethod {
    LINEAR = 0
};
// TODO: refactor this interpolation garbage with an interpolator that optionally holds a t already -> only copy!

struct ControlTrajectory {
    std::vector<f64> t;                       // time grid, monotonic increasing
    std::vector<std::vector<f64>> u;          // u[k][j] = value of k-th control at t[j]
    InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    // for repeated interpolation, cache last index
    mutable size_t last_index = 0;

    void interpolate_at(f64 t_query, f64* interpolation_values) const;
    void interpolate_at_linear(f64 t_query, f64* interpolation_values) const;
};

// given some data trajectories t, x(t), u(t), p
struct Trajectory {
    // x[i][j] = x_i(t_j) = x_i at time t[j] 
    std::vector<f64> t;
    std::vector<std::vector<f64>> x;
    std::vector<std::vector<f64>> u;
    std::vector<f64> p;
    InterpolationMethod interpolation = InterpolationMethod::LINEAR;

    Trajectory() = default;

    Trajectory(std::vector<f64> t, std::vector<std::vector<f64>> x, std::vector<std::vector<f64>> u,
               std::vector<f64> p, InterpolationMethod interpolation = InterpolationMethod::LINEAR)
        : t(t), x(x), u(u), p(p), interpolation(interpolation) {
    }

    Trajectory(const Trajectory& other)
        : t(other.t),
          x(other.x),
          u(other.u),
          p(other.p),
          interpolation(other.interpolation) {}

    // checks if a trajectory has the same nodes as a mesh + collocation (valid Trajectory + fLGR)
    bool compatible_with_mesh(const Mesh& mesh) const;

    // create new trajectories based on mesh & collocation
    Trajectory interpolate_onto_mesh(const Mesh& mesh) const;
    Trajectory interpolate_onto_mesh_linear(const Mesh& mesh) const;

    Trajectory interpolate_polynomial_from_mesh_onto_grid(const Mesh& mesh,
                                                          const std::vector<f64>& time_grid);

    // extract + copy information from the trajectory
    ControlTrajectory copy_extract_controls() const;
    FixedVector<f64> extract_initial_states() const; 

    // compare with other trajectory
    FixedVector<f64> state_errors(const Trajectory& other, Linalg::Norm norm) const;
    FixedVector<f64> state_errors_inf_norm(const Trajectory& other) const;
    FixedVector<f64> state_errors_2_norm(const Trajectory& other) const;
    FixedVector<f64> state_errors_1_norm(const Trajectory& other) const;

    // dumps
    void print();
    int to_csv(const std::string& filename) const;
};

// dual trajectory for [costates_f, costates_g]_{ij} constraints, costates_r constraints
struct CostateTrajectory {
    // costates_f[i][j] = costates_f_i(t_j) = costates_f_i at time t[j] 
    std::vector<f64> t;
    std::vector<std::vector<f64>> costates_f;
    std::vector<std::vector<f64>> costates_g;
    std::vector<f64> costates_r;
    InterpolationMethod interpolation;

    CostateTrajectory() = default;

    CostateTrajectory(std::vector<f64> t, std::vector<std::vector<f64>> costates_f, std::vector<std::vector<f64>> costates_g,
                   std::vector<f64> costates_r, InterpolationMethod interpolation = InterpolationMethod::LINEAR)
        : t(t), costates_f(costates_f), costates_g(costates_g), costates_r(costates_r), interpolation(interpolation) {
    }

    CostateTrajectory(const CostateTrajectory& other)
        : t(other.t),
          costates_f(other.costates_f),
          costates_g(other.costates_g),
          costates_r(other.costates_r),
          interpolation(other.interpolation) {}

    // checks if a dual trajectory has the same nodes as a mesh + collocation (valid CostateTrajectory + fLGR)
    bool compatible_with_mesh(const Mesh& mesh) const;

    // create new trajectories based on mesh & collocation
    CostateTrajectory interpolate_onto_mesh(const Mesh& mesh) const;
    CostateTrajectory interpolate_onto_mesh_linear(const Mesh& mesh) const;

    // dumps
    void print();
    int to_csv(const std::string& filename) const;
};

struct PrimalDualTrajectory {
    std::unique_ptr<Trajectory> primals;          // unscaled primal variables: x
    std::unique_ptr<CostateTrajectory> costates;  // transformed constraint multipliers / costates: \hat{lambda_g}
    std::unique_ptr<Trajectory> lower_costates;   // transformed (lb) variable bound multipliers / costates: \hat{x_L}
    std::unique_ptr<Trajectory> upper_costates;   // transformed (ub) variable bound multipliers / costates: \hat{x_U}

    // full constructor
    PrimalDualTrajectory(std::unique_ptr<Trajectory> primals_,
                         std::unique_ptr<CostateTrajectory> costates_,
                         std::unique_ptr<Trajectory> lower_costates_,
                         std::unique_ptr<Trajectory> upper_costates_)
        : primals(std::move(primals_)),
          costates(std::move(costates_)),
          lower_costates(std::move(lower_costates_)),
          upper_costates(std::move(upper_costates_)) {}

    // primals + costates constructor
    PrimalDualTrajectory(std::unique_ptr<Trajectory> primals_,
                         std::unique_ptr<CostateTrajectory> costates_)
        : primals(std::move(primals_)),
          costates(std::move(costates_)),
          lower_costates(nullptr),
          upper_costates(nullptr) {}

    // primals only constructor
    PrimalDualTrajectory(std::unique_ptr<Trajectory> primals_)
        : primals(std::move(primals_)),
          costates(nullptr),
          lower_costates(nullptr),
          upper_costates(nullptr) {}
};

// === shared helpers for Trajectory and CostateTrajectory ===

std::vector<f64> interpolate_polynomial_from_mesh_onto_grid_single(
    const Mesh& mesh,
    const std::vector<f64>& values,
    const std::vector<f64>& time_grid);

std::vector<std::vector<f64>> interpolate_polynomial_from_mesh_onto_grid_multiple(
    const Mesh& mesh,
    const std::vector<std::vector<f64>>& values,
    const std::vector<f64>& time_grid);

std::vector<f64> interpolate_linear_single(
    const std::vector<f64>& t,
    const std::vector<f64>& values,
    const std::vector<f64>& new_t);

std::vector<std::vector<f64>> interpolate_linear_multiple(
    const std::vector<f64>& t,
    const std::vector<std::vector<f64>>& values,
    const std::vector<f64>& new_t);

bool check_time_compatibility(
    const std::vector<f64>& t_vec,
    const std::vector<std::vector<std::vector<f64>>>& fields_to_check,
    const Mesh& mesh);

int write_trajectory_csv(
    const std::string& filename,
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field);

void print_trajectory(
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field);

#endif // OPT_TRAJECTORY_H
