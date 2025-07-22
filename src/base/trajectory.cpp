#include "trajectory.h"

bool Trajectory::compatible_with_mesh(const Mesh& mesh, const Collocation& collocation) const {
    // check size of time vector: expect mesh.node_count + 1 (t = 0.0 included)
    return check_time_compatibility(t, {x, u}, mesh, /* include_initial_time */ true);
}

Trajectory Trajectory::interpolate_onto_mesh(const Mesh& mesh, const Collocation& collocation) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return interpolate_onto_mesh_linear(mesh, collocation);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

Trajectory Trajectory::interpolate_onto_mesh_linear(const Mesh& mesh, const Collocation& collocation) const {
    Trajectory new_traj;

    std::vector<f64> new_t = {0};
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            f64 new_time = mesh.grid[i] + mesh.delta_t[i] * collocation.c[mesh.nodes[i]][j];
            new_t.push_back(new_time);
        }
    }

    new_traj.t = new_t;
    interpolate_linear_multiple(t, x, new_t, new_traj.x);
    interpolate_linear_multiple(t, u, new_t, new_traj.u);
    new_traj.p = p;

    return new_traj;
}

void Trajectory::print() {
    print_trajectory(t, {
        {"x", x},
        {"u", u}
    }, "p", p);
}

int Trajectory::to_csv(const std::string& filename) const {
    return write_trajectory_csv(filename, t, {
        {"x", x},
        {"u", u}
    }, "p", p);
}

FixedVector<f64> Trajectory::extract_initial_states() const {
    FixedVector<f64> x0(x.size());

    for (size_t i = 0; i < x.size(); i++) {
        x0[i] = x[i][0];
    }

    return x0;
}

FixedVector<f64> Trajectory::state_errors_inf_norm(const Trajectory& other) const {
    FixedVector<f64> max_abs_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            throw std::runtime_error("State trajectory length mismatch in Trajectory::state_errors_inf_norm.");
        }

        f64* max_err = &max_abs_errors[x_idx];
        for (size_t t_idx = 0; t_idx < x_traj_1.size(); t_idx++) {
            f64 diff = std::abs(x_traj_1[t_idx] - x_traj_2[t_idx]);
            if (diff > *max_err) {
                *max_err = diff;
            }
        }
    }

    return max_abs_errors;
}

FixedVector<f64> Trajectory::state_errors_2_norm(const Trajectory& other) const {
    FixedVector<f64> norm_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            throw std::runtime_error("State trajectory length mismatch in Trajectory::state_errors_2_norm.");
        }

        f64 sum_sq_diff = 0.0;
        for (size_t t_idx = 0; t_idx < x_traj_1.size(); t_idx++) {
            f64 diff = x_traj_1[t_idx] - x_traj_2[t_idx];
            sum_sq_diff += diff * diff;
        }
        norm_errors[x_idx] = std::sqrt(sum_sq_diff);
    }

    return norm_errors;
}

FixedVector<f64> Trajectory::state_errors_1_norm(const Trajectory& other) const {
    FixedVector<f64> norm_errors(x.size());

    for (size_t x_idx = 0; x_idx < x.size(); x_idx++) {
        const auto& x_traj_1 = x[x_idx];
        const auto& x_traj_2 = other.x[x_idx];

        if (x_traj_1.size() != x_traj_2.size()) {
            throw std::runtime_error("State trajectory length mismatch in Trajectory::state_errors_2_norm.");
        }

        f64 sum_abs_diff = 0.0;
        for (size_t t_idx = 0; t_idx < x_traj_1.size(); t_idx++) {
            sum_abs_diff += std::abs(x_traj_1[t_idx] - x_traj_2[t_idx]);
        }
        norm_errors[x_idx] = sum_abs_diff;
    }

    return norm_errors;
}

FixedVector<f64> Trajectory::state_errors(const Trajectory& other, Linalg::Norm norm) const {
    if (this->t.size() != other.t.size()) {
        throw std::runtime_error("Time vector size mismatch in Trajectory::state_errors.");
    }

    switch (norm) {
        case Linalg::Norm::NORM_INF:
            return state_errors_inf_norm(other);
        case Linalg::Norm::NORM_2:
            return state_errors_2_norm(other);
        case Linalg::Norm::NORM_1:
            return state_errors_1_norm(other);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

// === Control Trajectory ===

ControlTrajectory Trajectory::copy_extract_controls() const {
    ControlTrajectory controls_copy;
    controls_copy.t = t;                         // copies t from Trajectory
    controls_copy.u = u;                         // copies u from Trajectory
    controls_copy.interpolation = interpolation; // copy interpolation
    controls_copy.last_index = 0;                // initialize last_index as 0

    return controls_copy;
}

void ControlTrajectory::interpolate_at_linear(f64 t_query, f64* interpolation_values) const {
    const size_t t_len = t.size();

    // out of bounds cases
    if (t_query <= t.front()) {
        for (size_t k = 0; k < u.size(); ++k) {
            interpolation_values[k] = u[k][0];
        }
        last_index = 0;
        return;
    }
    if (t_query >= t.back()) {
        for (size_t k = 0; k < u.size(); ++k) {
            interpolation_values[k] = u[k].back();
        }
        last_index = t_len - 2; // safe last segment
        return;
    }

    // search for the correct interval using the last_index
    size_t i = last_index;
    while (i + 1 < t_len && t_query > t[i + 1]) {
        ++i;
    }
    while (i > 0 && t_query < t[i]) {
        --i;
    }

    // interval [t[i], t[i+1]]
    f64 t1 = t[i];
    f64 t2 = t[i + 1];
    f64 alpha = (t_query - t1) / (t2 - t1);

    for (size_t k = 0; k < u.size(); ++k) {
        f64 u1 = u[k][i];
        f64 u2 = u[k][i + 1];
        interpolation_values[k] = u1 + alpha * (u2 - u1);
    }

    last_index = i; // keep last_index for next call
    return;
}

void ControlTrajectory::interpolate_at(f64 t_query, f64* interpolation_values) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            interpolate_at_linear(t_query, interpolation_values);
            return;
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

// === Dual Trajectory ===

bool CostateTrajectory::compatible_with_mesh(const Mesh& mesh, const Collocation& collocation) const {
    // For duals, time grid has no t=0. So expect mesh.node_count entries (not +1)
    return check_time_compatibility(t, {costates_f, costates_g}, mesh, /* include_initial_time */ false);
}

CostateTrajectory CostateTrajectory::interpolate_onto_mesh(const Mesh& mesh, const Collocation& collocation) const {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            return interpolate_onto_mesh_linear(mesh, collocation);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

CostateTrajectory CostateTrajectory::interpolate_onto_mesh_linear(const Mesh& mesh, const Collocation& collocation) const {
    CostateTrajectory new_dual;

    std::vector<f64> new_t;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            f64 new_time = mesh.grid[i] + mesh.delta_t[i] * collocation.c[mesh.nodes[i]][j];
            new_t.push_back(new_time);
        }
    }

    new_dual.t = new_t;
    interpolate_linear_multiple(t, costates_f, new_t, new_dual.costates_f);
    interpolate_linear_multiple(t, costates_g, new_t, new_dual.costates_g);
    new_dual.costates_r = costates_r;

    return new_dual;
}

void CostateTrajectory::print() {
    print_trajectory(t, {
        {"costates_f", costates_f},
        {"costates_g", costates_g}
    }, "costates_r", costates_r);
}

int CostateTrajectory::to_csv(const std::string& filename) const {
    return write_trajectory_csv(filename, t, {
        {"costates_f", costates_f},
        {"costates_g", costates_g}
    }, "costates_r", costates_r);
}

// === helpers for Dual and standard Trajectory ===

bool check_time_compatibility(
    const std::vector<f64>& t_vec,
    const std::vector<std::vector<std::vector<f64>>>& fields_to_check,
    const Mesh& mesh,
    bool include_initial_time /* true for Trajectory, false for CostateTrajectory */ )
{
    const f64 tol = 1e-12;

    int expected_size = mesh.node_count + (include_initial_time ? 1 : 0);
    if ((int)t_vec.size() != expected_size) {
        return false;
    }

    int time_idx = include_initial_time ? 1 : 0;
    for (int i = 0; i < mesh.intervals; i++) {
        for (int j = 0; j < mesh.nodes[i]; j++) {
            if (std::abs(t_vec[time_idx++] - mesh.t[i][j]) > tol) {
                return false;
            }
        }
    }

    for (const auto& vect : fields_to_check) {
        for (const auto& field : vect) {
            if ((int)field.size() != (int)t_vec.size()) {
                return false;
            }
        }
    }

    return true;
}

void interpolate_linear_single(
    const std::vector<f64>& old_t,
    const std::vector<f64>& values,
    const std::vector<f64>& new_t,
    std::vector<f64>& out_values)
{
    out_values.resize(new_t.size());

    for (int i = 0; i < int(new_t.size()); i++) {
        f64 t_new = new_t[i];
        auto it = std::lower_bound(old_t.begin(), old_t.end(), t_new);

        if (it == old_t.begin()) {
            out_values[i] = values[0];
        } else if (it == old_t.end()) {
            out_values[i] = values.back();
        } else {
            int idx = std::distance(old_t.begin(), it);
            f64 t1 = old_t[idx - 1];
            f64 t2 = old_t[idx];
            f64 y1 = values[idx - 1];
            f64 y2 = values[idx];
            out_values[i] = y1 + (t_new - t1) * (y2 - y1) / (t2 - t1);
        }
    }
}

void interpolate_linear_multiple(
    const std::vector<f64>& old_t,
    const std::vector<std::vector<f64>>& values,
    const std::vector<f64>& new_t,
    std::vector<std::vector<f64>>& out_values)
{
    int fields = int(values.size());
    out_values.resize(fields);
    for (int field = 0; field < fields; field++) {
        interpolate_linear_single(old_t, values[field], new_t, out_values[field]);
    }
}

void print_trajectory(
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field)
{
    auto print_vector = [](const std::string& name, const std::vector<f64>& vec) {
        std::cout << name << " = [";
        for (size_t i = 0; i < vec.size(); ++i) {
            std::cout << vec[i];
            if (i + 1 < vec.size()) std::cout << ", ";
        }
        std::cout << "]\n";
    };

    auto print_matrix = [](const std::string& name, const std::vector<std::vector<f64>>& mat) {
        std::cout << name << " = [\n";
        for (const auto& row : mat) {
            std::cout << "  [";
            for (size_t j = 0; j < row.size(); ++j) {
                std::cout << row[j];
                if (j + 1 < row.size()) std::cout << ", ";
            }
            std::cout << "],\n";
        }
        std::cout << "]\n";
    };

    print_vector("t", t);
    for (const auto& [name, mat] : fields) {
        print_matrix(name, mat);
    }
    print_vector(static_name, static_field);
}

int write_trajectory_csv(
    const std::string& filename,
    const std::vector<f64>& t,
    const std::vector<std::pair<std::string, std::vector<std::vector<f64>>>>& fields,
    const std::string& static_name,
    const std::vector<f64>& static_field) // static_field should ideally have only one value or a set of static values
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[Warning] Failed to open file for writing: " << filename << "\n";
        return -1;
    }

    // header
    file << "time";
    for (const auto& [name, mat] : fields) {
        for (size_t i = 0; i < mat.size(); ++i) {
            file << "," << name << "[" << i << "]";
        }
    }

    // static field header(s)
    for (size_t i = 0; i < static_field.size(); ++i) {
        file << "," << static_name;
        if (static_field.size() > 1) {
            file << "[" << i << "]";
        }
    }
    file << "\n";

    file << std::setprecision(16);

    for (size_t k = 0; k < t.size(); ++k) {
        file << t[k];
        for (const auto& [_, mat] : fields) {
            for (const auto& series : mat) {
                if (k < series.size()) {
                    file << "," << series[k];
                } else {
                    file << ",";
                }
            }
        }

        for (size_t i = 0; i < static_field.size(); ++i) {
            file << "," << static_field[i];
        }
        file << "\n";
    }

    file.close();
    return 0;
}
