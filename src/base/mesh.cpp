#include "mesh.h"

/** 
 * @brief Create a mesh with fixed (h, p) for final time tf
 *
 * @param intervals  Number of Intervals
 * @param tf         Final Time
 * @param stages     Number of Nodes for each Interval
 * @return Mesh      Mesh
 */
Mesh Mesh::create_equidistant_fixed_stages(double tf, int intervals, int stages, Collocation& collocation) {
    FixedVector<double> grid(intervals + 1);
    FixedVector<double> delta_t(intervals);
    FixedVector<int> nodes(intervals);
    FixedField<int, 2> acc_nodes(intervals, stages);
    FixedField<double, 2> t(intervals, stages);

    double h = tf / intervals;
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

FixedField<int, 2> Mesh::create_acc_offset_xu(int off_x, int off_xu) {
    FixedField<int, 2> off_acc_xu(intervals);
    int off = off_x;
    for (int i = 0; i < intervals; i++) {
        off_acc_xu[i] = FixedVector<int>(nodes[i]);
        for (int j = 0; j < nodes[i]; j++) {
            off_acc_xu[i][j] = off;
            off += off_xu;
        }
    }
    return off_acc_xu;
}

FixedField<int, 2> Mesh::create_acc_offset_fg(int off_fg) {
    FixedField<int, 2> acc_fg = acc_nodes;
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
            return linear_interpolation(mesh, collocation);
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}

Trajectory Trajectory::linear_interpolation(Mesh& mesh, Collocation& collocation) {
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
    for (int k = 0; k < int_size(x); k++) {
        new_guess.x[k].resize(new_t.size());
        for (int i = 0; i < int_size(new_t); i++) {
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
    for (int k = 0; k < int_size(u); k++) {
        new_guess.u[k].resize(new_t.size());
        for (int i = 0; i < int_size(new_t); i++) {
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

void Trajectory::print() {
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
    print_matrix("x", x);
    print_matrix("u", u);
    print_vector("p", p);
}

void Trajectory::to_csv(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "[Warning] Failed to open file for writing: " << filename << "\n";
        return;
    }

    // header
    file << "time";
    for (size_t i = 0; i < x.size(); ++i)
        file << ",x[" << i << "]";
    for (size_t i = 0; i < u.size(); ++i)
        file << ",u[" << i << "]";
    file << "\n";

    size_t num_rows = t.size();
    file << std::setprecision(16);

    for (size_t k = 0; k < num_rows; ++k) {
        file << t[k];
        for (size_t i = 0; i < x.size(); ++i)
            file << "," << x[i][k];
        for (size_t i = 0; i < u.size(); ++i)
            file << "," << u[i][k];
        file << "\n";
    }

    file.close();
}

ControlTrajectory Trajectory::copy_extract_controls() {
    ControlTrajectory controls_copy;
    controls_copy.t = t;                         // copies t from Trajectory
    controls_copy.u = u;                         // copies u from Trajectory
    controls_copy.interpolation = interpolation; // copy interpolation
    controls_copy.last_index = 0;                // initialize last_index as 0

    return controls_copy;
}

void ControlTrajectory::linear_interpolation(f64 t_query, f64* interpolation_values) const {
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

void ControlTrajectory::interpolate(f64 t_query, f64* interpolation_values) {
    switch (interpolation) {
        case InterpolationMethod::LINEAR:
            linear_interpolation(t_query, interpolation_values);
            return;
        default:
            throw std::runtime_error("Unknown interpolation method!");
    }
}
