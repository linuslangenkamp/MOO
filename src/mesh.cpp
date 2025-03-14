#include "mesh.h"

/** 
 * @brief Create a mesh with fixed (h, p) for final time tf
 *
 * @param intervals  Number of vntervals
 * @param tf         Final Time
 * @param p          Fixed Polynomial Degree
 * @return Mesh      Mesh
 */
Mesh Mesh::createEquidistantMeshFixedDegree(int intervals, double tf, int p) {
    std::vector<double> grid(intervals + 1);
    std::vector<double> deltaT(intervals);
    std::vector<int> steps(intervals);
    std::vector<int> cum_steps(intervals);

    double h = tf / intervals;
    for (int i = 0; i < intervals; i++) {
        grid[i] = i * h;
    }
    grid[intervals - 1] = tf;

    for (int i = 0; i < intervals; i++) {
        deltaT[i] = h;
        steps[i] = p;
        cum_steps[i] = p * (i + 1);
    }

    return {intervals, tf, std::move(grid), std::move(deltaT), std::move(steps), std::move(cum_steps)};
}

std::vector<std::vector<int>> Mesh::createAccOffsetXU(int off_x, int off_xu) {
    std::vector<std::vector<int>> off_acc_xu(intervals);
    int off = off_x;
    for (int i = 0; i < intervals; i++) {
        std::vector<int> local; 
        for (int j = 0; j < steps[i]; j++) {
            local.push_back(off);
            off += off_xu;
        }
        off_acc_xu.push_back(local);
    }
    return off_acc_xu;
}
