#ifndef OPT_MESH_H
#define OPT_MESH_H

#include <vector>

struct Mesh {
    Mesh(int intervals, double tf, std::vector<double> grid, std::vector<double> deltaT, std::vector<int> steps, std::vector<int> cum_steps)
        : intervals(intervals), tf(tf), grid(std::move(grid)), deltaT(std::move(deltaT)), steps(std::move(steps)), cum_steps(std::move(cum_steps))  {
    }

    int intervals;
    double tf;
    std::vector<double> grid;   // grid base points
    std::vector<double> deltaT; // step size h for each interval
    std::vector<int> steps;     // number of collocation nodes p for each interval
    std::vector<int> cum_steps; // cum_steps[k] = sum_{i=0}^{k} steps[i] -> cum_steps[-1] = # collocation nodes in total

    static Mesh createEquidistantMeshFixedDegree(int intervals, double tf, int p);
    std::vector<std::vector<int>> createAccOffsetXU(int off_x, int off_xu);
};

#endif  // OPT_MESH_H
