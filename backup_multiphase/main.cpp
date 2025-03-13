#include <iostream>

#include "integrator.h"
#include "util.h"
#include "problem.h"

//  cmake --build build && cd build && ./gdopt_experimental && cd ..

int main() {
    const Integrator fLGR = Integrator();
    auto fs = std::make_shared<FullSweep>();
    fs->xSize = 5;
    auto phase = Phase();
    phase.fullSweep = fs;
    auto phase2 = Phase();
    phase2.fullSweep = fs;
    phase2.fullSweep->xSize++;
    std::cout << phase.fullSweep->xSize << "\n";
    return 0;
}