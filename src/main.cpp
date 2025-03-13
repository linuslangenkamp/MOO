#include <iostream>

#include "collocation.h"
#include "util.h"
#include "problem.h"

//  cmake --build build && cd build && ./gdopt_experimental && cd ..

int main() {
    const Collocation fLGR = Collocation();
    auto fs = FullSweep();
    auto model = Problem();
    model.fullSweep = fs;
    auto model2 = Problem();
    model2.fullSweep = fs;
    model2.xSize++;
    std::cout << model.xSize << "\n";
    return 0;
}