#include <iostream>

#include "collocation.h"
#include "util.h"
#include "nlp.h"
#include "problem.h"

//  cmake --build build && cd build && ./gdopt_experimental && cd ..

int main() {
    const Collocation fLGR = Collocation();
    Mesh mesh = Mesh::createEquidistantMeshFixedDegree(10, 1, 3);
    Problem model = Problem();
    model.xSize = 3;
    model.uSize = 2;
    model.pSize = 1;
    std::cout << "Y" << std::endl;
    NLP nlp = NLP(model, fLGR, mesh);
    for (const auto& v : nlp.off_acc_xu) {
    std::cout << Util::vectorToString(v) << std::endl;

    }
    std::cout << nlp.number_vars << std::endl;
    return 0;
}