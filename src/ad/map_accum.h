// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_MAP_ACCUM_H
#define MOO_AD_MAP_ACCUM_H

#include "function.h"
#include "map.h"

#include <vector>

namespace ad {

struct MapAccumResult {
    Vec final_carry;
    Vec carry_trajectory;
    std::vector<Vec> outputs;
};

MapAccumResult map_accum(const Function &step,
                         const Vec &carry_input,
                         const Vec &initial_carry,
                         int reps,
                         std::vector<Binding> sequence_bindings);

Vec fold(const Function &step,
         const Vec &carry_input,
         const Vec &initial_carry,
         int reps,
         std::vector<Binding> sequence_bindings);

Vec scan(const Function &step,
         const Vec &carry_input,
         const Vec &initial_carry,
         int reps,
         std::vector<Binding> sequence_bindings);

} // namespace ad

#endif // MOO_AD_MAP_ACCUM_H
