// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_DETAIL_MAPPING_H
#define MOO_AD_DETAIL_MAPPING_H

#include "../vec.h"
#include "node.h"

#include <memory>
#include <vector>

namespace ad::detail {

struct FunctionCore;

int mapped_index(const MappedBindingNode &binding, int reps, int rep, int local_component);
std::vector<int> materialize_indices(const MappedBindingNode &binding, int reps);

Vec mapped_call_from_bindings(std::shared_ptr<const FunctionCore> function,
                              std::vector<MappedBindingNode> bindings,
                              int reps,
                              MappedOutput output = MappedOutput::concat());

} // namespace ad::detail

#endif // MOO_AD_DETAIL_MAPPING_H
