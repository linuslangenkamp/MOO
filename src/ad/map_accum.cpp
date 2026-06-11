// SPDX-License-Identifier: LGPL-3.0-or-later
#include "map_accum.h"

#include "detail/function_core.h"
#include "detail/mapping.h"
#include "detail/node.h"
#include "detail/simplify.h"

#include <set>
#include <stdexcept>
#include <utility>

namespace ad {
namespace {

const detail::FunctionCore &require_step_core(const Function &step) {
    const std::shared_ptr<const detail::FunctionCore> &core = detail::function_core(step);
    if (!core) {
        throw std::runtime_error("map-accum step function must be valid");
    }
    return *core;
}

int find_carry_input_index(const detail::FunctionCore &core, const Vec &carry_input) {
    if (!carry_input.valid()) {
        throw std::runtime_error("map-accum carry input must be a valid vector");
    }
    for (int i = 0; i < core.info.input_count; ++i) {
        if (core.inputs[static_cast<std::size_t>(i)].node_id() == carry_input.node_id()) {
            return i;
        }
    }
    throw std::runtime_error("map-accum carry input must be one of the step inputs");
}

const Binding *find_binding_for_input(const std::vector<Binding> &bindings, NodeId input_node_id) {
    const Binding *found = nullptr;
    for (const Binding &binding : bindings) {
        if (!binding.valid()) {
            throw std::runtime_error("map-accum binding must be valid");
        }
        if (binding.local_input().node_id() == input_node_id) {
            if (found != nullptr) {
                throw std::runtime_error("map-accum step input is bound more than once");
            }
            found = &binding;
        }
    }
    return found;
}

std::vector<int> stored_indices(const Map &map) {
    if (map.stores_expanded_indices()) {
        return map.indices();
    }
    return {};
}

detail::MappedBindingNode mapped_binding_from_public(const Vec &local_input, const Map &map) {
    detail::MappedBindingNode binding;
    binding.local_input = detail::vec_node(local_input);
    binding.source = detail::vec_node(map.source());
    binding.reps = map.reps();
    binding.local_size = map.local_size();
    binding.indices = stored_indices(map);
    binding.map_kind = map.kind();
    binding.base = map.base();
    binding.rep_stride = map.rep_stride();
    binding.component_stride = map.component_stride();
    binding.shift = map.shift();
    binding.offsets = map.offsets();
    return binding;
}

std::vector<detail::MappedBindingNode> ordered_sequence_bindings(const detail::FunctionCore &core,
                                                                 int carry_input_index,
                                                                 const std::vector<Binding> &bindings,
                                                                 int reps) {
    if (static_cast<int>(bindings.size()) != core.info.input_count - 1) {
        throw std::runtime_error("map-accum must bind every non-carry step input exactly once");
    }

    std::set<NodeId> input_node_ids;
    for (const FunctionInputInfo &input : core.info.inputs) {
        input_node_ids.insert(input.node_id);
    }
    for (const Binding &binding : bindings) {
        if (!binding.valid()) {
            throw std::runtime_error("map-accum binding must be valid");
        }
        if (binding.local_input().node_id() == core.inputs[static_cast<std::size_t>(carry_input_index)].node_id()) {
            throw std::runtime_error("map-accum carry input must not be sequence-bound");
        }
        if (input_node_ids.find(binding.local_input().node_id()) == input_node_ids.end()) {
            throw std::runtime_error("map-accum binding local input is not a step input group");
        }
    }

    std::vector<detail::MappedBindingNode> ordered;
    ordered.reserve(bindings.size());
    for (int i = 0; i < core.info.input_count; ++i) {
        if (i == carry_input_index) {
            continue;
        }
        const Vec &input = core.inputs[static_cast<std::size_t>(i)];
        const FunctionInputInfo &input_info = core.info.inputs[static_cast<std::size_t>(i)];
        const Binding *binding = find_binding_for_input(bindings, input.node_id());
        if (binding == nullptr) {
            throw std::runtime_error("map-accum is missing a binding for a non-carry step input");
        }
        if (binding->map().local_size() != input_info.size) {
            throw std::runtime_error("map-accum binding local size does not match step input size");
        }
        if (binding->map().reps() != reps) {
            throw std::runtime_error("map-accum binding repetition count does not match map-accum reps");
        }
        ordered.push_back(mapped_binding_from_public(input, binding->map()));
    }
    return ordered;
}

int emit_offset(const detail::FunctionCore &core, int reps, int output_index) {
    const int carry_size = core.info.outputs[0].size;
    int offset = carry_size + reps * carry_size;
    for (int i = 1; i < output_index; ++i) {
        offset += reps * core.info.outputs[static_cast<std::size_t>(i)].size;
    }
    return offset;
}

} // namespace

MapAccumResult map_accum(const Function &step,
                         const Vec &carry_input,
                         const Vec &initial_carry,
                         int reps,
                         std::vector<Binding> sequence_bindings) {
    if (reps <= 0) {
        throw std::runtime_error("map-accum repetition count must be positive");
    }
    const detail::FunctionCore &core = require_step_core(step);
    const int carry_input_index = find_carry_input_index(core, carry_input);
    std::vector<detail::MappedBindingNode> ordered = ordered_sequence_bindings(core, carry_input_index, sequence_bindings, reps);
    Vec flat = detail::map_accum_call_from_bindings(detail::function_core(step),
                                                    carry_input_index,
                                                    detail::simplify_vec(initial_carry),
                                                    std::move(ordered),
                                                    reps);

    const int carry_size = core.info.outputs[0].size;
    MapAccumResult result;
    result.final_carry = flat.slice(0, carry_size);
    result.carry_trajectory = flat.slice(carry_size, reps * carry_size);
    result.outputs.reserve(core.info.outputs.size() > 0 ? core.info.outputs.size() - 1 : 0);
    for (int i = 1; i < core.info.output_count; ++i) {
        const int size = core.info.outputs[static_cast<std::size_t>(i)].size;
        result.outputs.push_back(flat.slice(emit_offset(core, reps, i), reps * size));
    }
    return result;
}

Vec fold(const Function &step,
         const Vec &carry_input,
         const Vec &initial_carry,
         int reps,
         std::vector<Binding> sequence_bindings) {
    return map_accum(step, carry_input, initial_carry, reps, std::move(sequence_bindings)).final_carry;
}

Vec scan(const Function &step,
         const Vec &carry_input,
         const Vec &initial_carry,
         int reps,
         std::vector<Binding> sequence_bindings) {
    return map_accum(step, carry_input, initial_carry, reps, std::move(sequence_bindings)).carry_trajectory;
}

} // namespace ad
