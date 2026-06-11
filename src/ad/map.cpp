// SPDX-License-Identifier: LGPL-3.0-or-later
#include "map.h"

#include "detail/function_core.h"
#include "detail/mapping.h"
#include "detail/node.h"
#include "detail/simplify.h"
#include "detail/traversal.h"

#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ad {
namespace {

void require_valid_source(const Vec &source) {
    if (!source.valid()) {
        throw std::runtime_error("map source must be a valid vector");
    }
}

void require_positive_reps(int reps) {
    if (reps <= 0) {
        throw std::runtime_error("map repetition count must be positive");
    }
}

void require_nonnegative_local_size(int local_size) {
    if (local_size < 0) {
        throw std::runtime_error("map local size must be non-negative");
    }
}

std::size_t table_size(int reps, int local_size) {
    require_positive_reps(reps);
    require_nonnegative_local_size(local_size);
    if (local_size != 0 && static_cast<std::size_t>(reps) > std::numeric_limits<std::size_t>::max() / static_cast<std::size_t>(local_size)) {
        throw std::runtime_error("map index table size overflow");
    }
    return static_cast<std::size_t>(reps) * static_cast<std::size_t>(local_size);
}

int checked_int(long long value, const char *message) {
    if (value < static_cast<long long>(std::numeric_limits<int>::min()) ||
        value > static_cast<long long>(std::numeric_limits<int>::max())) {
        throw std::runtime_error(message);
    }
    return static_cast<int>(value);
}

int checked_product(int lhs, int rhs, const char *message) {
    return checked_int(static_cast<long long>(lhs) * static_cast<long long>(rhs), message);
}

int checked_sum(long long lhs, long long rhs, const char *message) {
    return checked_int(lhs + rhs, message);
}

int map_index_from_descriptor(MapKind kind,
                              int reps,
                              int local_size,
                              const std::vector<int> &indices,
                              int base,
                              int rep_stride,
                              int component_stride,
                              int shift,
                              const std::vector<int> &offsets,
                              int rep,
                              int local_component) {
    if (rep < 0 || rep >= reps || local_component < 0 || local_component >= local_size) {
        throw std::out_of_range("map index coordinates are out of range");
    }

    switch (kind) {
        case MapKind::Single:
            if (rep != 0 || indices.size() != static_cast<std::size_t>(local_size)) {
                throw std::runtime_error("single map descriptor is invalid");
            }
            return indices[static_cast<std::size_t>(local_component)];
        case MapKind::Blocks:
            return checked_int(static_cast<long long>(rep) * static_cast<long long>(local_size) + local_component,
                               "block map index overflow");
        case MapKind::ShiftedBlocks: {
            const int shifted_rep = checked_sum(rep, shift, "shifted block map index overflow");
            return checked_int(static_cast<long long>(shifted_rep) * static_cast<long long>(local_size) + local_component,
                               "shifted block map index overflow");
        }
        case MapKind::Stride:
            return checked_int(static_cast<long long>(base) +
                                   static_cast<long long>(rep) * static_cast<long long>(rep_stride) +
                                   static_cast<long long>(local_component) * static_cast<long long>(component_stride),
                               "stride map index overflow");
        case MapKind::Stencil:
            if (offsets.size() != static_cast<std::size_t>(local_size)) {
                throw std::runtime_error("stencil map descriptor is invalid");
            }
            return checked_int(static_cast<long long>(base) +
                                   static_cast<long long>(rep) * static_cast<long long>(rep_stride) +
                                   static_cast<long long>(offsets[static_cast<std::size_t>(local_component)]),
                               "stencil map index overflow");
        case MapKind::ExplicitIndices:
        case MapKind::Table: {
            const std::size_t expected_size = table_size(reps, local_size);
            if (indices.size() != expected_size) {
                throw std::runtime_error("map index table size does not match reps * local_size");
            }
            return indices[static_cast<std::size_t>(rep * local_size + local_component)];
        }
        case MapKind::Invalid:
            throw std::runtime_error("map descriptor is invalid");
    }

    throw std::runtime_error("unsupported map descriptor");
}

std::vector<int> materialize_indices_from_descriptor(MapKind kind,
                                                     int reps,
                                                     int local_size,
                                                     const std::vector<int> &indices,
                                                     int base,
                                                     int rep_stride,
                                                     int component_stride,
                                                     int shift,
                                                     const std::vector<int> &offsets) {
    std::vector<int> out;
    out.reserve(table_size(reps, local_size));
    for (int rep = 0; rep < reps; ++rep) {
        for (int component = 0; component < local_size; ++component) {
            out.push_back(map_index_from_descriptor(kind, reps, local_size, indices, base, rep_stride, component_stride, shift, offsets, rep, component));
        }
    }
    return out;
}

void validate_descriptor_indices(const Vec &source,
                                 int reps,
                                 int local_size,
                                 const std::vector<int> &indices,
                                 MapKind kind,
                                 int base,
                                 int rep_stride,
                                 int component_stride,
                                 int shift,
                                 const std::vector<int> &offsets) {
    (void)table_size(reps, local_size);
    if (kind == MapKind::Single && reps != 1) {
        throw std::runtime_error("single map must have exactly one repetition");
    }
    if (kind == MapKind::Single && indices.size() != static_cast<std::size_t>(local_size)) {
        throw std::runtime_error("single map index count does not match local size");
    }
    if ((kind == MapKind::ExplicitIndices || kind == MapKind::Table) && indices.size() != table_size(reps, local_size)) {
        throw std::runtime_error("map index table size does not match reps * local_size");
    }
    if (kind == MapKind::Stencil && offsets.size() != static_cast<std::size_t>(local_size)) {
        throw std::runtime_error("stencil map offsets size does not match local size");
    }
    for (int rep = 0; rep < reps; ++rep) {
        for (int component = 0; component < local_size; ++component) {
            const int index = map_index_from_descriptor(kind, reps, local_size, indices, base, rep_stride, component_stride, shift, offsets, rep, component);
            if (index < 0 || index >= source.size()) {
                throw std::runtime_error("map index is out of source vector range");
            }
        }
    }
}

const detail::FunctionCore &require_function_core(const Function &function) {
    const std::shared_ptr<const detail::FunctionCore> &core = detail::function_core(function);
    if (!core) {
        throw std::runtime_error("cannot map an invalid function");
    }
    return *core;
}

const Binding *find_binding_for_input(const std::vector<Binding> &bindings, NodeId input_node_id) {
    const Binding *found = nullptr;
    for (const Binding &binding : bindings) {
        if (!binding.valid()) {
            throw std::runtime_error("mapped function binding must be valid");
        }
        if (!binding.local_input().valid()) {
            throw std::runtime_error("mapped function binding local input must be valid");
        }
        if (binding.local_input().node_id() == input_node_id) {
            if (found != nullptr) {
                throw std::runtime_error("mapped function input is bound more than once");
            }
            found = &binding;
        }
    }
    return found;
}

void reject_bindings_to_non_inputs(const detail::FunctionCore &function, const std::vector<Binding> &bindings) {
    std::set<NodeId> input_node_ids;
    for (const FunctionInputInfo &input : function.info.inputs) {
        input_node_ids.insert(input.node_id);
    }

    for (const Binding &binding : bindings) {
        if (!binding.valid()) {
            throw std::runtime_error("mapped function binding must be valid");
        }
        if (input_node_ids.find(binding.local_input().node_id()) == input_node_ids.end()) {
            throw std::runtime_error("mapped function binding local input is not a callee input group");
        }
    }
}

int mapped_output_size(const MappedOutput &output, int reps, int local_output_size) {
    switch (output.mode()) {
        case OutputMode::Concat:
            return checked_product(reps, local_output_size, "mapped concat output size overflow");
        case OutputMode::Scatter:
            if (output.output_size() < 0) {
                throw std::runtime_error("mapped scatter output size must be non-negative");
            }
            if (output.indices().size() != table_size(reps, local_output_size)) {
                throw std::runtime_error("mapped scatter output index count must equal reps * local_output_size");
            }
            for (int index : output.indices()) {
                if (index < 0 || index >= output.output_size()) {
                    throw std::runtime_error("mapped scatter output index is out of output range");
                }
            }
            return output.output_size();
        case OutputMode::Sum:
            return local_output_size;
        case OutputMode::WeightedSum:
            if (output.weights().size() != static_cast<std::size_t>(reps)) {
                throw std::runtime_error("mapped weighted sum weight count must equal reps");
            }
            return local_output_size;
    }

    throw std::runtime_error("unsupported mapped output mode");
}

} // namespace

MappedOutput MappedOutput::concat() {
    return MappedOutput{};
}

MappedOutput MappedOutput::scatter(std::vector<int> indices, int output_size) {
    MappedOutput output;
    output.mode_ = OutputMode::Scatter;
    output.output_size_ = output_size;
    output.indices_ = std::move(indices);
    return output;
}

MappedOutput MappedOutput::sum() {
    MappedOutput output;
    output.mode_ = OutputMode::Sum;
    return output;
}

MappedOutput MappedOutput::weighted_sum(std::vector<double> weights) {
    MappedOutput output;
    output.mode_ = OutputMode::WeightedSum;
    output.weights_ = std::move(weights);
    return output;
}

OutputMode MappedOutput::mode() const {
    return mode_;
}

int MappedOutput::output_size() const {
    return output_size_;
}

const std::vector<int> &MappedOutput::indices() const {
    return indices_;
}

const std::vector<double> &MappedOutput::weights() const {
    return weights_;
}

Map::Map(Vec source,
         int reps,
         int local_size,
         std::vector<int> indices,
         MapKind kind,
         int base,
         int rep_stride,
         int component_stride,
         int shift,
         std::vector<int> offsets)
    : source_(std::move(source)),
      reps_(reps),
      local_size_(local_size),
      indices_(std::move(indices)),
      kind_(kind),
      base_(base),
      rep_stride_(rep_stride),
      component_stride_(component_stride),
      shift_(shift),
      offsets_(std::move(offsets)) {
    require_valid_source(source_);
    source_ = detail::simplify_vec(source_);
    validate_descriptor_indices(source_, reps_, local_size_, indices_, kind_, base_, rep_stride_, component_stride_, shift_, offsets_);
}

Map Map::single(const Vec &source, std::vector<int> indices) {
    const int local_size = static_cast<int>(indices.size());
    return Map(source, 1, local_size, std::move(indices), MapKind::Single, 0, 0, 1, 0, {});
}

Map Map::blocks(const Vec &source, int reps, int block_size) {
    require_positive_reps(reps);
    if (block_size <= 0) {
        throw std::runtime_error("block map block size must be positive");
    }
    return Map(source, reps, block_size, {}, MapKind::Blocks, 0, block_size, 1, 0, {});
}

Map Map::shifted_blocks(const Vec &source, int reps, int block_size, int shift) {
    require_positive_reps(reps);
    if (block_size <= 0) {
        throw std::runtime_error("block map block size must be positive");
    }
    return Map(source, reps, block_size, {}, MapKind::ShiftedBlocks, 0, block_size, 1, shift, {});
}

Map Map::stride(const Vec &source, int reps, int local_size, int base, int rep_stride, int component_stride) {
    return Map(source,
               reps,
               local_size,
               {},
               MapKind::Stride,
               base,
               rep_stride,
               component_stride,
               0,
               {});
}

Map Map::stencil(const Vec &source, int reps, int base, int rep_stride, std::vector<int> offsets) {
    const int local_size = static_cast<int>(offsets.size());
    return Map(source, reps, local_size, {}, MapKind::Stencil, base, rep_stride, 1, 0, std::move(offsets));
}

Map Map::explicit_indices(const Vec &source, int reps, int local_size, std::vector<int> flat_indices) {
    return Map(source, reps, local_size, std::move(flat_indices), MapKind::ExplicitIndices, 0, 0, 1, 0, {});
}

Map Map::table(const Vec &source, int reps, int local_size, std::vector<int> flat_table) {
    return Map(source, reps, local_size, std::move(flat_table), MapKind::Table, 0, 0, 1, 0, {});
}

bool Map::valid() const {
    if (!source_.valid() || reps_ <= 0 || local_size_ < 0 || kind_ == MapKind::Invalid) {
        return false;
    }
    try {
        validate_descriptor_indices(source_, reps_, local_size_, indices_, kind_, base_, rep_stride_, component_stride_, shift_, offsets_);
    } catch (const std::exception &) {
        return false;
    }
    return true;
}

int Map::reps() const {
    return reps_;
}

int Map::local_size() const {
    return local_size_;
}

MapKind Map::kind() const {
    return kind_;
}

int Map::base() const {
    return base_;
}

int Map::rep_stride() const {
    return rep_stride_;
}

int Map::component_stride() const {
    return component_stride_;
}

int Map::shift() const {
    return shift_;
}

const std::vector<int> &Map::offsets() const {
    return offsets_;
}

bool Map::stores_expanded_indices() const {
    return kind_ == MapKind::Single || kind_ == MapKind::ExplicitIndices || kind_ == MapKind::Table;
}

std::vector<int> Map::indices() const {
    if (!valid()) {
        throw std::runtime_error("cannot materialize invalid map indices");
    }
    return materialize_indices_from_descriptor(kind_, reps_, local_size_, indices_, base_, rep_stride_, component_stride_, shift_, offsets_);
}

int Map::index(int rep, int local_component) const {
    if (!valid()) {
        throw std::runtime_error("cannot read invalid map index");
    }
    if (rep < 0 || rep >= reps_ || local_component < 0 || local_component >= local_size_) {
        throw std::out_of_range("map index coordinates are out of range");
    }
    return map_index_from_descriptor(kind_, reps_, local_size_, indices_, base_, rep_stride_, component_stride_, shift_, offsets_, rep, local_component);
}

const Vec &Map::source() const {
    if (!source_.valid()) {
        throw std::runtime_error("invalid map has no source vector");
    }
    return source_;
}

Binding::Binding(Vec local_input, Map map)
    : local_input_(std::move(local_input)),
      map_(std::move(map)) {}

bool Binding::valid() const {
    return local_input_.valid() && map_.valid();
}

const Vec &Binding::local_input() const {
    if (!local_input_.valid()) {
        throw std::runtime_error("invalid binding has no local input");
    }
    return local_input_;
}

const Map &Binding::map() const {
    if (!map_.valid()) {
        throw std::runtime_error("invalid binding has no map");
    }
    return map_;
}

Binding bind(const Vec &local_input, const Map &map) {
    if (!local_input.valid()) {
        throw std::runtime_error("binding local input must be a valid vector");
    }
    if (!map.valid()) {
        throw std::runtime_error("binding map must be valid");
    }
    return Binding(local_input, map);
}

Vec map(const Function &function, std::vector<Binding> bindings, OutputMode mode) {
    switch (mode) {
        case OutputMode::Concat:
            return map(function, std::move(bindings), MappedOutput::concat());
        case OutputMode::Sum:
            return map(function, std::move(bindings), MappedOutput::sum());
        case OutputMode::Scatter:
            throw std::runtime_error("mapped scatter output requires MappedOutput::scatter(indices, output_size)");
        case OutputMode::WeightedSum:
            throw std::runtime_error("mapped weighted sum output requires MappedOutput::weighted_sum(weights)");
    }

    throw std::runtime_error("unsupported mapped output mode");
}

Vec map(const Function &function, std::vector<Binding> bindings, const MappedOutput &output) {
    const detail::FunctionCore &core = require_function_core(function);
    if (static_cast<int>(bindings.size()) != core.info.input_count) {
        throw std::runtime_error("mapped function must bind every callee input exactly once");
    }

    reject_bindings_to_non_inputs(core, bindings);

    int reps = -1;
    std::vector<detail::MappedBindingNode> ordered_bindings;
    ordered_bindings.reserve(bindings.size());

    for (std::size_t i = 0; i < core.inputs.size(); ++i) {
        const Vec &input = core.inputs[i];
        const FunctionInputInfo &input_info = core.info.inputs[i];
        const Binding *binding = find_binding_for_input(bindings, input.node_id());
        if (binding == nullptr) {
            throw std::runtime_error("mapped function is missing a binding for a callee input");
        }
        if (binding->map().local_size() != input_info.size) {
            throw std::runtime_error("mapped function binding local size does not match callee input size");
        }
        if (reps < 0) {
            reps = binding->map().reps();
        } else if (binding->map().reps() != reps) {
            throw std::runtime_error("mapped function bindings must all use the same repetition count");
        }

        detail::MappedBindingNode mapped;
        mapped.local_input = detail::vec_node(input);
        mapped.source = detail::vec_node(binding->map().source());
        mapped.reps = binding->map().reps();
        mapped.local_size = binding->map().local_size();
        mapped.indices = binding->map().indices_;
        mapped.map_kind = binding->map().kind_;
        mapped.base = binding->map().base_;
        mapped.rep_stride = binding->map().rep_stride_;
        mapped.component_stride = binding->map().component_stride_;
        mapped.shift = binding->map().shift_;
        mapped.offsets = binding->map().offsets_;
        ordered_bindings.push_back(std::move(mapped));
    }

    if (reps <= 0) {
        throw std::runtime_error("mapped function repetition count must be positive");
    }

    return detail::mapped_call_from_bindings(detail::function_core(function), std::move(ordered_bindings), reps, output);
}

namespace detail {

int mapped_index(const MappedBindingNode &binding, int reps, int rep, int local_component) {
    return map_index_from_descriptor(binding.map_kind,
                                     reps,
                                     binding.local_size,
                                     binding.indices,
                                     binding.base,
                                     binding.rep_stride,
                                     binding.component_stride,
                                     binding.shift,
                                     binding.offsets,
                                     rep,
                                     local_component);
}

std::vector<int> materialize_indices(const MappedBindingNode &binding, int reps) {
    return materialize_indices_from_descriptor(binding.map_kind,
                                               reps,
                                               binding.local_size,
                                               binding.indices,
                                               binding.base,
                                               binding.rep_stride,
                                               binding.component_stride,
                                               binding.shift,
                                               binding.offsets);
}

Vec mapped_call_from_bindings(std::shared_ptr<const FunctionCore> function,
                              std::vector<MappedBindingNode> bindings,
                              int reps,
                              MappedOutput output) {
    if (!function) {
        throw std::runtime_error("cannot create mapped call for invalid function");
    }
    require_positive_reps(reps);
    if (static_cast<int>(bindings.size()) != function->info.input_count) {
        throw std::runtime_error("mapped call must bind every function input exactly once");
    }

    for (std::size_t i = 0; i < bindings.size(); ++i) {
        MappedBindingNode &binding = bindings[i];
        if (!binding.local_input) {
            throw std::runtime_error("mapped call binding is missing local input");
        }
        if (!binding.source) {
            throw std::runtime_error("mapped call binding is missing source");
        }
        binding.source = vec_node(simplify_vec(vec_from_node(binding.source)));
        const Vec &expected_input = function->inputs[i];
        if (binding.local_input->id != expected_input.node_id()) {
            throw std::runtime_error("mapped call binding local input does not match function input order");
        }
        if (binding.local_size != function->info.inputs[i].size) {
            throw std::runtime_error("mapped call binding local size does not match function input size");
        }
        if (binding.reps != reps) {
            throw std::runtime_error("mapped call binding repetition count does not match mapped call");
        }
        (void)table_size(reps, binding.local_size);
        for (int rep = 0; rep < reps; ++rep) {
            for (int component = 0; component < binding.local_size; ++component) {
                const int index = mapped_index(binding, reps, rep, component);
                if (index < 0 || index >= binding.source->size) {
                throw std::runtime_error("mapped call source index is out of range");
                }
            }
        }
    }

    const int output_size = mapped_output_size(output, reps, function->info.output_size);

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::MappedFunctionCall;
    node->size = output_size;
    node->function = std::move(function);
    node->reps = reps;
    node->mapped_bindings = std::move(bindings);
    node->mapped_output = std::move(output);
    return vec_from_node(node);
}

Vec map_accum_call_from_bindings(std::shared_ptr<const FunctionCore> function,
                                 int carry_input_index,
                                 Vec initial_carry,
                                 std::vector<MappedBindingNode> bindings,
                                 int reps) {
    if (!function) {
        throw std::runtime_error("cannot create map-accum call for invalid function");
    }
    require_positive_reps(reps);
    if (!initial_carry.valid()) {
        throw std::runtime_error("map-accum initial carry must be a valid vector");
    }
    initial_carry = simplify_vec(initial_carry);
    if (carry_input_index < 0 || carry_input_index >= function->info.input_count) {
        throw std::runtime_error("map-accum carry input is not a function input");
    }
    if (function->info.output_count < 1) {
        throw std::runtime_error("map-accum step function must have at least one output");
    }

    const int carry_size = function->info.inputs[static_cast<std::size_t>(carry_input_index)].size;
    if (carry_size <= 0) {
        throw std::runtime_error("map-accum carry input must be non-empty");
    }
    if (initial_carry.size() != carry_size) {
        throw std::runtime_error("map-accum initial carry size does not match carry input size");
    }
    if (function->info.outputs[0].size != carry_size) {
        throw std::runtime_error("map-accum first step output must match carry input size");
    }
    if (static_cast<int>(bindings.size()) != function->info.input_count - 1) {
        throw std::runtime_error("map-accum must bind every non-carry step input exactly once");
    }

    std::size_t binding_index = 0;
    for (int input_index = 0; input_index < function->info.input_count; ++input_index) {
        if (input_index == carry_input_index) {
            continue;
        }
        if (binding_index >= bindings.size()) {
            throw std::runtime_error("map-accum binding count is inconsistent");
        }
        MappedBindingNode &binding = bindings[binding_index++];
        if (!binding.local_input) {
            throw std::runtime_error("map-accum binding is missing local input");
        }
        if (!binding.source) {
            throw std::runtime_error("map-accum binding is missing source");
        }
        binding.source = vec_node(simplify_vec(vec_from_node(binding.source)));
        const Vec &expected_input = function->inputs[static_cast<std::size_t>(input_index)];
        if (binding.local_input->id != expected_input.node_id()) {
            throw std::runtime_error("map-accum binding local input does not match function input order");
        }
        if (binding.local_size != function->info.inputs[static_cast<std::size_t>(input_index)].size) {
            throw std::runtime_error("map-accum binding local size does not match function input size");
        }
        if (binding.reps != reps) {
            throw std::runtime_error("map-accum binding repetition count does not match map-accum call");
        }
        (void)table_size(reps, binding.local_size);
        for (int rep = 0; rep < reps; ++rep) {
            for (int component = 0; component < binding.local_size; ++component) {
                const int index = mapped_index(binding, reps, rep, component);
                if (index < 0 || index >= binding.source->size) {
                    throw std::runtime_error("map-accum source index is out of range");
                }
            }
        }
    }

    const int output_size = checked_product(reps, function->info.output_size, "map-accum output size overflow");
    const int full_size = checked_sum(carry_size, output_size, "map-accum output size overflow");

    auto node = std::make_shared<VecNode>();
    node->kind = GraphNodeKind::MapAccumCall;
    node->size = full_size;
    node->lhs = vec_node(initial_carry);
    node->function = std::move(function);
    node->reps = reps;
    node->carry_input_index = carry_input_index;
    node->mapped_bindings = std::move(bindings);
    return vec_from_node(node);
}

} // namespace detail

} // namespace ad
