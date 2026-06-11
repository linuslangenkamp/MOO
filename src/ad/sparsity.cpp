// SPDX-License-Identifier: LGPL-3.0-or-later
#include "sparsity.h"

#include "expr.h"
#include "function.h"
#include "vec.h"
#include "detail/function_core.h"
#include "detail/matrix_ops.h"
#include "detail/mapping.h"
#include "detail/node.h"

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

namespace ad {
namespace {

struct SparsityContext {
    Vars wrt;
    std::map<VarId, int> column;
};

SparsityContext make_context(const Vars &wrt) {
    SparsityContext context;
    context.wrt = wrt;

    std::set<VarId> seen;
    for (int i = 0; i < wrt.size(); ++i) {
        const Var &var = wrt[i];
        if (!var.valid()) {
            throw std::runtime_error("sparsity variable layout contains invalid variable");
        }
        if (!seen.insert(var.id()).second) {
            throw std::runtime_error("sparsity variable layout contains duplicate variable IDs");
        }
        context.column.emplace(var.id(), i);
    }

    return context;
}

class PatternBuilder {
public:
    PatternBuilder(int rows, int cols)
        : rows_(rows),
          cols_(cols) {
        if (rows < 0 || cols < 0) {
            throw std::runtime_error("sparsity dimensions must be non-negative");
        }
    }

    void add(int row, int col) {
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::runtime_error("sparsity entry out of bounds");
        }
        entries_.emplace_back(row, col);
    }

    void add_pattern(const SparsityPattern &pattern, int row_offset = 0) {
        if (pattern.cols() != cols_) {
            throw std::runtime_error("cannot combine sparsity patterns with different column counts");
        }
        for (const auto &[row, col] : pattern.entries()) {
            add(row + row_offset, col);
        }
    }

    SparsityPattern build() {
        return SparsityPattern(rows_, cols_, std::move(entries_));
    }

private:
    int rows_ = 0;
    int cols_ = 0;
    std::vector<std::pair<int, int>> entries_;
};

bool is_zero(const std::shared_ptr<const detail::ScalarNode> &node) {
    return node && node->kind == GraphNodeKind::ScalarConstant && node->value == 0.0;
}

bool is_zero(const std::shared_ptr<const detail::VecNode> &node) {
    if (!node || node->kind != GraphNodeKind::VectorConstant) {
        return false;
    }
    for (double value : node->constants) {
        if (value != 0.0) {
            return false;
        }
    }
    return true;
}

bool literal_zero_entry(const std::shared_ptr<const detail::VecNode> &node, int index) {
    if (!node || node->kind != GraphNodeKind::VectorConstant) {
        return false;
    }
    if (index < 0 || index >= node->size) {
        throw std::runtime_error("literal zero entry index out of range");
    }
    return node->constants[static_cast<std::size_t>(index)] == 0.0;
}

SparsityPattern scalar_sparsity(const std::shared_ptr<const detail::ScalarNode> &node, const SparsityContext &context);
SparsityPattern vec_sparsity(const std::shared_ptr<const detail::VecNode> &node, const SparsityContext &context);

SparsityPattern empty_pattern(int rows, int cols) {
    return SparsityPattern(rows, cols, {});
}

SparsityPattern union_patterns(int rows, int cols, const std::vector<SparsityPattern> &patterns) {
    PatternBuilder builder(rows, cols);
    for (const SparsityPattern &pattern : patterns) {
        builder.add_pattern(pattern);
    }
    return builder.build();
}

SparsityPattern scalar_to_row(const SparsityPattern &scalar, int row, int rows) {
    PatternBuilder builder(rows, scalar.cols());
    for (const auto &[unused_row, col] : scalar.entries()) {
        (void)unused_row;
        builder.add(row, col);
    }
    return builder.build();
}

SparsityPattern select_row(const SparsityPattern &pattern, int selected_row) {
    PatternBuilder builder(1, pattern.cols());
    for (const auto &[row, col] : pattern.entries()) {
        if (row == selected_row) {
            builder.add(0, col);
        }
    }
    return builder.build();
}

SparsityPattern collapse_to_scalar(const SparsityPattern &pattern) {
    PatternBuilder builder(1, pattern.cols());
    for (const auto &[row, col] : pattern.entries()) {
        (void)row;
        builder.add(0, col);
    }
    return builder.build();
}

SparsityPattern repeat_scalar_by_rows(const SparsityPattern &scalar, int rows) {
    PatternBuilder builder(rows, scalar.cols());
    for (int row = 0; row < rows; ++row) {
        for (const auto &[unused_row, col] : scalar.entries()) {
            (void)unused_row;
            builder.add(row, col);
        }
    }
    return builder.build();
}

std::vector<std::vector<int>> row_support(const SparsityPattern &pattern) {
    std::vector<std::vector<int>> rows(static_cast<std::size_t>(pattern.rows()));
    for (const auto &[row, col] : pattern.entries()) {
        rows[static_cast<std::size_t>(row)].push_back(col);
    }
    return rows;
}

struct InputColumn {
    int group = -1;
    int component = -1;
};

InputColumn input_column(const detail::FunctionCore &function, int col) {
    int offset = 0;
    for (int group = 0; group < function.info.input_count; ++group) {
        const int size = function.info.inputs[static_cast<std::size_t>(group)].size;
        if (col >= offset && col < offset + size) {
            return {group, col - offset};
        }
        offset += size;
    }
    throw std::runtime_error("callee sparsity contains invalid local input column");
}

SparsityPattern local_jacobian(const detail::FunctionCore &function) {
    return sparsity(function.output, function.input_vars);
}

void add_argument_row(PatternBuilder &builder, const std::vector<std::vector<int>> &arg_rows, int output_row, int arg_row) {
    for (int col : arg_rows[static_cast<std::size_t>(arg_row)]) {
        builder.add(output_row, col);
    }
}

SparsityPattern function_call_sparsity(const std::shared_ptr<const detail::VecNode> &node, const SparsityContext &context) {
    if (!node->function) {
        throw std::runtime_error("function call node missing callee while computing sparsity");
    }

    const SparsityPattern local = local_jacobian(*node->function);
    std::vector<SparsityPattern> arg_patterns;
    std::vector<std::vector<std::vector<int>>> arg_rows;
    arg_patterns.reserve(node->arguments.size());
    arg_rows.reserve(node->arguments.size());
    for (const auto &argument : node->arguments) {
        arg_patterns.push_back(vec_sparsity(argument, context));
        arg_rows.push_back(row_support(arg_patterns.back()));
    }

    PatternBuilder builder(node->size, context.wrt.size());
    for (const auto &[local_row, local_col] : local.entries()) {
        InputColumn input = input_column(*node->function, local_col);
        add_argument_row(builder, arg_rows[static_cast<std::size_t>(input.group)], local_row, input.component);
    }
    return builder.build();
}

SparsityPattern mapped_function_call_sparsity(const std::shared_ptr<const detail::VecNode> &node, const SparsityContext &context) {
    if (!node->function) {
        throw std::runtime_error("mapped function call node missing callee while computing sparsity");
    }
    if (static_cast<int>(node->mapped_bindings.size()) != node->function->info.input_count) {
        throw std::runtime_error("mapped function call binding count does not match callee input count while computing sparsity");
    }

    const SparsityPattern local = local_jacobian(*node->function);
    std::vector<SparsityPattern> source_patterns;
    std::vector<std::vector<std::vector<int>>> source_rows;
    source_patterns.reserve(node->mapped_bindings.size());
    source_rows.reserve(node->mapped_bindings.size());

    for (std::size_t group = 0; group < node->mapped_bindings.size(); ++group) {
        const detail::MappedBindingNode &binding = node->mapped_bindings[group];
        if (!binding.source) {
            throw std::runtime_error("mapped function call binding is missing source while computing sparsity");
        }
        if (binding.local_size != node->function->info.inputs[group].size) {
            throw std::runtime_error("mapped function call binding local size does not match callee input while computing sparsity");
        }
        source_patterns.push_back(vec_sparsity(binding.source, context));
        source_rows.push_back(row_support(source_patterns.back()));
    }

    PatternBuilder builder(node->size, context.wrt.size());
    const int local_output_size = node->function->info.output_size;
    for (const auto &[local_row, local_col] : local.entries()) {
        InputColumn input = input_column(*node->function, local_col);
        const detail::MappedBindingNode &binding = node->mapped_bindings[static_cast<std::size_t>(input.group)];
        for (int rep = 0; rep < node->reps; ++rep) {
            const int source_row = detail::mapped_index(binding, node->reps, rep, input.component);
            if (source_row < 0 || source_row >= static_cast<int>(source_rows[static_cast<std::size_t>(input.group)].size())) {
                throw std::runtime_error("mapped function call source index is out of range while computing sparsity");
            }
            int global_output_row = -1;
            if (node->mapped_output.mode() == OutputMode::Concat) {
                global_output_row = rep * local_output_size + local_row;
            } else if (node->mapped_output.mode() == OutputMode::Scatter) {
                const int output_index = rep * local_output_size + local_row;
                if (output_index < 0 || output_index >= static_cast<int>(node->mapped_output.indices().size())) {
                    throw std::runtime_error("mapped scatter output index table is invalid while computing sparsity");
                }
                global_output_row = node->mapped_output.indices()[static_cast<std::size_t>(output_index)];
            } else if (node->mapped_output.mode() == OutputMode::Sum) {
                global_output_row = local_row;
            } else if (node->mapped_output.mode() == OutputMode::WeightedSum) {
                if (node->mapped_output.weights().size() != static_cast<std::size_t>(node->reps)) {
                    throw std::runtime_error("mapped weighted sum weights are invalid while computing sparsity");
                }
                if (node->mapped_output.weights()[static_cast<std::size_t>(rep)] == 0.0) {
                    continue;
                }
                global_output_row = local_row;
            }
            for (int col : source_rows[static_cast<std::size_t>(input.group)][static_cast<std::size_t>(source_row)]) {
                builder.add(global_output_row, col);
            }
        }
    }
    return builder.build();
}

std::vector<std::set<int>> row_dependency_sets(const SparsityPattern &pattern, int rows) {
    std::vector<std::set<int>> deps(static_cast<std::size_t>(rows));
    for (const auto &[row, col] : pattern.entries()) {
        deps[static_cast<std::size_t>(row)].insert(col);
    }
    return deps;
}

void add_dependencies(PatternBuilder &builder, int row, const std::set<int> &deps) {
    for (int col : deps) {
        builder.add(row, col);
    }
}

int input_offset(const detail::FunctionCore &core, int input_index) {
    int offset = 0;
    for (int i = 0; i < input_index; ++i) {
        offset += core.info.inputs[static_cast<std::size_t>(i)].size;
    }
    return offset;
}

int map_accum_emit_offset(const detail::FunctionCore &core, int reps, int output_index) {
    const int carry_size = core.info.outputs[0].size;
    int offset = carry_size + reps * carry_size;
    for (int i = 1; i < output_index; ++i) {
        offset += reps * core.info.outputs[static_cast<std::size_t>(i)].size;
    }
    return offset;
}

SparsityPattern map_accum_call_sparsity(const std::shared_ptr<const detail::VecNode> &node, const SparsityContext &context) {
    if (!node->function) {
        throw std::runtime_error("map-accum call node is missing callee while computing sparsity");
    }
    if (node->carry_input_index < 0 || node->carry_input_index >= node->function->info.input_count) {
        throw std::runtime_error("map-accum call has invalid carry input index while computing sparsity");
    }

    const int carry_size = node->function->info.outputs[0].size;
    PatternBuilder builder(node->size, context.wrt.size());
    std::vector<std::set<int>> carry_deps = row_dependency_sets(vec_sparsity(node->lhs, context), carry_size);
    std::vector<std::vector<std::set<int>>> source_deps;
    source_deps.reserve(node->mapped_bindings.size());
    for (const detail::MappedBindingNode &binding : node->mapped_bindings) {
        source_deps.push_back(row_dependency_sets(vec_sparsity(binding.source, context), binding.source->size));
    }

    const SparsityPattern local = detail::function_from_core(node->function).jacobian_sparsity();
    for (int rep = 0; rep < node->reps; ++rep) {
        std::vector<std::set<int>> input_deps(static_cast<std::size_t>(node->function->info.input_size));
        std::size_t sequence_binding = 0;
        for (int input_index = 0; input_index < node->function->info.input_count; ++input_index) {
            const int offset = input_offset(*node->function, input_index);
            const int size = node->function->info.inputs[static_cast<std::size_t>(input_index)].size;
            if (input_index == node->carry_input_index) {
                for (int i = 0; i < size; ++i) {
                    input_deps[static_cast<std::size_t>(offset + i)] = carry_deps[static_cast<std::size_t>(i)];
                }
            } else {
                const detail::MappedBindingNode &binding = node->mapped_bindings[sequence_binding];
                for (int component = 0; component < size; ++component) {
                    const int source_row = detail::mapped_index(binding, node->reps, rep, component);
                    input_deps[static_cast<std::size_t>(offset + component)] = source_deps[sequence_binding][static_cast<std::size_t>(source_row)];
                }
                ++sequence_binding;
            }
        }

        std::vector<std::set<int>> output_deps(static_cast<std::size_t>(node->function->info.output_size));
        for (const auto &[row, col] : local.entries()) {
            output_deps[static_cast<std::size_t>(row)].insert(input_deps[static_cast<std::size_t>(col)].begin(),
                                                              input_deps[static_cast<std::size_t>(col)].end());
        }

        const FunctionOutputInfo &carry_output = node->function->info.outputs[0];
        std::vector<std::set<int>> next_carry_deps(static_cast<std::size_t>(carry_size));
        for (int i = 0; i < carry_size; ++i) {
            next_carry_deps[static_cast<std::size_t>(i)] = output_deps[static_cast<std::size_t>(carry_output.offset + i)];
            add_dependencies(builder, carry_size + rep * carry_size + i, next_carry_deps[static_cast<std::size_t>(i)]);
        }
        for (int output_index = 1; output_index < node->function->info.output_count; ++output_index) {
            const FunctionOutputInfo &output = node->function->info.outputs[static_cast<std::size_t>(output_index)];
            const int target_offset = map_accum_emit_offset(*node->function, node->reps, output_index) + rep * output.size;
            for (int i = 0; i < output.size; ++i) {
                add_dependencies(builder, target_offset + i, output_deps[static_cast<std::size_t>(output.offset + i)]);
            }
        }
        carry_deps = std::move(next_carry_deps);
    }

    for (int i = 0; i < carry_size; ++i) {
        add_dependencies(builder, i, carry_deps[static_cast<std::size_t>(i)]);
    }
    return builder.build();
}

SparsityPattern scalar_sparsity(const std::shared_ptr<const detail::ScalarNode> &node, const SparsityContext &context) {
    if (!node) {
        throw std::runtime_error("invalid scalar graph node while computing sparsity");
    }

    switch (node->kind) {
        case GraphNodeKind::ScalarConstant:
        case GraphNodeKind::ScalarParameter:
            return empty_pattern(1, context.wrt.size());
        case GraphNodeKind::ScalarVariable: {
            PatternBuilder builder(1, context.wrt.size());
            const auto found = context.column.find(node->var.id());
            if (found != context.column.end()) {
                builder.add(0, found->second);
            }
            return builder.build();
        }
        case GraphNodeKind::ScalarUnary:
            return scalar_sparsity(node->lhs, context);
        case GraphNodeKind::ScalarBinary:
            if (node->binary == detail::ScalarBinaryOp::Mul && (is_zero(node->lhs) || is_zero(node->rhs))) {
                return empty_pattern(1, context.wrt.size());
            }
            return union_patterns(1, context.wrt.size(), {scalar_sparsity(node->lhs, context), scalar_sparsity(node->rhs, context)});
        case GraphNodeKind::VectorElement:
            return select_row(vec_sparsity(node->vec, context), node->index);
        case GraphNodeKind::Sum:
            return collapse_to_scalar(vec_sparsity(node->vec, context));
        case GraphNodeKind::Dot:
            return union_patterns(1, context.wrt.size(), {
                collapse_to_scalar(vec_sparsity(node->vec_lhs, context)),
                collapse_to_scalar(vec_sparsity(node->vec_rhs, context)),
            });
        default:
            throw std::runtime_error("unsupported scalar graph node while computing sparsity");
    }
}

SparsityPattern vector_variable_sparsity(const std::shared_ptr<const detail::VecNode> &node, const SparsityContext &context) {
    PatternBuilder builder(node->size, context.wrt.size());
    for (int row = 0; row < node->size; ++row) {
        const auto found = context.column.find(node->vars[static_cast<std::size_t>(row)].id());
        if (found != context.column.end()) {
            builder.add(row, found->second);
        }
    }
    return builder.build();
}

SparsityPattern vec_sparsity(const std::shared_ptr<const detail::VecNode> &node, const SparsityContext &context) {
    if (!node) {
        throw std::runtime_error("invalid vector graph node while computing sparsity");
    }

    switch (node->kind) {
        case GraphNodeKind::VectorVariable:
            return vector_variable_sparsity(node, context);
        case GraphNodeKind::VectorParameter:
        case GraphNodeKind::VectorConstant:
            return empty_pattern(node->size, context.wrt.size());
        case GraphNodeKind::VectorFromElements: {
            PatternBuilder builder(node->size, context.wrt.size());
            for (int row = 0; row < node->size; ++row) {
                builder.add_pattern(scalar_to_row(scalar_sparsity(node->elements[static_cast<std::size_t>(row)], context), row, node->size));
            }
            return builder.build();
        }
        case GraphNodeKind::VectorUnary:
            return vec_sparsity(node->lhs, context);
        case GraphNodeKind::VectorBinary:
            if (node->binary == detail::VecBinaryOp::Mul && (is_zero(node->lhs) || is_zero(node->rhs))) {
                return empty_pattern(node->size, context.wrt.size());
            }
            return union_patterns(node->size, context.wrt.size(), {vec_sparsity(node->lhs, context), vec_sparsity(node->rhs, context)});
        case GraphNodeKind::VectorScale: {
            if (is_zero(node->scale) || is_zero(node->lhs)) {
                return empty_pattern(node->size, context.wrt.size());
            }
            return union_patterns(node->size, context.wrt.size(), {
                vec_sparsity(node->lhs, context),
                repeat_scalar_by_rows(scalar_sparsity(node->scale, context), node->size),
            });
        }
        case GraphNodeKind::DenseMatVec: {
            const SparsityPattern child = vec_sparsity(node->lhs, context);
            const auto child_rows = row_support(child);
            PatternBuilder builder(node->size, context.wrt.size());
            for (int row = 0; row < node->dense.rows; ++row) {
                for (int inner = 0; inner < node->dense.cols; ++inner) {
                    if (node->dense(row, inner) == 0.0) {
                        continue;
                    }
                    for (int col : child_rows[static_cast<std::size_t>(inner)]) {
                        builder.add(row, col);
                    }
                }
            }
            return builder.build();
        }
        case GraphNodeKind::SparseMatVec: {
            const SparsityPattern child = vec_sparsity(node->lhs, context);
            const auto child_rows = row_support(child);
            PatternBuilder builder(node->size, context.wrt.size());
            for (int k = 0; k < node->sparse.nnz(); ++k) {
                if (node->sparse.values[static_cast<std::size_t>(k)] == 0.0) {
                    continue;
                }
                const int row = node->sparse.row[static_cast<std::size_t>(k)];
                const int inner = node->sparse.col[static_cast<std::size_t>(k)];
                for (int col : child_rows[static_cast<std::size_t>(inner)]) {
                    builder.add(row, col);
                }
            }
            return builder.build();
        }
        case GraphNodeKind::SymbolicMatVec: {
            const SparsityPattern matrix = vec_sparsity(node->lhs, context);
            const SparsityPattern rhs = vec_sparsity(node->rhs, context);
            const auto matrix_rows = row_support(matrix);
            const auto rhs_rows = row_support(rhs);
            PatternBuilder builder(node->size, context.wrt.size());
            for (int row = 0; row < node->mat_lhs_rows; ++row) {
                for (int inner = 0; inner < node->mat_lhs_cols; ++inner) {
                    const int matrix_index = detail::matrix_flat_index(row, inner, node->mat_lhs_rows, node->mat_lhs_cols, node->mat_lhs_layout);
                    if (!literal_zero_entry(node->rhs, inner)) {
                        for (int col : matrix_rows[static_cast<std::size_t>(matrix_index)]) {
                            builder.add(row, col);
                        }
                    }
                    if (!literal_zero_entry(node->lhs, matrix_index)) {
                        for (int col : rhs_rows[static_cast<std::size_t>(inner)]) {
                            builder.add(row, col);
                        }
                    }
                }
            }
            return builder.build();
        }
        case GraphNodeKind::SymbolicMatMul: {
            const SparsityPattern lhs = vec_sparsity(node->lhs, context);
            const SparsityPattern rhs = vec_sparsity(node->rhs, context);
            const auto lhs_rows = row_support(lhs);
            const auto rhs_rows = row_support(rhs);
            PatternBuilder builder(node->size, context.wrt.size());
            for (int row = 0; row < node->mat_lhs_rows; ++row) {
                for (int col_index = 0; col_index < node->mat_rhs_cols; ++col_index) {
                    const int output_row = detail::matrix_flat_index(row, col_index, node->mat_lhs_rows, node->mat_rhs_cols, node->mat_result_layout);
                    for (int inner = 0; inner < node->mat_lhs_cols; ++inner) {
                        const int lhs_index = detail::matrix_flat_index(row, inner, node->mat_lhs_rows, node->mat_lhs_cols, node->mat_lhs_layout);
                        const int rhs_index = detail::matrix_flat_index(inner, col_index, node->mat_rhs_rows, node->mat_rhs_cols, node->mat_rhs_layout);
                        if (!literal_zero_entry(node->rhs, rhs_index)) {
                            for (int pattern_col : lhs_rows[static_cast<std::size_t>(lhs_index)]) {
                                builder.add(output_row, pattern_col);
                            }
                        }
                        if (!literal_zero_entry(node->lhs, lhs_index)) {
                            for (int pattern_col : rhs_rows[static_cast<std::size_t>(rhs_index)]) {
                                builder.add(output_row, pattern_col);
                            }
                        }
                    }
                }
            }
            return builder.build();
        }
        case GraphNodeKind::OuterProduct: {
            const SparsityPattern lhs = vec_sparsity(node->lhs, context);
            const SparsityPattern rhs = vec_sparsity(node->rhs, context);
            const auto lhs_rows = row_support(lhs);
            const auto rhs_rows = row_support(rhs);
            PatternBuilder builder(node->size, context.wrt.size());
            for (int row = 0; row < node->mat_lhs_rows; ++row) {
                for (int col_index = 0; col_index < node->mat_rhs_cols; ++col_index) {
                    const int output_row = detail::matrix_flat_index(row, col_index, node->mat_lhs_rows, node->mat_rhs_cols, node->mat_result_layout);
                    if (!literal_zero_entry(node->rhs, col_index)) {
                        for (int pattern_col : lhs_rows[static_cast<std::size_t>(row)]) {
                            builder.add(output_row, pattern_col);
                        }
                    }
                    if (!literal_zero_entry(node->lhs, row)) {
                        for (int pattern_col : rhs_rows[static_cast<std::size_t>(col_index)]) {
                            builder.add(output_row, pattern_col);
                        }
                    }
                }
            }
            return builder.build();
        }
        case GraphNodeKind::LinearSolve: {
            const SparsityPattern matrix = vec_sparsity(node->lhs, context);
            const SparsityPattern rhs = vec_sparsity(node->rhs, context);
            PatternBuilder builder(node->size, context.wrt.size());
            for (const auto &[unused_row, col] : matrix.entries()) {
                (void)unused_row;
                for (int row = 0; row < node->size; ++row) {
                    builder.add(row, col);
                }
            }
            for (const auto &[unused_row, col] : rhs.entries()) {
                (void)unused_row;
                for (int row = 0; row < node->size; ++row) {
                    builder.add(row, col);
                }
            }
            return builder.build();
        }
        case GraphNodeKind::Slice: {
            const SparsityPattern child = vec_sparsity(node->lhs, context);
            PatternBuilder builder(node->size, context.wrt.size());
            for (const auto &[row, col] : child.entries()) {
                if (row >= node->start && row < node->start + node->size) {
                    builder.add(row - node->start, col);
                }
            }
            return builder.build();
        }
        case GraphNodeKind::ScatterSlice: {
            const SparsityPattern child = vec_sparsity(node->lhs, context);
            PatternBuilder builder(node->size, context.wrt.size());
            for (const auto &[row, col] : child.entries()) {
                builder.add(row + node->start, col);
            }
            return builder.build();
        }
        case GraphNodeKind::Gather: {
            const SparsityPattern child = vec_sparsity(node->lhs, context);
            const auto child_rows = row_support(child);
            PatternBuilder builder(node->size, context.wrt.size());
            for (int row = 0; row < node->size; ++row) {
                const int source_row = node->indices[static_cast<std::size_t>(row)];
                for (int col : child_rows[static_cast<std::size_t>(source_row)]) {
                    builder.add(row, col);
                }
            }
            return builder.build();
        }
        case GraphNodeKind::ScatterAdd: {
            const SparsityPattern child = vec_sparsity(node->lhs, context);
            PatternBuilder builder(node->size, context.wrt.size());
            for (const auto &[row, col] : child.entries()) {
                builder.add(node->indices[static_cast<std::size_t>(row)], col);
            }
            return builder.build();
        }
        case GraphNodeKind::Concat: {
            PatternBuilder builder(node->size, context.wrt.size());
            builder.add_pattern(vec_sparsity(node->lhs, context));
            builder.add_pattern(vec_sparsity(node->rhs, context), node->lhs->size);
            return builder.build();
        }
        case GraphNodeKind::FunctionCall:
            return function_call_sparsity(node, context);
        case GraphNodeKind::MappedFunctionCall:
            return mapped_function_call_sparsity(node, context);
        case GraphNodeKind::MapAccumCall:
            return map_accum_call_sparsity(node, context);
        default:
            throw std::runtime_error("unsupported vector graph node while computing sparsity");
    }
}

} // namespace

SparsityPattern::SparsityPattern(int rows, int cols, std::vector<std::pair<int, int>> entries)
    : rows_(rows),
      cols_(cols),
      entries_(std::move(entries)) {
    if (rows_ < 0 || cols_ < 0) {
        throw std::runtime_error("sparsity dimensions must be non-negative");
    }
    for (const auto &[row, col] : entries_) {
        if (row < 0 || row >= rows_ || col < 0 || col >= cols_) {
            throw std::runtime_error("sparsity entry out of bounds");
        }
    }
    std::sort(entries_.begin(), entries_.end());
    entries_.erase(std::unique(entries_.begin(), entries_.end()), entries_.end());
}

int SparsityPattern::rows() const {
    return rows_;
}

int SparsityPattern::cols() const {
    return cols_;
}

const std::vector<std::pair<int, int>> &SparsityPattern::entries() const {
    return entries_;
}

bool SparsityPattern::contains(int row, int col) const {
    return std::binary_search(entries_.begin(), entries_.end(), std::make_pair(row, col));
}

bool SparsityPattern::empty() const {
    return entries_.empty();
}

int SparsityPattern::nnz() const {
    return static_cast<int>(entries_.size());
}

SparsityPattern sparsity(const Expr &expr, const Vars &wrt) {
    if (!expr.valid()) {
        throw std::runtime_error("cannot compute sparsity of invalid scalar expression");
    }
    return scalar_sparsity(detail::scalar_node(expr), make_context(wrt));
}

SparsityPattern sparsity(const Vec &vec, const Vars &wrt) {
    if (!vec.valid()) {
        throw std::runtime_error("cannot compute sparsity of invalid vector expression");
    }
    return vec_sparsity(detail::vec_node(vec), make_context(wrt));
}

} // namespace ad
