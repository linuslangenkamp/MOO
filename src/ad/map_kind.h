// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_MAP_KIND_H
#define MOO_AD_MAP_KIND_H

#include <vector>

namespace ad {

enum class MapKind {
    Invalid,
    Single,
    Blocks,
    ShiftedBlocks,
    Stride,
    Stencil,
    ExplicitIndices,
    Table
};

enum class OutputMode {
    Concat,
    Scatter,
    Sum,
    WeightedSum
};

class MappedOutput {
public:
    MappedOutput() = default;

    static MappedOutput concat();
    static MappedOutput scatter(std::vector<int> indices, int output_size);
    static MappedOutput sum();
    static MappedOutput weighted_sum(std::vector<double> weights);

    OutputMode mode() const;
    int output_size() const;
    const std::vector<int> &indices() const;
    const std::vector<double> &weights() const;

private:
    OutputMode mode_ = OutputMode::Concat;
    int output_size_ = 0;
    std::vector<int> indices_;
    std::vector<double> weights_;
};

} // namespace ad

#endif // MOO_AD_MAP_KIND_H
