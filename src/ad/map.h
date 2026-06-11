// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_MAP_H
#define MOO_AD_MAP_H

#include "function.h"
#include "map_kind.h"

#include <vector>

namespace ad {

class Binding;

class Map {
public:
    Map() = default;

    static Map single(const Vec &source, std::vector<int> indices);
    static Map blocks(const Vec &source, int reps, int block_size);
    static Map shifted_blocks(const Vec &source, int reps, int block_size, int shift);
    static Map stride(const Vec &source, int reps, int local_size, int base, int rep_stride, int component_stride);
    static Map stencil(const Vec &source, int reps, int base, int rep_stride, std::vector<int> offsets);
    static Map explicit_indices(const Vec &source, int reps, int local_size, std::vector<int> flat_indices);
    static Map table(const Vec &source, int reps, int local_size, std::vector<int> flat_table);

    bool valid() const;
    int reps() const;
    int local_size() const;
    MapKind kind() const;
    int base() const;
    int rep_stride() const;
    int component_stride() const;
    int shift() const;
    const std::vector<int> &offsets() const;
    bool stores_expanded_indices() const;
    std::vector<int> indices() const;
    int index(int rep, int local_component) const;
    const Vec &source() const;

private:
    Map(Vec source,
        int reps,
        int local_size,
        std::vector<int> indices,
        MapKind kind,
        int base,
        int rep_stride,
        int component_stride,
        int shift,
        std::vector<int> offsets);

    Vec source_;
    int reps_ = 0;
    int local_size_ = 0;
    std::vector<int> indices_;
    MapKind kind_ = MapKind::Invalid;
    int base_ = 0;
    int rep_stride_ = 0;
    int component_stride_ = 1;
    int shift_ = 0;
    std::vector<int> offsets_;

    friend Vec map(const Function &function, std::vector<Binding> bindings, const MappedOutput &output);
};

class Binding {
public:
    Binding() = default;

    bool valid() const;
    const Vec &local_input() const;
    const Map &map() const;

private:
    Binding(Vec local_input, Map map);

    Vec local_input_;
    Map map_;

    friend Binding bind(const Vec &local_input, const Map &map);
};

Binding bind(const Vec &local_input, const Map &map);

Vec map(const Function &function, std::vector<Binding> bindings, OutputMode mode = OutputMode::Concat);
Vec map(const Function &function, std::vector<Binding> bindings, const MappedOutput &output);

} // namespace ad

#endif // MOO_AD_MAP_H
