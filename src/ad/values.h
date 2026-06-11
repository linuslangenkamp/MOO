// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_VALUES_H
#define MOO_AD_VALUES_H

#include "symbol.h"

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

namespace ad {

class Expr;
class Function;
class Vec;

namespace detail {
struct EvalWorkspaceAccess;
} // namespace detail

class Values {
public:
    Values() = default;

    void clear();

    void set(const Var &var, double value);
    void set(const Param &param, double value);
    void set(const Expr &symbol, double value);
    void set(const Vec &symbol_group, const double *values, int size);
    void set(const Vec &symbol_group, const std::vector<double> &values);

    bool has(const Var &var) const;
    bool has(const Param &param) const;
    double get(const Var &var) const;
    double get(const Param &param) const;

private:
    std::unordered_map<VarId, double> vars_;
    std::unordered_map<ParamId, double> params_;
};

class EvalWorkspace {
public:
    EvalWorkspace() = default;

    void clear();
    void reserve(const Function &function);
    void reserve_doubles(std::size_t count);

private:
    struct VarFrame {
        const Vars *vars = nullptr;
        const double *values = nullptr;
        int size = 0;
    };

    struct ParamFrame {
        const Params *params = nullptr;
        const double *values = nullptr;
        int size = 0;
    };

    struct Buffer {
        std::unique_ptr<double[]> data;
        int capacity = 0;
    };

    struct PointerBuffer {
        std::unique_ptr<double *[]> data;
        int capacity = 0;
    };

    std::vector<Buffer> buffers_;
    std::vector<PointerBuffer> pointer_buffers_;
    std::size_t used_buffers_ = 0;
    std::size_t used_pointer_buffers_ = 0;
    std::vector<VarFrame> var_frames_;
    std::vector<ParamFrame> param_frames_;

    friend struct detail::EvalWorkspaceAccess;
};

} // namespace ad

#endif // MOO_AD_VALUES_H
