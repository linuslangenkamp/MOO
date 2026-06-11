// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_FUNCTION_H
#define MOO_AD_FUNCTION_H

#include "sparsity.h"
#include "vec.h"

#include <memory>
#include <string>
#include <vector>

namespace ad {

class Function;
class Values;
class EvalWorkspace;

namespace detail {
struct FunctionCore;
const std::shared_ptr<const FunctionCore> &function_core(const Function &function);
Function function_from_core(std::shared_ptr<const FunctionCore> core);
}

struct FunctionInputInfo {
    std::string label;
    int size = 0;
    NodeId node_id = 0;
};

struct FunctionOutputInfo {
    std::string label;
    int offset = 0;
    int size = 0;
    NodeId node_id = 0;
    GraphInfo graph;
};

struct FunctionInfo {
    int input_count = 0;
    int input_size = 0;
    int output_count = 0;
    int output_size = 0;
    int parameter_count = 0;
    std::vector<FunctionInputInfo> inputs;
    std::vector<FunctionOutputInfo> outputs;
    GraphInfo output_graph;
};

class Function {
public:
    Function() = default;
    Function(std::vector<Vec> inputs, Vec output);
    Function(std::vector<Vec> inputs, std::vector<Vec> outputs);
    Function(std::vector<Vec> inputs, Params parameters, Vec output);
    Function(std::vector<Vec> inputs, Params parameters, std::vector<Vec> outputs);

    bool valid() const;
    int input_count() const;
    int input_size() const;
    int output_count() const;
    int output_size() const;

    const std::vector<Vec> &inputs() const;
    const Vec &output() const;
    const std::vector<Vec> &outputs() const;
    const Vec &output(int index) const;
    int output_offset(int index) const;

    Vars input_vars() const;
    Params parameters() const;

    FunctionInfo info() const;
    Vec call(std::vector<Vec> arguments) const;
    std::vector<Vec> call_outputs(std::vector<Vec> arguments) const;
    Function forward_function() const;
    Function reverse_function() const;
    Vec forward(const Vec &seed) const;
    Vec reverse(const Vec &lambda) const;
    SparsityPattern jacobian_sparsity() const;
    void eval(const double *input, int input_size,
              const double *params, int param_size,
              EvalWorkspace &workspace,
              double *output, int output_size) const;
    void eval(const Values &values,
              EvalWorkspace &workspace,
              double *output, int output_size) const;

private:
    explicit Function(std::shared_ptr<const detail::FunctionCore> core);

    std::shared_ptr<const detail::FunctionCore> core_;

    friend const std::shared_ptr<const detail::FunctionCore> &detail::function_core(const Function &function);
    friend Function detail::function_from_core(std::shared_ptr<const detail::FunctionCore> core);
};

} // namespace ad

#endif // MOO_AD_FUNCTION_H
