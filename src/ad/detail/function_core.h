// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_DETAIL_FUNCTION_CORE_H
#define MOO_AD_DETAIL_FUNCTION_CORE_H

#include "../function.h"

#include <mutex>

namespace ad::detail {

struct FunctionCore {
    std::vector<Vec> inputs;
    Vec output;
    Vars input_vars;
    Params parameters;
    FunctionInfo info;
    mutable std::mutex transform_mutex;
    mutable std::shared_ptr<const FunctionCore> forward_cache;
    mutable std::shared_ptr<const FunctionCore> reverse_cache;
};

} // namespace ad::detail

#endif // MOO_AD_DETAIL_FUNCTION_CORE_H
