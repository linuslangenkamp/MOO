// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_CODEGEN_OPTIONS_H
#define MOO_AD_CODEGEN_OPTIONS_H

#include <string>

namespace ad {

struct CodegenOptions {
    std::string function_name = "ad_function";
    std::string scalar_type = "double";
    bool emit_sparsity_metadata = true;
};

} // namespace ad

#endif // MOO_AD_CODEGEN_OPTIONS_H
