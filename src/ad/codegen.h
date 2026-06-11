// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_CODEGEN_H
#define MOO_AD_CODEGEN_H

#include "codegen_options.h"
#include "function.h"

#include <string>

namespace ad {

struct CodegenResult {
    std::string source;
    std::string entry_name;
    int input_size = 0;
    int parameter_size = 0;
    int output_size = 0;
};

CodegenResult generate_c(const Function &function,
                         const CodegenOptions &options = {});

} // namespace ad

#endif // MOO_AD_CODEGEN_H
