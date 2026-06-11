// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_DETAIL_TRAVERSAL_H
#define MOO_AD_DETAIL_TRAVERSAL_H

#include "../expr.h"
#include "../vec.h"

#include <string>
#include <vector>

namespace ad::detail {

struct VecVariableGroupInfo {
    std::string label;
    SymbolGroupId group_id = 0;
    NodeId node_id = 0;
    int size = 0;
    std::vector<Var> vars;
};

struct VecParameterGroupInfo {
    std::string label;
    SymbolGroupId group_id = 0;
    NodeId node_id = 0;
    int size = 0;
    std::vector<Param> params;
};

Vars collect_variables(const Expr &expr);
Vars collect_variables(const Vec &vec);
Params collect_parameters(const Expr &expr);
Params collect_parameters(const Vec &vec);

bool is_vector_variable_group(const Vec &vec);
bool is_vector_parameter_group(const Vec &vec);
VecVariableGroupInfo vector_variable_group_info(const Vec &vec);
VecParameterGroupInfo vector_parameter_group_info(const Vec &vec);

} // namespace ad::detail

#endif // MOO_AD_DETAIL_TRAVERSAL_H
