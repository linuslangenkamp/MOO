#ifndef OPT_OM_PRINTS_H
#define OPT_OM_PRINTS_H

#include "simulation_data.h"
#include "simulation/solver/model_help.h"

#include <base/util.h>
#include <base/nlp_structs.h>


namespace OpenModelica {

void print_real_var_names(DATA* data);
void print_parameters(DATA* data);
void print_real_var_names_values(DATA* data);
void print_jacobian_sparsity(const JACOBIAN* jac, bool print_pattern, const char* name);
void print_bounds_fixed_vector(FixedVector<Bounds>& vec);

} // namespace OpenModelica

#endif // OPT_OM_PRINTS_H