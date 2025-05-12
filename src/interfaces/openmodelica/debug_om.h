#ifndef OPT_OM_PRINTS
#define OPT_OM_PRINTS

#include "simulation_data.h"

void print_real_var_names(DATA* data);
void print_real_var_names_values(DATA* data);
void print_jacobian_sparsity(const JACOBIAN* jac, bool print_pattern, const char* name);

#endif // OPT_OM_PRINTS