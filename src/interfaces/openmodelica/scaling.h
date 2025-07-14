#ifndef OPT_OM_SCALING_GDOP_H
#define OPT_OM_SCALING_GDOP_H

#include <nlp/instances/gdop/gdop.h>

#include "info_gdop.h"

namespace OpenModelica {

NominalScaling create_gdop_nominal_scaling(GDOP::GDOP& gdop, InfoGDOP& info);

} // namespace OpenModelica

#endif // OPT_OM_SCALING_GDOP_H
