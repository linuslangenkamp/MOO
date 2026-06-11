// SPDX-License-Identifier: LGPL-3.0-or-later
#ifndef MOO_AD_OPTIMIZE_H
#define MOO_AD_OPTIMIZE_H

#include "expr.h"
#include "function.h"
#include "vec.h"

namespace ad {

Expr optimize(const Expr &expr);
Vec optimize(const Vec &vec);
Function optimize(const Function &function);

} // namespace ad

#endif // MOO_AD_OPTIMIZE_H
