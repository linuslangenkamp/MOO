// SPDX-License-Identifier: LGPL-3.0-or-later
//
// This file is part of MOO - Modelica / Model Optimizer
// Copyright (C) 2026 University of Applied Sciences and Arts
// Bielefeld, Faculty of Engineering and Mathematics
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef MOO_NLP_C_STRUCTURES_H
#define MOO_NLP_C_STRUCTURES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <interfaces/gdop/structures.h>

typedef struct c_nlp_callbacks_t {
    void (*eval_all)(const f64* x, const f64* rp, f64* out, void* user_data);
    void (*jac_all)(const f64* x, const f64* rp, f64* out, void* user_data);
    void (*hes_all)(const f64* x, const f64* rp, const f64* lambda, f64 obj_factor, f64* out, void* user_data);
} c_nlp_callbacks_t;

typedef struct c_nlp_problem_t {
    const int x_size;
    const int rp_size;
    const int g_size;

    f64* rp;
    f64* x0;
    bounds_t* x_bounds;
    bounds_t* g_bounds;

    f64 obj_nominal;
    f64* x_nominal;
    f64* g_nominal;

    coo_t* obj_jac;
    coo_t* g_jac;
    coo_t* hes;

    c_nlp_callbacks_t* callbacks;
    solver_ctx_t* solver_ctx;
    void* user_data;
} c_nlp_problem_t;

#ifdef __cplusplus
}
#endif

#endif
