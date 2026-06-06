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

#ifndef MOO_INIT_C_STRUCTURES_H
#define MOO_INIT_C_STRUCTURES_H

#ifdef __cplusplus
extern "C" {
#endif

#include <interfaces/gdop/structures.h>

typedef struct c_init_callbacks_t {
    void (*eval_all)(const f64* z, const f64* rp, f64* out, void* user_data);
    void (*jac_all)(const f64* z, const f64* rp, f64* out, void* user_data);
    void (*hes_all)(const f64* z, const f64* rp, const f64* lambda, f64 obj_factor, f64* out, void* user_data);
} c_init_callbacks_t;

typedef struct c_init_problem_t {
    const int y_size;
    const int p_size;
    const int rp_size;
    const int f_size;
    const int g_size;

    f64* rp;
    f64* y0;
    f64* p0;
    f64* dp0;

    bounds_t* y_bounds;
    bounds_t* p_bounds;
    bounds_t* dp_bounds;
    bounds_t* g_bounds;

    f64* y_nominal;
    f64* p_nominal;
    f64* dp_nominal;
    f64* f_nominal;
    f64* g_nominal;
    f64 obj_nominal;

    coo_t* obj_jac;
    coo_t* f_jac;
    coo_t* g_jac;
    coo_t* hes;

    c_init_callbacks_t* callbacks;
    solver_ctx_t* solver_ctx;
    void* user_data;
} c_init_problem_t;

#ifdef __cplusplus
}
#endif

#endif
