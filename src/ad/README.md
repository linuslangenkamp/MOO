# MOO AD

MOO AD is the in-tree automatic differentiation and C-code-emission layer used
by the Python modeling frontend. It is implemented by the C++17 header
`ad.h` and optional pybind11 bindings exposed as `moo.ad`.

The library builds symbolic expression graphs for vector-valued functions,
transforms those graphs with forward and reverse AD, evaluates them with a VM,
queries structural sparsity, and emits C kernels that can be linked into MOO
problem interfaces.

## What It Provides

- Expression graph construction for vector-valued functions.
- Named input groups and named parameter groups.
- Graph optimization: constant folding, algebraic simplification, common
  subexpression elimination, and dead-code elimination.
- Symbolic forward mode with `forward_diff`.
- Symbolic reverse mode with `reverse_diff`.
- Higher-order derivatives by composing derivative graphs.
- Forward-over-reverse Hessian-vector products.
- VM evaluation.
- Staged VM evaluation for repeated HVP directions.
- One-shot C code emission with `to_c`.
- Staged C code emission with separate `_prepare` and `_apply` functions.
- Direct sparse derivative C emission for selected Jacobian/Hessian entries.
- Greedy column coloring for compressed derivative evaluation.
- Structural Jacobian and Hessian sparsity from optimized derivative graphs.
- Python bindings for graph construction, AD transforms, evaluation, sparsity,
  and C emission.

## C++ Quick Start

```cpp
#include "ad.h"

using namespace ad;

Graph g;
auto x = g.inputs("x", 2);

std::vector<Expr> y = {
    sin(x[0]),
    x[0] * x[0] + exp(x[1]),
};

GraphFunction F = function_from(std::move(g), x, y);
```

A `GraphFunction` is a first-class value. You can evaluate it, differentiate
it, emit it as C, or use it for sparsity analysis.

## Differentiation

```cpp
// Directional derivative dF/dx * v.
GraphFunction JVP = forward_diff(F, "x", "v");

// Gradient of lambda^T F with respect to x.
GraphFunction Grad = reverse_diff(F, "lambda", "x");

// Hessian-vector product of lambda^T F.
GraphFunction HVP = forward_diff(Grad, "x", "v");
```

The group names matter:

- `"x"` is the differentiated input group;
- `"v"` is the JVP direction parameter group;
- `"lambda"` is the reverse-mode seed parameter group.

## C++ VM Evaluation

```cpp
VM vm(HVP);

double xv[2] = {1.0, 2.0};
double lambda[2] = {3.0, 4.0};
double v[2] = {0.1, -0.2};
double out[2];

EvalEnv env;
env.input("x", xv)
   .param("lambda", lambda)
   .param("v", v);

vm.evaluate(env, out);
```

`EvalEnv` binds named input and parameter groups to raw arrays. The VM checks
that every required group is present during evaluation.

## C++ Staged VM Evaluation

Staged evaluation is useful for Hessian-vector products. If the primal point
and multipliers are fixed while the direction `v` changes, prepare once and
apply many times:

```cpp
StagedVM staged(HVP, "v");

EvalEnv prep_env;
prep_env.input("x", xv)
        .param("lambda", lambda);

auto prepared = staged.prepare(prep_env);

for (auto& direction : directions) {
    prepared.apply(direction.data(), out);
}
```

`prepare(...)` computes all nodes independent of the chosen direction group.
`apply(...)` evaluates only the direction-dependent nodes.

## C Code Emission

One-shot emission:

```cpp
std::cout << to_c(F, "value");
```

emits:

```c
#include <math.h>

void value(const double* x, double* out);
```

For derivative graphs with parameters, the emitted function signature includes
those groups in graph order:

```c
void hvp_eval(const double* x,
              const double* lambda,
              const double* v,
              double* out);
```

Staged emission:

```cpp
std::cout << to_staged_c(HVP, "hvp", "v");
```

emits a cache type and two functions:

```c
typedef struct { ... } hvp_cache_t;

void hvp_prepare(const double* x,
                 const double* lambda,
                 hvp_cache_t* cache);

void hvp_apply(const hvp_cache_t* cache,
               const double* v,
               double* out);
```

MOO uses this staged form for colored compressed Hessian callbacks and for the
debug basis fallback. Direct sparse Hessian callbacks use coefficient
extraction instead and emit the requested sparse values directly.

## Structural Sparsity

Sparsity is computed on optimized graphs. Structural zeros that simplify to
constant zero are removed before sparsity is queried.

```cpp
auto J = jacobian_sparsity(F, "x");
auto H_lower = hessian_sparsity(HVP, "v");
auto H_full = hessian_sparsity_full(HVP, "v");
```

You can also call the lower-level reachability query:

```cpp
auto pattern = structural_sparsity(JVP, "v");
```

For derivative graphs, probe the linear seed or direction variable:

| Graph | Probe group | Gives |
| --- | --- | --- |
| `F` | input group `"x"` | Jacobian sparsity of `F` with respect to `x` |
| `JVP = forward_diff(F, "x", "v")` | parameter group `"v"` | Jacobian sparsity |
| `Grad = reverse_diff(F, "lambda", "x")` | parameter group `"lambda"` | transposed Jacobian sparsity |
| `HVP = forward_diff(Grad, "x", "v")` | parameter group `"v"` | Hessian sparsity |

`hessian_sparsity(...)` returns lower-triangular `(row, col)` pairs with
`row >= col`. `hessian_sparsity_full(...)` returns the symmetric full pattern.

The `SparsityPattern` API includes:

```cpp
pat.nnz()
pat.to_pairs()
pat.contract_outputs()
SparsityPattern::lower_triangular(pairs)
SparsityPattern::symmetrize(pairs)
```

## Direct Sparse Derivative Emission

For sparse solver callbacks, emit derivative buffers directly:

```cpp
auto J = jacobian_sparsity(F, "x").to_pairs();
auto Grad = reverse_diff(F, "lambda", "x");
auto HVP = forward_diff(Grad, "x", "v");
auto H = hessian_sparsity(HVP, "v");

std::cout << to_sparse_jacobian_c(F, "x", J, "jac_values");
std::cout << to_sparse_hessian_c(HVP, "v", H, "hes_values");
```

These functions extract coefficients from graphs that are linear in a
direction/seed parameter. The generated C outputs exactly the requested sparse
entries in the requested COO order. This is the preferred path for MOO's
generated NLP, Init, and GDOP callbacks.

## Coloring

For large sparse derivative buffers, MOO can also evaluate compressed
directions. The AD layer exposes greedy column coloring over a sparsity pattern:

```cpp
auto J = jacobian_sparsity(F, "x").to_pairs();
auto colors = greedy_column_coloring(J, F.input_size("x"));
```

The Python binding exposes the same operation through `GraphFunction.coloring`.
MOO's high-level code generator uses this metadata for
`model.codegen("colored")` and may choose it automatically for large sparse
blocks. Colored callbacks seed all columns of one color, run one JVP or staged
HVP apply, and scatter the requested sparse entries through static metadata
arrays. Direct sparse coefficient kernels remain the default for small and
medium blocks because they produce smaller, faster C in those cases.

## Python Bindings

Build the optional bindings with:

```bash
cmake -S . -B build -DMOO_WITH_PYTHON_BINDINGS=ON
cmake --build build --target _ad
PYTHONPATH=python python3 examples/moo/ad_bindings.py
```

The Python API mirrors the C++ graph API:

```python
from moo import ad

g = ad.GraphBuilder()
x = g.inputs("x", 3)
f = g.function(x, [2.0 * x[0] - x[1] + 3.0, x[2] * x[2] + 1.0])

grad = f.reverse_diff("lambda", "x")
hvp = grad.forward_diff("x", "v")
```

Evaluate with the VM:

```python
print(f.evaluate(inputs={"x": [1.0, 2.0, 3.0]}))
print(
    grad.evaluate(
        inputs={"x": [1.0, 2.0, 3.0]},
        params={"lambda": [1.0, 1.0]},
    )
)
print(
    hvp.evaluate(
        inputs={"x": [1.0, 2.0, 3.0]},
        params={"lambda": [1.0, 1.0], "v": [0.0, 0.0, 1.0]},
    )
)
```

Use staged VM evaluation when only the direction changes:

```python
prepared = hvp.prepare(
    inputs={"x": [1.0, 2.0, 3.0]},
    params={"lambda": [1.0, 1.0]},
)

print(prepared.apply([0.0, 0.0, 1.0]))
```

A one-shot staged helper is also available:

```python
print(
    hvp.staged_apply(
        [0.0, 0.0, 1.0],
        inputs={"x": [1.0, 2.0, 3.0]},
        params={"lambda": [1.0, 1.0]},
    )
)
```

Query sparsity and emit C:

```python
print(f.jacobian_sparsity("x"))
print(hvp.hessian_sparsity("v"))
print(f.to_c("demo_value"))
print(hvp.to_staged_c("demo_hvp", "v"))
print(f.to_sparse_jacobian_c("x", f.jacobian_sparsity("x"), "demo_jac"))
print(hvp.to_sparse_hessian_c("v", hvp.hessian_sparsity("v"), "demo_hes"))
print(f.coloring(f.jacobian_sparsity("x"), 3))
```

Python dictionaries are keyed by group name. Values must have the exact size
of the corresponding input or parameter group.

## Supported Operations

The expression layer currently supports:

```text
+  -  *  /  unary -  sin  cos  tan  exp  log  pow_const
```

In Python, exponentiation by a numeric constant uses `pow_const` internally:

```python
y = x[0] ** 2
```

## Building C++ Examples

`ad.h` is header-only:

```bash
g++ -std=c++17 -O2 -Isrc/ad src/ad/example_hvp.cpp -o example_hvp
./example_hvp
```

No third-party library is needed for the C++ header itself. pybind11 is only
needed for the optional Python module.
