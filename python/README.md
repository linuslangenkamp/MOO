# MOO Python Frontend

The `moo` Python package is a modeling and code-generation frontend for MOO.
It lets you describe optimization problems in Python, generate C callbacks with
MOO's in-tree AD layer, compile the generated problem against `libmoo`, run the
solver, and read the emitted result CSV files.

## Supported Problems

Use the factory that matches the problem you want to solve:

```python
from moo import gdop_model, nlp_model
```

- `gdop_model(...)`: direct-collocation dynamic optimization problems.
- `nlp_model(...)`: standard finite-dimensional NLPs.

`Model` and `GDOPModel` remain importable for compatibility, but new code
should prefer `gdop_model(...)`.

## Build Requirements

The generated C files link against a built MOO library. The Python codegen path
also requires the optional AD Python bindings:

```bash
cmake -S . -B build -DMOO_WITH_PYTHON_BINDINGS=ON
cmake --build build --target moo mooad _ad
PYTHONPATH=python python3 examples/moo/hello.py
```

The source checkout helper in `python/moo/__init__.py` automatically adds
`build/python/moo` to the package path when it exists, so `PYTHONPATH=python`
is enough for local development.

The Python package requires Python 3.10 or newer. `matplotlib` is a required
Python dependency for plotting. `pandas` is optional and is only needed when
requesting DataFrames from result tables.

For a local editable install in a virtual environment:

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -e .
```

## Common Workflow

For the usual generate-compile-solve workflow, use `run(...)`:

```python
out = "build/moo/my_problem"
run = model.run(out)
```

The lower-level steps are still available for codegen and benchmarking examples:

```python
c_path, h_path = model.generate(out)
exe_path = model.compile(out)
run = model.optimize(out)
```

`run()` is a convenience alias for `optimize(..., generate=True)`.
`optimize()` expects generated files to already exist unless `generate=True`
is passed.

Solver output is streamed to the Python process by default. Use `capture=True`
when you want to inspect stdout and stderr from Python:

```python
run = model.run("build/moo/my_problem", capture=True)
print(run.process.stdout)
print(run.returncode)
```

## Solver Settings

Every model supports a shared solver settings API:

```python
model.solver(
    backend="Ipopt",
    tolerance=1e-10,
    iterations=5000,
    derivative_test=False,
    linear_solver="MUMPS",
    hessian="Exact",
    jacobian="Exact",
    gradient="Exact",
    uno_preset="ipopt",
    cpu_time=3600.0,
    warm_start=False,
    qp=False,
)
```

You can override the solver per run:

```python
run = model.run("build/moo/my_problem", solver="Uno")
```

Additional command-line options can be forwarded with `solver_args`.
Explicit options in `solver_args` override values from `model.solver(...)`:

```python
run = model.run(
    "build/moo/my_problem",
    solver="Ipopt",
    solver_args=["--IpoptDerivativeTest=true"],
)
```

Uno selection is wired through the shared solver interface. Solver-backend
availability still depends on how MOO was configured and built.

## Codegen Strategy

Generated derivative callbacks use direct sparse AD graph functions by default. This
means the generated C writes Jacobian and Hessian sparse buffers directly
instead of evaluating a full JVP/HVP once for every sparse entry. For very
large sparse blocks, `auto` can switch to colored compressed JVP/HVP evaluation
when the color count is much smaller than the requested sparse buffer.

```python
model.codegen("auto")     # choose direct or colored derivative callbacks
model.codegen("direct")   # force direct sparse derivative callbacks
model.codegen("colored")  # force colored compressed JVP/HVP callbacks
```

Explicit mapped NLP blocks remain C loops because they are model structure, not
a codegen mode. The derivative strategy only controls how each local graph
function fills Jacobian/Hessian buffers.
For tiny local blocks, `auto` usually chooses direct sparse graph functions; for larger
sparse local blocks, forced or automatic colored graph functions use compressed JVP/HVP
evaluation inside the loop.

Each `generate(...)` call writes `codegen_report.txt` next to the generated C
file. The report includes derivative nnz, coloring counts, selected strategy,
and generated C size.

## Shared Local Matrix And Vector Expressions

NLP and GDOP use the same local expression layer. Matrix and vector
expressions are therefore available in raw NLP constraints and GDOP
Lagrange/dynamics/path/Mayer/boundary graph functions:

```python
from moo import matrix, vec

M = matrix([[1.0, 2.0], [3.0, 4.0]])
y = M @ vec([x0, x1])
```

Use `sparse_matrix(...)` for sparse constant operators, and
`matrix(D).otimes_eye(nx) @ blockvec(...)` for block-state operators.

## Structured NLP Blocks

For large NLPs, write repeated structure explicitly so generated C keeps C
loops around repeated local graph functions:

```python
from moo import nlp_model

m = nlp_model("banded")
x = m.add_variables("x", 500, lb=-10.0, ub=10.0)

m.minimize_sum(range(500), lambda i: (x[i] - 1.0) ** 2)
m.add_constraints(
    range(499),
    lambda i: x[i] + x[i + 1],
    lb=-5.0,
    ub=5.0,
    name="band",
)
m.codegen("direct")
```

Variable vectors support normal Python views:

```python
x[5:]
x[::2]
x[5:20:3]
```

Inside mapped blocks, use `vec`, `blockvec`, `vector`, and `matrix` for local
vector algebra. `matrix(D).otimes_eye(nx)` is represented as a structured
Kronecker operator, so it does not materialize the zero-heavy expanded matrix
in Python:

```python
from moo import blockvec, matrix, vec

D_otimes_I = matrix(D).otimes_eye(nx)

def nodes(k):
    return blockvec([x_prev(k)] + [x_stage(k, j) for j in range(stages)])

def defect(k, j):
    X = nodes(k)
    dX = D_otimes_I @ X
    return dX[j + 1] - h * f(X[j + 1], u_stage(k, j))

m.add_constraints(
    range(intervals),
    lambda k: [value for j in range(stages) for value in defect(k, j)],
    eq=0.0,
    name="collocation",
)
```

This is intentionally generic. MOO does not require a collocation-specific
helper; Radau, finite-difference, shooting, and hand-written sparse NLPs can use
the same expression, matrix, block-vector, and mapped-constraint primitives.

## GDOP Example

```python
from moo import gdop_model

m = gdop_model("hello")
x = m.add_state("x", start=1.0, final=0.0)
u = m.add_control("u")

m.set_time_fixed(t0=0.0, tf=1.0)
m.mesh(intervals=25, nodes=3)
m.set_dynamics(x, x + u)
m.add_lagrange(u * u + x * x)
m.solver(tolerance=1e-12, derivative_test=True)

run = m.run("build/moo/hello")
print(run.result.states["x"])
print(run.result.controls["u"])
```

GDOP models support:

- states, controls, static parameters, and runtime parameters;
- state bounds, control bounds, parameter bounds, guesses, and nominals;
- fixed endpoint state values;
- fixed or free initial/final time;
- Lagrange and Mayer objective terms;
- dynamics, path constraints, and boundary constraints;
- sparse Jacobian and staged Hessian callback codegen.

Use explicit time setup:

```python
m.set_time_fixed(t0=0.0, tf=1.0)
```

or:

```python
m.set_time_free(
    t0_guess=0.0,
    tf_guess=1.0,
    t0_bounds=(0.0, 2.0),
    tf_bounds=(1.0, 5.0),
)
```

`m.time` is available in Lagrange/dynamics/path expressions. `m.t0` and
`m.tf` are available in Mayer and boundary expressions.

Mesh setup is:

```python
m.mesh(intervals=25, nodes=3)
m.mesh_refinement(l2bn_p1_it=0, l2bn_p2_it=0, l2bn_p2_lvl=0.0)
```

`set_time(...)` and `set_mesh(...)` remain as deprecated compatibility
wrappers.

## Standard NLP Example

```python
from moo import nlp_model

m = nlp_model("qp")
x = m.add_variable("x", lb=-10.0, ub=10.0, guess=0.0)
y = m.add_variable("y", lb=-10.0, ub=10.0, guess=0.0)

m.minimize((x - 1.0) ** 2 + (y - 2.0) ** 2)
m.add_constraint(x + y, eq=3.0, name="sum")

run = m.run("build/moo/nlp_qp")
print(run.result.variables)
print(run.result.constraints)
```

NLP models support variables, bounds, guesses, nominals, runtime parameters,
one objective, bounded constraints, exact sparse Jacobians, and exact staged
Hessian callbacks.

For large repeated NLPs, use iterable structured helpers so MOO generates
compact C loops instead of fully unrolled scalar code:

```python
from moo import matrix, nlp_model, vec, vector

m = nlp_model("banded")
x = m.add_variables("x", 500, lb=-10.0, ub=10.0, guess=0.1)

m.minimize_sum(range(500), lambda i: (x[i] - 1.0) ** 2, name="tracking")
m.minimize_sum(range(499), lambda i: 0.01 * x[i] * x[i + 1], name="coupling")
m.add_constraints(range(499), lambda i: x[i] + x[i + 1], lb=-5.0, ub=5.0, name="band")

m.codegen("auto").generate("build/moo/banded")
```

Alternatively an integer can be passed, e.g. `499` instead of `range(499)`.
You can also pass stepped ranges or explicit index lists:

```python
x_even = x[::2]
x_tail = x[5:]

x_even.fix(0, 1.0)
x_tail.set_nominal(slice(None), 10.0)

m.add_constraints(range(5, 100, 3), lambda i: x[i] - x[i - 1], eq=0.0)
m.add_constraints([0, 2, 7], lambda i: x_even[i] + x[i + 1], lb=-5.0, ub=5.0)
```

Mapped block bodies may return scalar expressions, lists/tuples, or expression
vectors. We also include tiny matrix/vector wrappers that allow for nice block-level operations.

```python
D = matrix([[-1.0, 1.0]])
w = vector([1.0])

def f(xb, ub):
    return vec([
        -(ub[0] + ub[0] * ub[0] / 2) * xb[0],
        ub[0] * xb[0],
    ])

m.add_constraints(
    range(n),
    lambda k: vec([
        (D @ vec([x1[k], x1[k + 1]]))[0],
        (D @ vec([x2[k], x2[k + 1]]))[0],
    ]) - h * f(vec([x1[k], x2[k]]), vec([u[k]])),
    eq=0.0,
    name="collocation",
)

m.minimize_sum(range(n), lambda k: h * (w @ vec([u[k] * u[k]])))
```

Radau IIA, i.e. rescaled flipped Legendre-Gauss-Radau, constants from `data/fLGR/` are available directly from the package for
stages 1 through 25:

```python
from moo import matrix, radauIIA, vector

r = radauIIA(3)
D = matrix(r.D)     # differentiation matrix, including the previous value column
B = vector(r.B)     # quadrature weights
C = r.C             # collocation nodes
C0 = r.C0           # nodes including 0
W = r.W             # barycentric weights
W0 = r.W0           # barycentric weights including 0
```

For vector-valued states, `D.otimes_eye(nx)` creates a structured
`D otimes I_nx` operator. It works directly on `blockvec(...)`, so the Python
model does not materialize the expanded zero-heavy Kronecker matrix:

```python
D_otimes_I = matrix(r.D).otimes_eye(2)
X = blockvec([x_previous, x_stage_0, x_stage_1, x_stage_2])
dX = D_otimes_I @ X
defect_j = dX[j + 1] - h * f(X[j + 1], u_j)
```

Mapped objective and constraint blocks emit one local AD graph function and C loops
for value, Jacobian, and Hessian accumulation. `codegen_report.txt` lists the
mapped block count, repetitions, global derivative nnz, and selected `loop`
strategy.

Multi-output mapped constraints are useful when one repeated block has a large
local sparse Jacobian:

```python
m.add_constraints_map(
    blocks,
    lambda i: [(x[i + j] + x[i + j + 1]) ** 2 for j in range(block_size - 1)],
    lb=-1.0,
    ub=25.0,
    name="band_block",
)
m.codegen("colored")
```

## Direct Collocation

As can be seen from the previous section, direct collocation can also be implemented with the native NLP.
Although quite comfortable, modeling problems and direct collocation with this raw NLP interface is significantly harder.
Examples of direct collocation with the native NLP can be found in `examples/moo/nlp_collocation_batch_reactor.py` and
`examples/moo/nlp_collocation_batch_reactor_piecewise_ctrls.py`. This interface gives full control over all details of the problem, e.g. if controls are
piecewise constant or if the problem is a multi-phase optimization problem, but also requires the entire NLP and not only the
model equations to be generated to C code. This is done very efficiently with the described AD implementation, but also takes additional time
that is saved using the GDOP layer. Furthermore, the GDOP layer implements several pluggable strategies, e.g. mesh refinement, which is not implemented
with the raw NLP layer. So the decision between NLP and GDOP really comes down to your use case.

## Results

`optimize()` and `run()` return an `OptimizationRun`:

```python
run = model.run("build/moo/hello")
print(run.returncode)
print(run.process.args)
```

It contains:

- `run.process`: the `subprocess.CompletedProcess`;
- `run.results`: all CSV files parsed into a `ResultSet`;
- `run.result`: a typed result view for the problem type;
- `run.read_results()`: reload CSV files from the run directory;
- `run.table(name)`: return a stdlib-backed result table;
- `run.dataframe(name)`: return a pandas DataFrame if pandas is installed.

Read CSV output directly with:

```python
from moo import read_results

results = read_results("build")
primals = results.table("primals_optimal_solution")
print(primals.column("time"))
```

Typed result views expose common data by model variable names:

```python
gdop = run.result
gdop.time
gdop.states["x"]
gdop.controls["u"]
gdop.costates

init = run.result
init.objective
init.variables
init.parameters
init.deltas
init.f
init.g

nlp = run.result
nlp.objective
nlp.variables
nlp.constraints
```

## Plotting

Result objects include a plotting accessor:

```python
fig = run.result.plot.all()
fig = run.result.plot.states()
fig = run.result.plot.controls()
fig = run.result.plot.costates()
fig = run.result.plot.variables()
fig = run.result.plot.constraints()
```

GDOP results use time-series plots. NLP results use bar plots for
flat variable and constraint data. Plot methods return matplotlib figures.
In a script, returning a figure does not open a window by itself. Use
`show=True` for an interactive matplotlib window:

```python
run.result.plot.all(show=True)
```

or save the figure to a file:

```python
run.result.plot.all(save="build/moo/hello/solution.png")
```

Both can be combined:

```python
run.result.plot.all(save="build/moo/hello/solution.png", show=True)
```

## AD And Code Generation

Python model codegen uses MOO's native AD bindings directly.

The generated C callbacks still compile into a standalone executable that
links against `libmoo` and enters the matching MOO C interface:

- `src/interfaces/c` for GDOP;
- `src/interfaces/nlp` for standard NLP.

The AD layer emits:

- value functions;
- direct sparse Jacobian/Hessian functions by coefficient extraction;
- JVP and staged HVP functions for colored compressed filling;
- structural sparsity used to allocate C callback buffers.

For direct Hessian callbacks, generated C evaluates only the requested sparse
Hessian coefficients. For colored callbacks, generated C prepares the staged
HVP cache once at the fixed primal/lambda point and reuses it for all
directions.

## Low-Level AD Bindings

The `moo.ad` module exposes graph construction, AD transforms, VM evaluation,
staged VM evaluation, sparsity queries, and C emission:

```python
from moo import ad

g = ad.GraphFunctionBuilder()
x = g.inputs("x", 3)
f = g.function(x, g.vector([2.0 * x[0] - x[1] + 3.0, x[2] * x[2] + 1.0]))

grad = f.reverse_diff("lambda", "x")
hvp = grad.forward_diff("x", "v")

print(f.jacobian_sparsity("x"))
print(hvp.hessian_sparsity("v"))
print(f.evaluate(inputs={"x": [1.0, 2.0, 3.0]}))
print(grad.evaluate(inputs={"x": [1.0, 2.0, 3.0]}, params={"lambda": [1.0, 1.0]}))

prepared = hvp.prepare(
    inputs={"x": [1.0, 2.0, 3.0]},
    params={"lambda": [1.0, 1.0]},
)
print(prepared.apply([0.0, 0.0, 1.0]))
print(f.to_c("demo_value"))
```

See `src/ad/README.md` for the C++ and Python AD API.
See `reference/gdop.md`, `reference/init.md`, and `reference/nlp.md` for the
lower-level C++ and C interface references.

## Examples

- `examples/moo/hello.py`: small GDOP with derivative checking.
- `examples/moo/representative.py`: GDOP using states, controls, a parameter,
  runtime parameters, Lagrange/Mayer objectives, path constraints, and boundary
  constraints.
- `examples/moo/free_time.py`: GDOP with a bounded free final time.
- `examples/moo/free_time_boundary.py`: focused free-time Mayer and boundary
  example.
- `examples/moo/linear_gdop.py`: linear GDOP regression; useful for checking
  zero Hessian sparsity.
- `examples/moo/nlp_qp.py`: minimal standard convex QP.
- `examples/moo/NLP/spectral_nlp_free_time.py`: spectral hypersensitive NLP
  with the final time as an NLP variable constrained by `tf >= 0.5`.
- `examples/moo/nlp_sparse_benchmark.py`: sparse NLP code-size comparison for
  auto, direct sparse, and colored compressed derivative codegen.
- `examples/moo/ad_bindings.py`: low-level AD binding demo.

## Current Scope

This frontend is usable for small and medium generated problems and for
developing MOO's native AD/codegen path. Sparse Jacobian and Hessian callbacks
are exact. The default path emits direct sparse derivative graph functions;
colored compressed callbacks are available through `model.codegen(...)`.

External data callbacks, richer solver-specific option objects, and more
advanced mesh-refinement controls are not yet exposed as high-level Python
APIs.
