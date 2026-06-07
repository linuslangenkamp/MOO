# MOO Python Frontend

The `moo` Python package is a modeling and code-generation frontend for MOO.
It lets you describe optimization problems in Python, generate C callbacks with
MOO's in-tree AD layer, compile the generated problem against `libmoo`, run the
solver, and read the emitted result CSV files.

## Supported Problems

Use the factory that matches the problem you want to solve:

```python
from moo import gdop_model, init_model, nlp_model
```

- `gdop_model(...)`: direct-collocation dynamic optimization problems.
- `init_model(...)`: flat consistent-initialization NLPs over variables and
  corrected parameters.
- `nlp_model(...)`: standard finite-dimensional NLPs.

`Model` and `GDOPModel` remain importable for compatibility, but new code
should prefer `gdop_model(...)`.

## Build Requirements

The generated C files link against a built MOO library. The Python codegen path
also requires the optional AD Python bindings:

```bash
cmake -S . -B build -DMOO_WITH_PYTHON_BINDINGS=ON
cmake --build build --target moo _ad
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

All model types share the same lifecycle:

```python
out = "build/moo/my_problem"
c_path, h_path = model.generate(out)
exe_path = model.compile(out)
run = model.optimize(out)
```

For the usual generate-compile-solve workflow, use `run(...)`:

```python
run = model.run("build/moo/my_problem")
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

Generated derivative callbacks use direct sparse AD kernels by default. This
means the generated C writes Jacobian and Hessian sparse buffers directly
instead of evaluating a full JVP/HVP once for every sparse entry. For very
large sparse blocks, `auto` can switch to colored compressed JVP/HVP evaluation
when the color count is much smaller than the requested sparse buffer.

```python
model.codegen("auto")          # structure=auto, local=auto
model.codegen("loop")          # structure=loop, local=auto
model.codegen("direct")        # structure=auto, local=direct
model.codegen("colored")       # structure=auto, local=colored
model.codegen("basis")         # structure=auto, local=basis
model.codegen("loop-direct")   # structure=loop, local=direct
model.codegen("loop-colored")  # structure=loop, local=colored
model.codegen(structure="loop", local="colored")
```

The structure axis controls whether repeated blocks remain C loops. The local
axis controls how each repeated block's local derivative kernel is generated.
For tiny local blocks, `auto` usually chooses direct sparse kernels; for larger
sparse local blocks, forced or automatic colored kernels use compressed JVP/HVP
evaluation inside the loop.

Each `generate(...)` call writes `codegen_report.txt` next to the generated C
file. The report includes derivative nnz, coloring counts, selected strategy,
and generated C size.

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

## Init Example

```python
from moo import init_model

m = init_model("simple_init")
y = m.add_variable("y", lb=-10.0, ub=10.0, guess=0.0)
p = m.add_parameter("p", lb=0.0, ub=5.0, base=1.0)
dp = m.delta(p)

m.add_f(y + p - 3.0)
m.add_g(y, lb=0.0)
m.set_objective((y - 1.0) ** 2 + (dp - 1.0) ** 2)

run = m.run("build/moo/init_simple")
print(run.result.variables)
print(run.result.parameters)
print(run.result.deltas)
```

Init models support variables `y`, effective parameters `p`, parameter
corrections `dp = p - base`, equality residuals `f`, bounded constraints `g`,
runtime parameters, exact sparse Jacobians, and exact staged Hessian
callbacks.

For compatibility, a parameter can still be destructured:

```python
p, dp = m.add_parameter("p", base=1.0)
```

New code should prefer:

```python
p = m.add_parameter("p", base=1.0)
dp = m.delta(p)
```

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

For large repeated NLPs, use the structured helpers so MOO generates compact C
loops instead of fully unrolled scalar code:

```python
from moo import nlp_model

m = nlp_model("banded")
x = m.add_variables("x", 500, lb=-10.0, ub=10.0, guess=0.1)

m.minimize_sum(500, lambda i: (x[i] - 1.0) ** 2, name="tracking")
m.minimize_sum(499, lambda i: 0.01 * x[i] * x[i + 1], name="coupling")
m.add_constraints(499, lambda i: x[i] + x[i + 1], lb=-5.0, ub=5.0, name="band")

m.codegen("auto").generate("build/moo/banded")
```

Mapped objective and constraint blocks emit one local AD kernel and C loops
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
m.codegen("loop-colored")
```

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

GDOP results use time-series plots. Init and NLP results use bar plots for
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
- `src/interfaces/init` for Init;
- `src/interfaces/nlp` for standard NLP.

The AD layer emits:

- value functions;
- direct sparse Jacobian/Hessian functions by coefficient extraction;
- JVP and staged HVP functions for colored compressed or debug fallback filling;
- structural sparsity used to allocate C callback buffers.

For direct Hessian callbacks, generated C evaluates only the requested sparse
Hessian coefficients. For colored or basis callbacks, generated C prepares the
staged HVP cache once at the fixed primal/lambda point and reuses it for all
directions.

## Low-Level AD Bindings

The `moo.ad` module exposes graph construction, AD transforms, VM evaluation,
staged VM evaluation, sparsity queries, and C emission:

```python
from moo import ad

g = ad.GraphBuilder()
x = g.inputs("x", 3)
f = g.function(x, [2.0 * x[0] - x[1] + 3.0, x[2] * x[2] + 1.0])

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
- `examples/moo/init_simple.py`: minimal Init problem with solution `y = 1`,
  `p = 2`.
- `examples/moo/nlp_qp.py`: minimal standard convex QP.
- `examples/moo/nlp_sparse_benchmark.py`: sparse NLP code-size comparison for
  auto, direct sparse, colored compressed, and basis-loop derivative codegen.
- `examples/moo/ad_bindings.py`: low-level AD binding demo.

## Current Scope

This frontend is usable for small and medium generated problems and for
developing MOO's native AD/codegen path. Sparse Jacobian and Hessian callbacks
are exact. The default path emits direct sparse derivative kernels; colored
compressed callbacks and the old basis-direction fallback are available through
`model.codegen(...)`.

External data callbacks, richer solver-specific option objects, and more
advanced mesh-refinement controls are not yet exposed as high-level Python
APIs.
