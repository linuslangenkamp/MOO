# Standard NLP

MOO's standard NLP interface solves finite-dimensional nonlinear programs of
the form:

```math
\begin{aligned}
\min_x\quad & f(x) \\
\text{s.t.}\quad & g_{\mathrm{lb}} \le g(x) \le g_{\mathrm{ub}}, \\
& x_{\mathrm{lb}} \le x \le x_{\mathrm{ub}}.
\end{aligned}
```

Use this problem type when the optimization variables are already finite
dimensional and no dynamic transcription or initialization-specific parameter
correction is needed.

## Python Usage

The Python frontend exposes standard NLPs through `nlp_model(...)`:

```python
from moo import nlp_model

m = nlp_model("qp")
x = m.add_variable("x", lb=-10.0, ub=10.0, guess=0.0, nominal=1.0)
y = m.add_variable("y", lb=-10.0, ub=10.0, guess=0.0, nominal=1.0)

m.minimize((x - 1.0) ** 2 + (y - 2.0) ** 2)
m.add_constraint(x + y, eq=3.0, name="sum")
m.solver(tolerance=1e-12, derivative_test=True)

run = m.run("build/moo/nlp_qp")
print(run.result.objective)
print(run.result.variables)
print(run.result.constraints)
```

Supported Python modeling features:

- variables with bounds, guesses, and nominals;
- runtime parameters;
- one scalar objective;
- bounded constraints or equality constraints with `eq=...`;
- objective nominal and constraint nominals;
- exact AD-generated objective gradient, constraint Jacobian, and Lagrangian
  Hessian callbacks;
- typed result access and plotting.

## C Interface

Generated Python problems target the C interface in `src/interfaces/nlp`.
Custom C frontends can use the same interface.

The central structure is `c_nlp_problem_t`:

```c
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
```

The callbacks are:

```c
void eval_all(const f64* x, const f64* rp, f64* out, void* user_data);
void jac_all(const f64* x, const f64* rp, f64* out, void* user_data);
void hes_all(const f64* x,
             const f64* rp,
             const f64* lambda,
             f64 obj_factor,
             f64* out,
             void* user_data);
```

`eval_all` writes `[objective, g...]`. `jac_all` writes all declared objective
gradient and constraint Jacobian values into one flat buffer. `hes_all` writes
the lower-triangular Hessian of the Lagrangian:

```math
\sigma \nabla^2 f(x) + \lambda^\mathsf{T}\nabla^2 g(x)
```

The `coo_t` structures define row indices, column indices, and buffer indices.
This lets a generated or custom frontend choose its own derivative buffer order
while MOO assembles the sparse NLP layout expected by the solver.

## Running

The C entry point is:

```c
int main_nlp(int argc, char** argv, c_nlp_problem_t* problem);
```

Generated Python code creates a small executable whose `main(...)` calls
`main_nlp(...)`. Solver options are passed as command-line arguments, for
example:

```bash
./build/moo/nlp_qp/nlp_qp_codegen --NLPSolver=Ipopt --Tolerance=1e-12
```

After optimization, the interface writes `nlp_optimal_solution.csv` in the run
directory. The Python frontend parses this into `NLPResult`.

## Results

`NLPResult` exposes:

- `objective`;
- `variables`, keyed by Python variable name;
- `constraints`, keyed by Python constraint name;
- `lambdas`;
- `bound_duals`;
- raw CSV access through `table(...)` and `dataframe(...)`.

Plot variables or constraints with:

```python
run.result.plot.variables()
run.result.plot.constraints()
run.result.plot.all()
```
