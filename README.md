# **MOO: Modelica / Model Optimizer - A Generic Framework for Dynamic Optimization**

TODO: refactor this and remove major errors
TODO: add installation?!

**MOO: Modelica / Model Optimizer** is a flexible and extensible framework for solving complex optimization problems. It is distributed under the **OSMC Public License (OSMC-PL) Version 1.2** and the **GNU Lesser General Public License (LGPL) Version 3**.

**MOO** provides a generic Nonlinear Programming (NLP) layer with built-in scaling support and a solver-agnostic interface. While primarily designed for dynamic optimization in the Modelica ecosystem (General Dynamic Optimization Problems (GDOPs), training of Physics-enhanced Neural ODEs (PeN-ODEs)), it is equally applicable to other domains, e.g. model predictive control (MPC).

## 0. Information

### 0.1 License

The **Modelica/Model Optimizer (MOO)** is distributed under the **OSMC Public License (OSMC-PL) Version 1.2** and the **GNU Lesser General Public License (LGPL) Version 3**.

Under OSMC-PL, the following usage modes are supported:  
- **GPL** – Any party may use and redistribute MOO under GPL Version 3.  
- **OSMC-Internal-EPL** – Level 1 members of OSMC may also use and redistribute under these conditions.  
- **OSMC-External-EPL** – Level 2 members of OSMC may also use and redistribute under these conditions.  

This program is distributed **without any warranty**, not even the implied warranty of merchantability or fitness for a particular purpose. See the full license text for detailed terms and conditions: [OSMC Public License 1.2](https://openmodelica.org/osmc-license)

Additionally, the author also redistributes **MOO** under the terms of the **GNU Lesser General Public License (LGPL) Version 3**, providing further flexibility for use, modification, and redistribution. See the full license text for detailed terms and conditions: [LGPL-3.0](https://www.gnu.org/licenses/lgpl-3.0.html)

### 0.2 Developers

The **Modelica/Model Optimizer (MOO)** is an open-source project. It is developed and maintained by the following contributors:

- Linus Langenkamp

Contributions from the community are welcome, and additional authors or maintainers may be listed as the project evolves. The development focuses on delivering a modular optimization runtime for OpenModelica but is designed for general applicability beyond this ecosystem.

## 1. Generic Nonlinear Programming (NLP)

The **MOO** framework represents optimization problems in a generic NLP form, enabling integration with a wide range of solvers. Its scaling mechanism improves numerical conditioning and solver performance. By decoupling problem definition from solver backends, **MOO** offers a clean API for dynamic optimization workflows, including mesh refinement and Physics-enhanced Neural ODE training.

### 1.1. NLP Problem Formulation

A general NLP problem is typically formulated as:

$$
\begin{aligned}
\min_{x} \,\, & f(x) \\
\text{s.t.} \,\, & g^{LB} \le g(x) \le g^{UB} \\
& x^{LB} \le x \le x^{UB}
\end{aligned}
$$

Where:

* $x$ is the vector of primal variables.

* $f(x)$ is the scalar objective function.

* $g(x)$ is the vector of constraint functions.

* $x^{LB}, x^{UB}$ are the lower and upper bounds on the primal variables.

* $g^{LB}, g^{UB}$ are the lower and upper bounds on the constraint functions.

The `NLP::NLP` class serves as the abstract base class for defining such a problem. Any specific NLP must derive from this class and implement its pure virtual methods.

### 1.2. NLP Solver Interface (`NLP::NLPSolver`)

The `NLP::NLPSolver` class defines a generic interface for any NLP solver to adhere to. This allows different numerical solvers to be easily integrated into the **MOO** framework.

* **`NLPSolver(NLP& nlp, NLPSolverSettings& solver_settings)`**: The constructor takes a reference to an `NLP` problem instance and solver-specific settings.

* **`virtual void optimize() = 0;`**: This pure virtual method is the core of the solver interface. Any concrete solver implementation must provide its own logic for initiating and running the optimization process when this method is called.

The `NLP::NLP` class itself acts as a contract between the problem definition and any NLP solver. It provides a standardized interface for solvers to query problem dimensions, bounds, sparsity patterns, and to request function and derivative evaluations.

The `NLP` class distinguishes between:

* **User Callbacks (Pure Virtual Methods)**: These are methods that a user *must* implement in their derived NLP class (e.g., `GDOP::GDOP`) to define the specific problem. They typically compute the unscaled values of the objective, constraints, and their derivatives.

* **Solver API (Public Methods prefixed with `solver_`)**: These are methods called by the NLP solver. They handle the necessary data transformations (like scaling/unscaling) before calling the user's callbacks and then transform the results back for the solver.

Key methods in the `NLP` interface include:

* **`get_sizes(int& number_vars, int& number_constraints)`**: User callback to specify the total number of primal variables and constraints.

* **`get_nnz(int& nnz_jac, int& nnz_hes)`**: User callback to provide the number of non-zero elements in the Jacobian and Hessian.

* **`get_bounds(FixedVector<f64>& x_lb, FixedVector<f64>& x_ub, FixedVector<f64>& g_lb, FixedVector<f64>& g_ub)`**: User callback to define the lower and upper bounds for variables and constraints.

* **`get_initial_guess(...)`**: User callback to provide an initial guess for primal variables, dual variables (Lagrange multipliers), and dual multipliers for variable bounds.

* **`get_jac_sparsity(FixedVector<int>& i_row_jac, FixedVector<int>& j_col_jac)`**: User callback to define the sparsity pattern (row and column indices in COO format) of the Jacobian of the constraints.

* **`get_hes_sparsity(FixedVector<int>& i_row_hes, FixedVector<int>& j_col_hes)`**: User callback to define the sparsity pattern (row and column indices in COO format) of the Hessian of the Lagrangian.

* **`eval_f(...)`, `eval_g(...)`, `eval_grad_f(...)`, `eval_jac_g(...)`, `eval_hes(...)`**: User callbacks for evaluating the objective function, constraint functions, objective gradient, constraint Jacobian, and Lagrangian Hessian, respectively.

* **`finalize_solution(...)`**: User callback to process the optimal solution returned by the solver.

The `solver_` prefixed methods (e.g., `solver_get_info`, `solver_eval_f`) internally manage the data flow and scaling, ensuring that the solver receives and provides data in its expected format, while the user only needs to implement the core mathematical evaluations.

### 1.3. IPOPT Solver Implementation

The `IpoptSolver::IpoptSolver` class is a concrete implementation of the `NLP::NLPSolver` interface, providing the bridge to the IPOPT (Interior Point OPTimizer) library.

* **`IpoptSolver(NLP::NLP& nlp, NLP::NLPSolverSettings& solver_settings)`**: The constructor initializes the solver with the specific NLP problem and settings.

* **`Ipopt::SmartPtr<IpoptAdapter> adapter;`**: This is a smart pointer to an `IpoptAdapter` instance. The `IpoptAdapter` is the crucial component that translates calls from the IPOPT library (which expects a `TNLP` interface) into calls to our generic `NLP::NLP` interface.

* **`Ipopt::SmartPtr<Ipopt::IpoptApplication> app;`**: This manages the IPOPT application instance, responsible for setting up the solver, configuring options, and initiating the optimization.

* **`void optimize() override;`**: This method, implemented by `IpoptSolver`, sets up the IPOPT application, configures its settings (via `set_settings`), and then calls `app->OptimizeTillConvergence()` using the `adapter` to solve the NLP.

#### 1.3.1. `IpoptAdapter` (Bridge to IPOPT's TNLP)

The `IpoptSolver::IpoptAdapter` class is the key component that enables `NLP::NLP` problems to be solved by IPOPT. It inherits from IPOPT's `Ipopt::TNLP` (Truncated Nonlinear Program) abstract base class and implements all its required virtual methods. Each of these methods acts as a wrapper, translating IPOPT's requests into calls to the corresponding `NLP::NLP`'s `solver_` methods.

Key methods implemented by `IpoptAdapter`:

* **`get_nlp_info(...)`**: Calls `nlp.solver_get_info(...)` to provide problem dimensions and sparsity information to IPOPT.

* **`get_bounds_info(...)`**: Calls `nlp.solver_get_bounds(...)` to provide variable and constraint bounds to IPOPT.

* **`get_starting_point(...)`**: Calls `nlp.solver_get_initial_guess(...)` to provide initial variable and multiplier guesses to IPOPT.

* **`get_scaling_parameters(...)`**: This method is called by IPOPT to query scaling information. Does nothing, since scaling is done by `NLP` anyway.

* **`eval_f(...)`**: Calls `nlp.solver_eval_f(...)` to evaluate the objective function.

* **`eval_grad_f(...)`**: Calls `nlp.solver_eval_grad_f(...)` to evaluate the objective gradient.

* **`eval_g(...)`**: Calls `nlp.solver_eval_g(...)` to evaluate the constraint functions.

* **`eval_jac_g(...)`**: Calls `nlp.solver_eval_jac(...)` to evaluate the Jacobian of the constraints.

* **`eval_h(...)`**: Calls `nlp.solver_eval_hes(...)` to evaluate the Hessian of the Lagrangian.

* **`finalize_solution(...)`**: Calls `nlp.solver_finalize_solution(...)` to pass the optimal solution back to the `NLP` problem for post-processing.

* **`intermediate_callback(...)`**: Allows for custom actions during the optimization process (e.g., logging, monitoring progress).

This adapter design ensures a clean separation of concerns: the `NLP::NLP` defines the problem, `NLP::Scaling` handles transformations, and `IpoptAdapter` translates between our generic NLP interface and IPOPT's specific requirements.

## 2. NLP Scaling

Numerical optimization problems, especially those arising from dynamic optimization, often benefit significantly from scaling. Scaling transforms the problem variables and functions so that they have similar magnitudes, which can improve the numerical stability and convergence rate of the solver.

### 2.1. Scaling Interface (`NLP::Scaling`)

The `NLP` framework provides a generic scaling mechanism through the `NLP::Scaling` abstract base class.

* **`NLP::get_scaling()`**: This virtual method allows the user to provide a custom `Scaling` object. If overridden, it should return a `std::shared_ptr<Scaling>` pointing to an instance of a class derived from `Scaling`. If not overridden or if it returns `nullptr`, the framework defaults to `NoScaling`, meaning no scaling is applied.

* **`Scaling` Base Class**: This abstract base class defines the interface for various scaling operations:

  * `inplace_scale_x(f64* x_unscaled)`: Scales primal variables in place.

  * `inplace_scale_g(f64* g_unscaled)`: Scales constraint values in place.

  * `unscale_x(const f64* x_scaled, f64* x_unscaled, int number_vars)`: Unscales primal variables.

  * `unscale_g(const f64* g_scaled, f64* g_unscaled, int number_constraints)`: Unscales constraint values.

  * `unscale_f(const f64* f_scaled, f64* f_unscaled)`: Unscales objective function value.

  * `scale_f(const f64* f_unscaled, f64* f_scaled)`: Scales objective function value.

  * `scale_g(const f64* g_unscaled, f64* g_scaled, int number_constraints)`: Scales constraint values.

  * `scale_grad_f(const f64* grad_unscaled, f64* grad_scaled, int number_vars)`: Scales objective gradient.

  * `scale_jac(const f64* jac_unscaled, f64* jac_scaled, int* i_row_jac, int* j_col_jac, int jac_nnz)`: Scales Jacobian elements.

  * `scale_hes(const f64* hes_unscaled, f64* hes_scaled, int* i_row_hes, int* j_col_hes, int hes_nnz)`: Scales Hessian elements.

These methods are designed to handle the transformation of variables, functions, and their derivatives between the unscaled (user-defined) and scaled (solver-facing) domains.

### 2.2. Scaling Implementations

* **`NoScaling`**: A concrete implementation of `Scaling` that performs no actual scaling (i.e., all scaling factors are 1.0). It primarily involves `memcpy` operations to copy data without modification. This is useful when the NLP is already well-scaled.

* **`NominalScaling`**: A more advanced concrete implementation of `Scaling` that applies scaling based on nominal values provided by the user. It calculates scaling factors as `1 / nominal_value` for variables and functions. It also provides methods to create specific scaling factors for gradients, Jacobians, and Hessians based on the nominal values of the corresponding variables and functions.

* **`Other`**: Since `Scaling` is completely generic, the user can provide an own implementation, which will be automatically called by **MOO**.

### 2.3. Automatic Scaling in NLP Solver API

The `NLP` class's `solver_` methods automatically handle the application of these scaling operations:

* **Input to User Callbacks**: Variables (`x`) passed to user `eval_` callbacks are *unscaled* from the solver's internal representation. Dual variables (`lambda`) and objective factor (`sigma_f`) are *scaled* before being passed to `eval_hes`.

* **Output from User Callbacks**: Results from user `eval_` callbacks (objective value, gradient, constraints, Jacobian, Hessian) are *scaled* before being returned to the solver.

* **Bounds and Initial Guess**: Variable and constraint bounds, as well as the initial guess for primal variables and dual multipliers of variable bounds, are *scaled* when provided to the solver. The initial guess for dual variables (`lambda`) is *unscaled* for the solver.

* **Final Solution**: The optimal solution components received from the solver are *unscaled* before being passed to `finalize_solution`.

This design ensures that the user's problem definition can focus purely on the unscaled mathematical model, while the `NLP` base class transparently manages the scaling interactions with the solver.

## 3. General Dynamic Optimization Problem (GDOP) Implementation

This section details how General Dynamic Optimization Problems (GDOPs) are implemented within the **MOO** framework, leveraging the generic NLP and solver interfaces.

### 3.1. GDOP Problem Formulation

The core dynamic optimization problem addressed by the GDOP framework is:

$$
\min_{u(t), p} M(x_0, x_f, u_f, p) + \int_{t_0}^{t_f} L(x, u, p, t) \, dt
$$


subject to:

$$
\begin{aligned}
    \frac{dx}{dt} &= f(x, u, p, t),  & \forall t &\in [t_0, t_f] \\
    g^{LB} &\leq g(x, u, p, t) \leq g^{UB},  & \forall t &\in [t_0, t_f] \\
    r^{LB} &\leq r(x_0, x_f, u_f, p) \leq r^{UB}
\end{aligned}
$$

Additionally, time-invariant box constraints are included, which can also be interpreted as path constraints.

# Pontryagin’s Minimum Principle for Optimal Control with Parameter Optimization

## Problem Statement

Consider the optimal control problem

$$
\begin{aligned}
&\min_{u(\cdot),\, p} && M(x(t_f), p) \\
&\text{subject to:} \\
& && \dot{x}(t) = f(x(t), u(t), p), \quad t \in [t_0, t_f] \\
& && h(x(t), u(t), p) \leq 0, \quad \forall t \in [t_0, t_f] \\
& && r(x(t_0), x(t_f), p) \leq 0
\end{aligned}
$$

where

- $x(t) \in \mathbb{R}^n$: state vector
- $u(t) \in \mathbb{R}^m$: control vector (measurable functions)
- $p \in \mathbb{R}^q$: parameters (decision variables)
- $M: \mathbb{R}^n \times \mathbb{R}^q \to \mathbb{R}$: terminal cost function
- $f: \mathbb{R}^n \times \mathbb{R}^m \times \mathbb{R}^q \to \mathbb{R}^n$: dynamics
- $h: \mathbb{R}^n \times \mathbb{R}^m \times \mathbb{R}^q \to \mathbb{R}^r$: path constraints
- $r: \mathbb{R}^n \times \mathbb{R}^n \times \mathbb{R}^q \to \mathbb{R}^s$: boundary constraints

---

## Pontryagin’s Minimum Principle (PMP)

### Hamiltonian

Define the Hamiltonian function

$$
H(x, u, \lambda, p) = \lambda^\top f(x, u, p)
$$

where $\lambda(t) \in \mathbb{R}^n$ is the costate (adjoint) vector.

---

### Necessary Conditions for Optimality

For an optimal solution $(x^*(\cdot), u^*(\cdot), p^*)$, there exist multipliers $\lambda(t)$, $\mu(t) \geq 0$ (for path constraints), and $\nu \geq 0$ (for boundary constraints) such that:

1. **State dynamics**

$$
\dot{x}(t) = \frac{\partial H}{\partial \lambda} = f(x(t), u(t), p)
$$

2. **Costate (adjoint) dynamics**

$$
\dot{\lambda}(t) = - \frac{\partial H}{\partial x} - \left(\frac{\partial h}{\partial x}\right)^\top \mu(t) = - \left( \frac{\partial f}{\partial x}(x,u,p) \right)^\top \lambda(t) - \left( \frac{\partial h}{\partial x}(x,u,p) \right)^\top \mu(t)
$$

3. **Hamiltonian minimization condition**

At each time $t \in [t_0, t_f]$, the control minimizes the Hamiltonian subject to path constraints:

$$
u(t) \in \arg\min_{v} \left\{ H(x(t), v, \lambda(t), p) \mid h(x(t), v, p) \leq 0 \right\}
$$

4. **Parameter stationarity condition**

Parameters $p$ must satisfy

$$
0 = \frac{\partial M}{\partial p}(x_f, p) + \left(\frac{\partial r}{\partial p}(x_0, x_f, p)\right)^\top \nu + \int_{t_0}^{t_f} \left[ \lambda(t)^\top \frac{\partial f}{\partial p}(x,u,p) + \mu(t)^\top \frac{\partial h}{\partial p}(x,u,p) \right] dt
$$

5. **Transversality (boundary) conditions**

$$
\lambda(t_f) = \frac{\partial M}{\partial x_f}(x_f, p) + \left(\frac{\partial r}{\partial x_f}(x_0, x_f, p)\right)^\top \nu
$$

$$
\lambda(t_0) = -\left(\frac{\partial r}{\partial x_0}(x_0, x_f, p)\right)^\top \nu
$$

6. **Complementarity slackness for path constraints**

For each component $i$ of the path constraints,

$$
\mu_i(t) \geq 0, \quad h_i(x(t), u(t), p) \leq 0, \quad \mu_i(t) h_i(x(t), u(t), p) = 0
$$

7. **Complementarity slackness for boundary constraints**

For each component $j$ of boundary constraints,

$$
\nu_j \geq 0, \quad r_j(x_0, x_f, p) \leq 0, \quad \nu_j r_j(x_0, x_f, p) = 0
$$

---

## Derivation Sketch

1. **Form Lagrangian**

Introduce multipliers $\lambda(t)$ for dynamics, $\mu(t)$ for path constraints, $\nu$ for boundary constraints, and form

$$
\begin{aligned}
\mathcal{L} = &\, M(x(t_f), p) + \nu^\top r(x(t_0), x(t_f), p) \\
& + \int_{t_0}^{t_f} \lambda(t)^\top [f(x,u,p) - \dot{x}(t)] dt + \int_{t_0}^{t_f} \mu(t)^\top h(x,u,p) dt
\end{aligned}
$$

2. **Perform first variation $\delta \mathcal{L} = 0$** with respect to $x(\cdot)$, $u(\cdot)$, and $p$.

3. **Apply integration by parts** on terms with $\dot{x}$ to move derivatives from $\delta x$ to $\lambda$.

4. **Obtain adjoint equation** from variations in $x$, leading to

$$
\dot{\lambda} = - \frac{\partial f}{\partial x}^\top \lambda - \frac{\partial h}{\partial x}^\top \mu
$$

with boundary conditions involving $\lambda(t_0), \lambda(t_f)$.

5. **Optimality in control** arises from stationarity w.r.t. $u$:

$$
0 = \frac{\partial}{\partial u} \left[ \lambda^\top f + \mu^\top h \right]
$$

subject to $h(x,u,p) \leq 0$.

6. **Optimality in parameters** $p$ follows similarly:

$$
0 = \frac{\partial M}{\partial p} + \left(\frac{\partial r}{\partial p}\right)^\top \nu + \int \left[ \lambda^\top \frac{\partial f}{\partial p} + \mu^\top \frac{\partial h}{\partial p} \right] dt
$$

7. **Complementarity and feasibility conditions** on $\mu(t)$, $\nu$, $h$, and $r$ come from the Karush-Kuhn-Tucker (KKT) conditions on inequality constraints.

---

## Remarks

- The PMP generalizes the Euler-Lagrange equations in calculus of variations to optimal control problems with constraints.
- The inclusion of parameters $p$ as decision variables adds integral terms in the stationarity conditions.
- Path and boundary constraints introduce multipliers $\mu$, $\nu$ with complementarity conditions.
- In practice, solving this system involves two-point boundary value problems plus algebraic complementarity conditions, often handled numerically.

---

If you want, I can also provide a **LaTeX source file** or an extended explanation of each step.

---

**End of document.**


### 3.2. Discretization

The problem is discretized using a **Radau IIA collocation scheme** of arbitrary order (1–199) with 1-100 stages.
The collocation scheme can differ from interval to interval. It is therefore possible to implement a hp-adaptive pseudospectral mesh refinement algorithm.

The `Mesh` and `fLGR` objects are fundamental to this discretization:

* **`fLGR`**: Defines the properties of the collocation method itself, such as the type (Radau IIA) and potentially the collocation points and weights.

* **`Mesh`**: Defines the structure of the time horizon, dividing it into intervals and specifying the number of collocation stages within each interval. The `hp-adaptive` capability implies that the order of the collocation polynomial (stages) and the number/length of intervals can be dynamically adjusted during the optimization process to improve efficiency and accuracy.

### 3.3. Problem Definition (User Implementation: `GDOP::Problem`, `FullSweep`, `BoundarySweep`)

To define a specific GDOP, you will need to implement concrete classes derived from `GDOP::FullSweep` and `GDOP::BoundarySweep`. This is the primary area where the user's specific problem dynamics, objectives, and constraints are encoded. These user-defined components are then aggregated into a `GDOP::Problem` object.

#### 3.3.1. Overview of Key Components

1. **`ProblemConstants`**:
   This struct encapsulates all constant parameters of the GDOP, including:

   * Sizes of state (`x_size`), control (`u_size`), and parameter (`p_size`) vectors.

   * Flags indicating the presence of Mayer (`has_mayer`) and Lagrange (`has_lagrange`) terms in the objective function.

   * Dimensions of dynamic (`f_size`), path (`g_size`), and boundary (`r_size`) constraints.

   * Index ranges for accessing `f`, `g`, and `r` terms within evaluation buffers.

   * Bounds for variables (`x_bounds`, `u_bounds`, `p_bounds`) and constraints (`r_bounds`, `g_bounds`).

   * Fixed initial (`x0_fixed`) and final (`xf_fixed`) state components.

   * References to the `Mesh` and `fLGR` objects, which define the discretization of the problem.

2. **`FullSweep` (Abstract Base Class)**:
   This class represents the "full sweep" of the dynamic optimization problem, encompassing
   the objective and constraint functions that apply across the entire time horizon
   (i.e., at each collocation point within every mesh interval).
   You **must derive from this class** and implement the following pure virtual methods:

   * `virtual void callback_eval(const f64* xu_nlp, const f64* p) = 0;`
     Calculates the values of the Lagrange term (L), dynamic equations (f), and path constraints (g).
     These values should be written into `buffers.eval`.

   * `virtual void callback_jac(const f64* xu_nlp, const f64* p) = 0;`
     Calculates the Jacobian (first derivatives) of L, f, and g with respect to `x`, `u`, and `p`.
     These values should be written into `buffers.jac`.

   * `virtual void callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, f64* lambda) = 0;`
     Calculates the augmented Hessian (second derivatives) of L, f, and g. This involves
     the second derivatives of the Lagrangian function (objective + constraints scaled by multipliers).
     The results should be written into `buffers.aug_hes` and `buffers.aug_pp_hes`.
     `lagrange_factors` are the exact factors for the Lagrange terms, and `lambda` are the exact multipliers
     for the constraints (f and g).

3. **`BoundarySweep` (Abstract Base Class)**:
   This class represents the "boundary sweep" of the dynamic optimization problem, dealing with
   objective terms and constraints applied only at the initial (`x(t0)`) and final (`x(tf)`)
   time points.
   You **must derive from this class** and implement the following pure virtual methods:

   * `virtual void callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) = 0;`
     Calculates the values of the Mayer term (M) and boundary constraints (r).
     These values should be written into `buffers.eval`.

   * `virtual void callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p) = 0;`
     Calculates the Jacobian of M and r with respect to `x0`, `xf`, `uf`, and `p`.
     These values should be written into `buffers.jac`.

   * `virtual void callback_aug_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, f64* lambda) = 0;`
     Calculates the augmented Hessian of M and r.
     The results should be written into `buffers.aug_hes`.
     `mayer_factor` is the exact multiplier for the Mayer term, and `lambda` are the exact multipliers
     for the boundary constraints (r).

#### 3.3.2. Data Structures for Sparsity and Buffers

The sparsity structures (e.g., `JacobianLFG::dx`, `AugmentedHessianLFG::dx_dx`) hold
`JacobianSparsity` or `HessianSparsity` objects. These objects contain a `buf_index`
which is crucial for writing to the `FullSweepBuffers` and `BoundarySweepBuffers`.

**One Function Per `FunctionLFG` / `FunctionMR`**:
Each `FunctionLFG` object in the `FullSweep::lfg` vector corresponds to a single
scalar function (e.g., the Lagrange objective term L, or one component of the
dynamic equations f, or one component of the path constraints g).
Similarly, each `FunctionMR` object in the `BoundarySweep::mr` vector corresponds
to a single scalar function (e.g., the Mayer objective term M, or one component
of the boundary constraints r).

**Important Note on `buf_index` and Buffer Mapping**:
The sparsity information (e.g., `FullSweep::lfg`, `FullSweep::aug_hes`, etc.)
defines the structure and *sorted order* of the non-zero elements within the
evaluation, Jacobian, and Hessian buffers.

* The `buf_index` provided within each `FunctionLFG` or `FunctionMR` struct,
  and within each `JacobianSparsity` or `HessianSparsity` entry, indicates the
  *precise linear index* in the corresponding `FixedVector<f64>` buffer where the
  calculated value should be written.

* **For `FullSweep` (L, f, g) calculations (per collocation node)**:

  * The `buffers.eval`, `buffers.jac`, and `buffers.aug_hes` are allocated as
    large contiguous blocks of memory, designed to hold data for *all*
    collocation nodes.

  * The `eval_size`, `jac_size`, `aug_hes_size` members of `FullSweepBuffers`
    represent the memory chunk size required for *a single collocation node*.

  * To correctly write results for a specific node `(i, j)` (interval `i`, node `j`),
    you must first calculate the starting address for that node's data chunk
    within the larger buffer. The inline helper methods `get_eval_buffer(i, j)`,
    `get_jac_buffer(i, j)`, and `get_aug_hes_buffer(i, j)` provide this base pointer
    for the current node's data chunk.

  * You then use the `buf_index` from the relevant `FunctionLFG` object or
    `JacobianSparsity`/`HessianSparsity` entry as an *offset from this base pointer*:
    `base_pointer_for_node + sparsity_entry.buf_index`.

* **For `BoundarySweep` (M, r) calculations (single evaluation)**:

  * The `buffers.eval`, `buffers.jac`, and `buffers.aug_hes` are single, fixed-size buffers.
    Boundary conditions are evaluated once for the problem, not per collocation node.

  * Therefore, the `buf_index` from `FunctionMR` or the respective `JacobianSparsity`/`HessianSparsity`
    entry directly gives the index into these single buffers (no additional node-based offset is needed).

**Hessian Structure (Lower Triangular)**:
For Hessian calculations (`callback_aug_hes`), you are expected to provide only the
**lower triangular part** of the symmetric Hessian matrices. This means that for any
`HessianSparsity` entry `(row, col)`, you should only write a value if `row >= col`.
The `HessianSparsity` structures will only contain entries for the lower triangular part.

#### 3.3.3. Mathematical Expectations for Callbacks

The callbacks expect the following mathematical quantities to be computed and written into the provided buffers:

##### `FullSweep` Callbacks:

* **`callback_eval(const f64* xu_nlp, const f64* p)`**:
    This callback is responsible for computing the scalar values of the Lagrange term ($L$), each dynamic equation ($f_i$), and each path constraint ($g_j$) at a given collocation point. These computed values should be written into the `buffers.eval` array using the `buf_index` provided by the corresponding `FunctionLFG` objects (e.g., `buffers.eval[lfg[idx].buf_index] = value;`). The order of calculation within this callback does not affect correctness, as the `buf_index` ensures proper placement.

* **`callback_jac(const f64* xu_nlp, const f64* p)`**:
    This callback computes the partial derivatives (Jacobian) of the Lagrange term ($L$), dynamic equations ($f$), and path constraints ($g$) with respect to $x$, $u$, and $p$ at a given collocation point. These derivative values should be written into the `buffers.jac` array, using the `buf_index` specified by the `JacobianSparsity` entries within each `FunctionLFG` object (e.g., `buffers.jac[lfg[idx].jac.dx[entry_idx].buf_index] = derivative_value;`).

* **`callback_aug_hes(const f64* xu_nlp, const f64* p, const FixedField<f64, 2>& lagrange_factors, f64* lambda)`**:
    This callback computes the lower triangular part of the augmented Hessian (second derivatives of the Lagrangian) with respect to $(x, u, p)$ at a given collocation point. The Lagrangian $\mathcal{L}$ for a single node is:

    $$
    \mathcal{L}(x, u, p, t, \lambda, \text{$lfactor$}) = \text{$lfactor$} \cdot L(x, u, p, t) + \sum_{k} \lambda_k \cdot f_k(x, u, p, t) + \sum_{k} \lambda_k \cdot g_k(x, u, p, t)
    $$

    The computed second derivative values should be written into the `buffers.aug_hes` array, using the `buf_index` specified by the `HessianSparsity` entries. For parameters, values are written to `buffers.aug_pp_hes`.

##### `BoundarySweep` Callbacks:

* **`callback_eval(const f64* x0_nlp, const f64* xuf_nlp, const f64* p)`**:
    This callback computes the scalar values of the Mayer term ($M$) and each boundary constraint ($r_i$). These values should be written into the `buffers.eval` array using the `buf_index` provided by the corresponding `FunctionMR` objects.

* **`callback_jac(const f64* x0_nlp, const f64* xuf_nlp, const f64* p)`**:
    This callback computes the partial derivatives (Jacobian) of the Mayer term ($M$) and boundary constraints ($r$) with respect to $x_0$, $x_f$, $u_f$, and $p$. These values should be written into the `buffers.jac` array using the `buf_index` specified by the `JacobianSparsity` entries within each `FunctionMR` object.

* **`callback_aug_hes(const f64* x0_nlp, const f64* xuf_nlp, const f64* p, const f64 mayer_factor, f64* lambda)`**:
    This callback computes the lower triangular part of the augmented Hessian (second derivatives of the boundary Lagrangian) with respect to $(x_0, x_f, u_f, p)$. The boundary Lagrangian $\mathcal{L}_B$ is:

    $$
    \mathcal{L}_B(x_0, x_f, u_f, p, \lambda, \lambda_M) = \lambda_M \cdot M(x_0, x_f, u_f, p) + \sum_{k} \lambda_k \cdot r_k(x_0, x_f, u_f, p)
    $$

    Here $\lambda_M$ denotes the `mayer_factor`. The computed second derivative values should be written into the `buffers.aug_hes` array using the `buf_index` specified by the `HessianSparsity` entries.

### 3.4. GDOP-Specific NLP Implementation (`GDOP::GDOP`)

The `GDOP::GDOP` class is the concrete implementation of the abstract `NLP::NLP` interface for a General Dynamic Optimization Problem. It acts as the bridge between the problem's continuous-time formulation (discretized by `Mesh` and `fLGR`) and the generic NLP solver.

`GDOP::GDOP` inherits from `NLP::NLP` and implements all the pure virtual methods defined in the `NLP` interface. Its primary responsibility is to translate the requests from the NLP solver (e.g., "give me the objective value at this point") into calls to the user-defined `GDOP::FullSweep` and `GDOP::BoundarySweep` methods, aggregating the results across all collocation points and intervals to form the complete NLP objective, constraints, and derivatives.

#### 3.4.1. Internal Structure and Variable/Constraint Layout

The `GDOP::GDOP` class manages the mapping of the continuous GDOP variables and constraints onto a flattened NLP vector. This involves careful calculation of offsets and accumulated indices:

* **`off_x`, `off_u`, `off_p`**: These represent the sizes of the state, control, and parameter vectors, respectively, as defined in `ProblemConstants`.

* **`off_xu`**: The sum of `off_x` and `off_u`, representing the number of variables (state and control) at a single collocation point.

* **`off_acc_xu`**: A `FixedField<int, 2>` that stores the starting index in the flattened NLP variable vector for the `(x, u)` block at each collocation node `(i, j)`. This allows efficient access to $x_{i,j}$ and $u_{i,j}$.

* **`off_last_xu`**: The starting index of the `(x, u)` block corresponding to the last collocation point in the NLP variable vector.

* **`off_xu_total`**: The total number of `(x, u)` variables across all collocation nodes, which also serves as the starting index for the parameter variables `p` in the NLP variable vector.

* **`off_fg_total`**: The total number of `(f, g)` constraints across all collocation nodes, which also serves as the starting index for the boundary constraints `r` in the NLP constraint vector.

* **`off_acc_fg`**: A `FixedField<int, 2>` that stores the starting index in the flattened NLP constraint vector for the `(f, g)` block at each collocation node `(i, j)`.

* **`off_acc_jac_fg`**: A `FixedVector<int>` that stores the accumulated non-zero counts for the Jacobian of `(f, g)` constraints up to a certain point, aiding in populating the Jacobian sparsity pattern efficiently.

These offsets are crucial for `GDOP::GDOP` to correctly interpret and populate the flattened NLP vectors for variables (`curr_x`), constraints (`curr_g`), and their derivatives (`curr_grad`, `curr_jac`, `curr_hes`).

#### 3.4.2. Implementation of `NLP::NLP` Virtual Methods

`GDOP::GDOP` overrides the pure virtual methods from `NLP::NLP` to provide the specific GDOP logic:

* **`get_sizes(int& number_vars, int& number_constraints)`**:
  This method calculates the total number of NLP variables (sum of all `x`, `u` at all nodes, plus `p`) and constraints (sum of all `f`, `g` at all nodes, plus `r`). It also initializes all the internal offset members (`off_x`, `off_u`, `off_p`, `off_xu`, `off_last_xu`, `off_xu_total`, `off_fg_total`, `off_acc_xu`, `off_acc_fg`).

* **`get_bounds(FixedVector<f64>& x_lb, FixedVector<f64>& x_ub, FixedVector<f64>& g_lb, FixedVector<f64>& g_ub)`**:
  This method populates the NLP variable and constraint bounds. It maps the `problem.pc`'s `x_bounds`, `u_bounds`, `p_bounds`, `g_bounds`, and `r_bounds` to the flattened NLP bounds. Special handling is included for fixed initial (`x0_fixed`) and final (`xf_fixed`) states, where their bounds are set to the fixed values.

* **`get_initial_guess(...)`**:
  This method provides an initial guess for the NLP variables and duals. It utilizes an internal `std::unique_ptr<PrimalDualTrajectory> initial_guess`. If the `initial_guess` is not compatible with the current `Mesh` and `fLGR`, it is automatically interpolated to match the current discretization. The primal variables are flattened into `x_init`, and the costates (dual variables) and duals for bounds are transformed and flattened into `lambda_init`, `z_lb_init`, and `z_ub_init`.

* **`get_nnz(int& nnz_jac, int& nnz_hes)`**:
  This method calculates the total number of non-zero elements in the NLP Jacobian and Hessian. It calls `init_jacobian_nonzeros` and `init_hessian_nonzeros` internally, which analyze the sparsity patterns defined in `GDOP::Problem`'s `FullSweep` and `BoundarySweep` components and aggregate them across the mesh.

* **`get_jac_sparsity(FixedVector<int>& i_row_jac, FixedVector<int>& j_col_jac)`**:
  This complex method constructs the full NLP Jacobian sparsity pattern (row and column indices in COO format). It iterates through each collocation node and boundary point, combining:
  * The constant derivative terms from the collocation matrix (`collocation.D`) for the dynamic equations.
  * The sparsity patterns from the user's `FullSweep` (for `f` and `g` functions) and `BoundarySweep` (for `r` functions) implementations.
  It carefully maps these local sparsity patterns to the global NLP Jacobian indices using the calculated offsets. It also populates `const_der_jac` for the constant parts of the Jacobian.

* **`get_hes_sparsity(FixedVector<int>& i_row_hes, FixedVector<int>& j_col_hes)`**:
  This method constructs the full NLP Hessian sparsity pattern (row and column indices in COO format). It leverages internal `BlockSparsity` objects (`hes_a` through `hes_h`) and `OrderedIndexSet`s (`A` through `H`) to efficiently manage the aggregation of sparsity from the `FullSweep` and `BoundarySweep` Hessian structures across the entire mesh. The total `nnz_hes` is computed based on these aggregated block sparsity patterns. It also allocates `lagrange_obj_factors` which are used in Hessian evaluations.

* **`eval_f(...)`, `eval_g(...)`, `eval_grad_f(...)`, `eval_jac_g(...)`, `eval_hes(...)`**:
  These methods serve as wrappers. They first call internal `callback_evaluation`, `callback_jacobian`, or `callback_hessian` to trigger the user's `FullSweep` and `BoundarySweep` callbacks, filling the problem-level buffers. Then, they call their respective `_internal` methods (e.g., `eval_f_internal`) to aggregate these results from the problem's buffers into the flattened NLP buffers expected by the solver. For `eval_hes`, it also calls `update_curr_lambda_obj_factors` to prepare the Lagrange multipliers and objective factors for the user's Hessian callback.

* **`finalize_solution(...)`**:
  After the solver finds an optimal solution, this method is called. It takes the optimal NLP primal variables (`opt_x`) and duals (`opt_lambda`, `opt_z_lb`, `opt_z_ub`), and transforms them back into `Trajectory` and `CostateTrajectory` objects, which are more convenient for analysis and post-processing in the context of dynamic optimization. These are stored in `optimal_solution`.

#### 3.4.3. Hessian Sparsity Layout

The Hessian of the Lagrangian for the entire NLP problem has a specific block structure due to the collocation scheme. Only the lower triangular part of this symmetric matrix is provided. The layout can be visualized as:

$$
\begin{array}{c|c|c|c|c|c|c|c}
\text{Variables} & x_{00} & x_{01} u_{01} & x_{02} u_{02} & \dots & x_{nm-1} u_{nm-1} & x_{nm} u_{nm} & p \\
\hline
x_{00} & A & & & & & & \\
\hline
x_{01} & & B & & & & & \\
u_{01} & & & & & & & \\
\hline
x_{02} & & & B & & & & \\
u_{02} & & & & & & & \\
\hline
\vdots & & & & \ddots & & & \\
\hline
x_{nm-1} & & & & & B & & \\
u_{nm-1} & & & & & & & \\
\hline
x_{nm} & C & & & & & D & \\
u_{nm} & & & & & & & \\
\hline
p & E & F & F & \dots & F & G & H \\
\end{array}
$$

Where:

* $x_{00}$: Refers to the initial state variables $x(t_0)$.

* $x_{ij} u_{ij}$: Refers to state and control variables at collocation node $j$ in interval $i$.

* $x_{nm} u_{nm}$: Refers to state and control variables at the final collocation node $m$ in the last interval $n$.

* $p$: Refers to the global parameters.

The blocks A-H represent sub-matrices derived from the problem's `AugmentedHessianMR` (boundary) and `AugmentedHessianLFG` (full sweep) structures:

* **A**: Lower triangular block for $x_0$ derivatives, derived from `problem.boundary->aug_hes->dx0_dx0`.

* **B**: Lower triangular block for $(x, u)$ derivatives at an internal collocation node, derived from `problem.full->aug_hes->dx_dx`, `du_dx`, `du_du`. This block is repeated for each internal interval.

* **C**: Rectangular block for derivatives of $x_{\text{final}}$ and $u_{\text{final}}$ with respect to $x_0$, derived from `problem.boundary->aug_hes->dxf_dx0`, `duf_dx0`.

* **D**: Lower triangular block for derivatives of $x_{\text{final}}$ and $u_{\text{final}}$ with respect to themselves, derived from `problem.boundary->aug_hes->dxf_dxf`, `duf_dxf`, `duf_duf`.

* **E**: Rectangular block for derivatives of parameters $p$ with respect to $x_0$, derived from `problem.boundary->aug_hes->dp_dx0`.

* **F**: Rectangular block for derivatives of parameters $p$ with respect to $(x, u)$ at an internal collocation node, derived from `problem.full->aug_hes->dp_dx`, `dp_du`. This block is repeated for each internal interval.

* **G**: Rectangular block for derivatives of parameters $p$ with respect to $x_{\text{final}}$ and $u_{\text{final}}$, derived from `problem.boundary->aug_hes->dp_dxf`, `dp_duf`.

* **H**: Lower triangular block for derivatives of parameters $p$ with respect to themselves, derived from `problem.boundary->aug_hes->dp_dp` and `problem.full->aug_pp_hes->dp_dp`.

This block structure is essential for efficiently populating the Hessian matrix by iterating through the relevant sparsity patterns and applying the correct offsets.

### 3.5. GDOP Orchestration and Strategies

The GDOP framework provides sophisticated mechanisms for managing the overall optimization process, especially for problems that benefit from adaptive mesh refinement. This is primarily handled by the `GDOP::Orchestrator` and its concrete implementation, `GDOP::MeshRefinementOrchestrator`, in conjunction with various pluggable `GDOP::Strategies`.

#### 3.5.1. `GDOP::Orchestrator` (Base Class)

The `GDOP::Orchestrator` serves as the abstract base class for managing the entire GDOP solution process. It holds references to the `GDOP` problem instance, a collection of `Strategies`, and the chosen `NLP::NLPSolver`. Its core responsibility is to define the high-level workflow for solving the dynamic optimization problem.

* **`GDOP& gdop`**: A reference to the concrete `GDOP` NLP problem instance (which itself implements the `NLP::NLP` interface).

* **`std::unique_ptr<Strategies> strategies`**: A unique pointer to an object aggregating various pluggable strategy components.

* **`NLP::NLPSolver& solver`**: A reference to the chosen NLP solver (e.g., `IpoptSolver`).

* **`virtual void optimize() = 0;`**: The pure virtual method that concrete orchestrators must implement to define the optimization flow.

#### 3.5.2. `GDOP::MeshRefinementOrchestrator` (Adaptive Solving)

The `GDOP::MeshRefinementOrchestrator` is the concrete implementation of `GDOP::Orchestrator` that specifically handles adaptive mesh refinement. It drives an iterative process of solving the NLP on a given mesh, analyzing the solution, adapting the mesh if necessary, and then re-solving on the refined mesh. This iterative loop continues until convergence or other termination criteria are met.

The `optimize()` method in `MeshRefinementOrchestrator` outlines this iterative process:

1. **Reset Strategies**: `strategies->reset(gdop);`
   Initializes or resets any internal state within the strategies (e.g., iteration counters for mesh refinement).

2. **Initial Guess Generation**: `auto initial_guess = strategies->get_initial_guess(gdop);`
   Obtains the initial guess for the primal variables for the first NLP solve. This is provided by the `Initialization` strategy.

3. **Scaling Setup**: `gdop.set_scaling_factory(strategies->scaling_factory);`
   Informs the `GDOP` NLP instance about the chosen `ScalingFactory` strategy, which will be used to create the appropriate `NLP::Scaling` object for the solver.

4. **Main Optimization Loop (`for(;;)`)**:

   * **Set Initial Guess**: `gdop.set_initial_guess(std::move(initial_guess));`
     Provides the current initial guess (either the very first one or an interpolated one from a previous iteration) to the `GDOP` NLP problem.

   * **Solve NLP**: `solver.optimize();`
     Invokes the underlying NLP solver (e.g., IPOPT) to solve the current NLP problem on the current mesh.

   * **Mesh Refinement Detection**: `auto mesh_update = strategies->detect(gdop.get_mesh(), gdop.get_collocation(), *gdop.get_optimal_solution());`
     Calls the `MeshRefinement` strategy to analyze the optimal solution from the current solve and determine if a mesh update is needed. If no update is proposed (`!mesh_update`), the loop breaks.

   * **Create Refined Mesh**: `auto refined_mesh = Mesh(std::move(mesh_update), gdop.get_collocation());`
     If a mesh update is detected, a new `Mesh` object is constructed based on the proposed changes.

   * **Interpolate for Warm Start**: `initial_guess = strategies->get_refined_initial_guess(gdop.get_mesh(), refined_mesh, gdop.get_collocation(), *gdop.get_optimal_solution());`
     The `RefinedInitialization` strategy is used to interpolate the optimal solution from the *old* mesh onto the *new, refined* mesh. This interpolated solution serves as a "warm start" for the next NLP solve, significantly improving convergence. The `NLP::Option::WarmStart` is enabled for the solver.

   * **Update GDOP Instance**: `gdop.update(std::move(refined_mesh));`
     The `GDOP` NLP problem instance is updated with the new mesh, causing its internal buffers and sparsity patterns to be re-allocated and re-initialized for the next iteration.

5. **Post-Optimization Actions**:

   * **Verification**: `strategies->verify(gdop, *gdop.get_optimal_solution());`
     Calls the `Verifier` strategy to perform post-optimization checks on the final optimal solution (e.g., by simulating the system with the optimal controls and comparing states).

   * **Emission**: `strategies->emit(*gdop.get_optimal_solution()->primals);`
     Calls the `Emitter` strategy to output the final optimal primal trajectory (e.g., to a CSV file).

   * **Costate Output**: The optimal costates and bound duals are also emitted to CSV files for detailed analysis.

This orchestration pattern allows for highly flexible and adaptive solution procedures, crucial for complex dynamic optimization problems.

#### 3.5.3. `GDOP::Strategies` (Pluggable Behaviors)

The `GDOP::Strategies` class is a central aggregation point for various pluggable behaviors that define how the `Orchestrator` operates. It holds `std::shared_ptr` to different abstract strategy interfaces, allowing users to inject custom implementations for each aspect of the optimization process. This design promotes modularity and reusability.

Each strategy interface defines a specific operation:

* **`Initialization`**:

  * **Purpose**: Creates the initial guess for primal variables for the *very first* NLP solve.

  * **`operator()(const GDOP& gdop)`**: Returns a `std::unique_ptr<PrimalDualTrajectory>`.

  * **Example Implementations**: `ConstantInitialization` (simple constant guess), `SimulationInitialization` (initial guess from a forward simulation).

* **`RefinedInitialization`**:

  * **Purpose**: Creates an initial guess for variables on a *new, refined mesh*, typically by interpolating a previously found optimal solution from an older mesh. This is crucial for warm-starting the NLP solver during mesh refinement.

  * **`operator()(const Mesh& old_mesh, const Mesh& new_mesh, const fLGR& collocation, const PrimalDualTrajectory& trajectory)`**: Returns a `std::unique_ptr<PrimalDualTrajectory>`.

  * **Example Implementations**: `InterpolationRefinedInitialization` (uses an `Interpolation` strategy to create the refined guess).

* **`Simulation`**:

  * **Purpose**: Simulates the full system over the entire time horizon given a control trajectory. Used for verification or generating initial guesses.

  * **`operator()(const ControlTrajectory& controls, int num_steps, f64 start_time, f64 stop_time, f64* x_start_values)`**: Returns a `std::unique_ptr<Trajectory>`.

  * **Example Implementations**: `NoSimulation` (placeholder), custom simulation logic.

* **`SimulationStep`**:

  * **Purpose**: Simulates a single step or segment of the system. Useful for finer-grained analysis or model checking.

  * **`operator()(const ControlTrajectory& controls, f64 start_time, f64 stop_time, f64* x_start_values)`**: Returns a `std::unique_ptr<Trajectory>`.

  * **Example Implementations**: `NoSimulationStep` (placeholder), custom single-step simulation.

* **`MeshRefinement`**:

  * **Purpose**: Analyzes the current optimal solution and proposes updates to the mesh structure (e.g., adding intervals, changing collocation order) to improve accuracy or efficiency.

  * **`reset(const GDOP& gdop)`**: Resets any internal state of the refinement strategy (e.g., iteration counters).

  * **`operator()(const Mesh& mesh, const fLGR& collocation, const PrimalDualTrajectory& trajectory)`**: Returns a `std::unique_ptr<MeshUpdate>` (a structure describing the proposed mesh changes) or `nullptr` if no refinement is needed.

  * **Example Implementations**: `NoMeshRefinement` (no refinement), `L2BoundaryNorm` (h-method: [Adaptively Refined Mesh for fLGR-Based Dynamic Optimization](https://www.researchgate.net/publication/390454809_Adaptively_Refined_Mesh_for_Collocation-Based_Dynamic_Optimization)).

* **`Interpolation`**:

  * **Purpose**: Interpolates a trajectory (e.g., states or controls) from one mesh to another. Used by `RefinedInitialization` and potentially other strategies.

  * **`operator()(const Mesh& old_mesh, const Mesh& new_mesh, const fLGR& collocation, const std::vector<f64>& values)`**: Returns a `std::vector<f64>` of interpolated values.

  * **Example Implementations**: `LinearInterpolation`, `PolynomialInterpolation` (uses collocation scheme for higher-order interpolation).

* **`Emitter`**:

  * **Purpose**: Handles the output of the final optimal trajectory (e.g., saving to files, logging).

  * **`operator()(const Trajectory& trajectory)`**: Returns an `int` status code.

  * **Example Implementations**: `NoEmitter` (no output), `CSVEmitter` (saves to CSV).

* **`Verifier`**:

  * **Purpose**: Verifies the optimality or quality of the final solution, typically by comparing the optimized trajectory against a high-fidelity simulation or known analytical solution.

  * **`operator()(const GDOP& gdop, const PrimalDualTrajectory& trajectory)`**: Returns a `bool` indicating success or failure.

  * **Example Implementations**: `NoVerifier` (no verification), `SimulationVerifier` (verifies by running a simulation and comparing states).

* **`ScalingFactory`**:

  * **Purpose**: Creates and provides an `NLP::Scaling` object to the `GDOP` NLP instance. This allows the user to define custom scaling strategies.

  * **`operator()(const GDOP& gdop)`**: Returns a `std::shared_ptr<NLP::Scaling>`.

  * **Example Implementations**: `NoScalingFactory` (no scaling), custom factories for `NominalScaling` or other scaling methods.

The `GDOP::Strategies` class itself provides convenience methods (e.g., `get_initial_guess`, `simulate`, `detect`) that simply forward the calls to the underlying concrete strategy objects. This makes it easy to interact with the entire set of strategies through a single interface.

#### 3.5.4. The Power of Pluggable Strategies

This strategy pattern is powerful because it allows:

* **Customization**: Users can implement their own specific initialization, refinement, simulation, or output logic without modifying the core GDOP framework.

* **Experimentation**: Different strategies can be easily swapped in and out to test their impact on solver performance, accuracy, and robustness.

* **Modularity**: Each component focuses on a single responsibility, leading to cleaner and more maintainable code.

* **Limited Interaction**: As highlighted, the orchestrator and strategies interact with the low-level `GDOP` problem primarily through high-level `Trajectory` and `PrimalDualTrajectory` objects, abstracting away the complexities of the flattened NLP representation. This makes the overall system very "nice" and user-friendly from a high-level perspective.

## 4. Usage

1. Define your specific GDOP by creating classes that inherit from `FullSweep` and `BoundarySweep`.

2. In these derived classes, implement the `callback_eval`, `callback_jac`, and `callback_aug_hes` methods
   to compute your problem's specific objective, dynamics, and constraints, along with their derivatives.
   You will use the `lfg` (`FixedVector<FunctionLFG>`) and `mr` (`FixedVector<FunctionMR>`) members,
   and their contained sparsity information, to determine the correct `buf_index` where the
   computed values must be written.

3. Instantiate `ProblemConstants` with your GDOP's dimensions and properties.

4. Create instances of your derived `FullSweep` and `BoundarySweep` classes.

5. Finally, construct a `GDOP::Problem` object using `std::unique_ptr`s to your `FullSweep`,
   `BoundarySweep`, and `ProblemConstants` instances. This `Problem` object will then be
   passed to the `GDOP::GDOP` instance, which in turn acts as the NLP for the solver.

By following this structure, you can define complex General Dynamic Optimization Problems and
leverage the GDOP framework for numerical solution.