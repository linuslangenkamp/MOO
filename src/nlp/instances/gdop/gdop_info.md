# General Dynamic Optimization Problem (GDOP)

This folder contains an NLP implementation of the following dynamic optimization problem:

$$
\min_{u(t), p} M(x_0, x_f, p) + \int_{t_0}^{t_f} L(x, u, p, t) \, dt
$$

subject to:

$$
\begin{aligned}
    \frac{dx}{dt} &= f(x, u, p, t),  & \forall t &\in [t_0, t_f] \\
    g^{LB} &\leq g(x, u, p, t) \leq g^{UB},  & \forall t &\in [t_0, t_f] \\
    r^{LB} &\leq r(x_0, x_f, p) \leq r^{UB}
\end{aligned}
$$

Additionally, time-invariant box constraints are included, which can also be interpreted as path constraints.

## Discretization

The problem is discretized using a **Radau IIA collocation scheme** of arbitrary order (1–199) with 1-100 stages.
The collocation scheme can differ from interval to interval. It is therefore possible to implement a hp-adaptive pseudospectral mesh refinement algorithm.
