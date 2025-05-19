F = [[0, 2], [0, 2], [1], [1, 2], [3]]

# Step 1: Build sparsity map: M[(v1, v2)] → list of function outputs (rows)
M = dict()
for f in range(len(F)):
    for v1 in F[f]:
        for v2 in F[f]:
            key = tuple(sorted((v1, v2)))
            if key not in M:
                M[key] = [f]
            else:
                M[key].append(f)

# Step 2: Define colors and reverse lookup
COLORS = [[0, 1], [2, 3]]

# Step 3: Build Q[(c1, c2)] → [ ([i1, i2, ...], nnz_index) ]
Q = dict()
nnz_lookup = dict()
NNZ = 0

for i in range(5):
    for j in range(i + 1):
        key = (j, i) if j <= i else (i, j)
        if key in M:
            print(key)
            nnz_lookup[key] = NNZ
            NNZ += 1

for c1 in range(len(COLORS)):
    for c2 in range(c1 + 1):  # symmetric
        pair_entries = []
        for i1 in COLORS[c1]:
            for i2 in COLORS[c2]:
                key = tuple(sorted((i1, i2)))
                if key in M:
                    row_indices = tuple(set(M[key]))
                    flat_idx = nnz_lookup[key]
                    pair_entries.append([row_indices, flat_idx])
        Q[(c1, c2)] = pair_entries

print("Q mapping (color pairs → [([rows], flat_idx)]):\n")
for key, val in Q.items():
    print(f"{key}: {val}")

print(f"\nTotal NNZ (Lower) entries: {NNZ}")

import numpy as np

def f(x):
    x0, x1, x2, x3 = x
    return np.array([
        x0**2 + x2**2,
        2*x0**2 + 3*x2**2,
        2*x1**2,
        x1**4 + x2**4 * x1**2,
        x3**3
    ])
def dense_hessian_f_weighted(x, lambd):
    x0, x1, x2, x3 = x
    H = np.zeros((4, 4))

    # f0 = x0^2 + x2^2
    H[0, 0] += 2 * lambd[0]
    H[2, 2] += 2 * lambd[0]

    # f1 = 2*x0^2 + 3*x2^2
    H[0, 0] += 4 * lambd[1]
    H[2, 2] += 6 * lambd[1]

    # f2 = 2*x1^2
    H[1, 1] += 4 * lambd[2]

    # f3 = x1^4 + x2^4 * x1^2
    # ∂²f3/∂x1² = 12*x1^2 + 2*x2^4
    H[1, 1] += (12 * x1**2 + 2 * x2**4) * lambd[3]

    # ∂²f3/∂x2² = 12*x2^2 * x1^2
    H[2, 2] += (12 * x2**2 * x1**2) * lambd[3]

    # ∂²f3/∂x1∂x2 = 8 * x2^3 * x1
    mix = 8 * x2**3 * x1 * lambd[3]
    H[1, 2] += mix
    H[2, 1] += mix  # symmetric

    # f4 = x3^3
    # ∂²f4/∂x3² = 6 * x3
    H[3, 3] += 6 * x3 * lambd[4]

    return H

def jvp(x, s):
    x0, x1, x2, x3 = x
    s0, s1, s2, s3 = s

    return np.array([
        2 * x0 * s0 + 2 * x2 * s2,
        4 * x0 * s0 + 6 * x2 * s2,
        4 * x1 * s1,
        4 * x1**3 * s1 + (4 * x2**3 * x1**2) * s2 + (2 * x2**4 * x1) * s1,
        3 * x3**2 * s3
    ])


def seed_vector(color):
    s = np.zeros(4)
    for i in color:
        s[i] = 1.0
    return s

def compute_hessian_entries(x, lambd, h=1e-4):
    NNZ = max(flat_idx for pairs in Q.values() for _, flat_idx in pairs) + 1
    Hflat = np.zeros(NNZ)

    # Precompute base JVP for each color seed
    J = {}
    for c, color in enumerate(COLORS):
        s = seed_vector(color)
        J[c] = jvp(x, s)

    # Loop over color pairs to do perturbed evaluations and finite differences
    for c1 in range(len(COLORS)):
        s1 = seed_vector(COLORS[c1])
        x1 = x + h * s1  # perturb x in direction of color c1

        for c2 in range(c1 + 1):
            s2 = seed_vector(COLORS[c2])
            J_pert = jvp(x1, s2)

            # For each nnz entry associated with color pair (c1,c2)
            for rows, nnz_index in Q[(c1, c2)]:
                out = 0
                for r in rows:
                    Hflat[nnz_index] += lambd[r] * (J_pert[r] -  J[c2][r]) / h

    return Hflat


def richardson_extrapolation(f, x, lambd, h, n=2, order=1):
    """
    Perform Richardson extrapolation for a function that returns arrays.
    
    Parameters:
    - f: function to extrapolate (returns np.array)
    - x: point at which to evaluate the function
    - h: initial step size
    - n: number of extrapolation steps
    
    Returns:
    - extrapolated array
    - extrapolation table (3D array: steps x steps x array_dimensions)
    """
    # First evaluation to determine array shape
    sample_output = f(x, lambd, h)
    array_shape = sample_output.shape
    
    # Initialize the extrapolation table as a 3D array
    table = np.zeros((n, n, *array_shape))
    
    # Fill first column with function evaluations at different h values
    for i in range(n):
        table[i, 0] = f(x, lambd, h / (2**i))
    
    # Perform Richardson extrapolation (array operations)
    for j in range(1, n):
        factor = 2**(order * j)
        for i in range(j, n):
            table[i, j] = (factor * table[i, j-1] - table[i-1, j-1]) / (factor - 1)
    
    return table[n-1, n-1], table


x = np.array([1, 2, 4, 5])
lambd = np.array([1, 2, 5, 1, -5])

H = compute_hessian_entries(x, lambd, 1e-8)
print(dense_hessian_f_weighted(x, lambd))
print(H)

A, B = richardson_extrapolation(compute_hessian_entries, x, lambd, 1e-4, n=2)

print(A)



\section{Numerical Hessian with Jacobian Column-Coloring and Partial Lambda Adjoint Product}
\begin{algorithm}[H]
	\caption{Numerical Hessian with Jacobian Column-Coloring, JVP, and Partial Lambda Product}
	\begin{algorithmic}[1]
		\REQUIRE 
		Vector-valued function $\v{f} : \mathbb{R}^n \to \mathbb{R}^m$, \\
		Adjoint vector $\v{\lambda} \in \mathbb{R}^m$, \\
		Jacobian-vector product operator $\mathrm{JVP}(\cdot)$, \\
		Column coloring $\mathcal{C} = \{c_1, \dots, c_k\}$, \\
		Step size $h$, \\
		Mapping $Q: \mathcal{C} \times \mathcal{C} \to (\mathcal{I}, \mathcal{N})$, where \(\mathcal{I}\) are indices for finite differences and \(\mathcal{N}\) are indices in flat Hessian, \\
		Flat Hessian buffer $H^{\mathrm{flat}} \in \mathbb{R}^{\mathrm{nnz}}$, initialized to zero
		\vspace{0.5em}
		
		\STATE Evaluate base JVPs: for all colors $c$, compute 
		\[
		J_{c} := \nabla \v{f}(\v{x}) \cdot \v{s}_c
		\]
		\STATE Precompute weighted base JVPs for each color:
		\[
		\hat{J}_c := \v{\lambda}^\top J_c = \sum_{r=1}^m \lambda_r J_c[r]
		\]
		
		\FOR{each color $c_1 \in \mathcal{C}$}
		\STATE Define seed vector $\v{s}_{c_1}$ with Algorithm 2
		\STATE Perturb point: 
		\[
		\v{x}_{c_1} = \v{x} + h \cdot \v{s}_{c_1}
		\]
		
		\FOR{each color $c_2 \in \mathcal{C}$ with $c_2 \leq c_1$} 
		\STATE Define seed vector $\v{s}_{c_2}$ with Algorithm 2
		\STATE Evaluate perturbed JVP: 
		\[
		J_{c_1, c_2} := \nabla \v{f}(\v{x}_{c_1}) \cdot \v{s}_{c_2}
		\]
		
		\FOR{all associated $(\mathrm{rows}, \mathrm{nz}) \in Q(c_1, c_2)$}
		\STATE Compute finite difference weighted by $\v{\lambda}$:
		\[
		H^{\mathrm{flat}}_{\mathrm{nz}} = \sum_{r \in \mathrm{rows}} \frac{\lambda_r \cdot J_{c_1,c_2}[r] - \hat{J}_{c_2}[r]}{h}
		\]
		
		\ENDFOR
		\ENDFOR
		\ENDFOR
		
		\RETURN $H^{\mathrm{flat}}$
	\end{algorithmic}
\end{algorithm}

\textbf{Explanation:}

This algorithm approximates the action of the adjoint-weighted Hessian $\v{\lambda}^\top \nabla^2 \v{f}(\v{x})$ using finite differences of Jacobian-vector products (JVPs), while exploiting sparsity through Jacobian column coloring. The method computes only the required components of the Hessian, avoiding its full construction.

\begin{itemize}
	\item Let $\v{f} : \mathbb{R}^n \to \mathbb{R}^m$ be a vector-valued function, and $\v{\lambda} \in \mathbb{R}^m$ the adjoint vector defining the scalar-valued function $\v{g}(\v{x}) := \v{\lambda}^\top \v{f}(\v{x})$.
	
	\item The gradient of $g$ is given by:
	\[
	\nabla \v{g}(\v{x}) = \nabla \v{f}(\v{x})^\top \v{\lambda} = J(\v{x})^\top \v{\lambda}
	\]
	and the Hessian is:
	\[
	\nabla^2 \v{g}(\v{x}) = \sum_{r=1}^m \lambda_r \cdot \nabla^2 f_r(\v{x})
	\]
	
	\item The column coloring partitions variables into independent groups $\mathcal{C} = \{c_1, \dots, c_k\}$ such that no two variables in the same color affect the same row of the Jacobian. This allows perturbing all variables in a color group $c$ simultaneously using a seed vector $\v{s}_c$.
	
	\item For each color pair $(c_1, c_2)$, the algorithm computes a finite-difference approximation to second-order directional derivatives in directions $\v{s}_{c_1}$ and $\v{s}_{c_2}$:
	\[
	\frac{1}{h} \left( \nabla \v{f}(\v{x} + h \v{s}_{c_1}) \cdot \v{s}_{c_2} - \nabla \v{f}(\v{x}) \cdot \v{s}_{c_2} \right)
	= \frac{J_{c_1,c_2} - J_{c_2}}{h}
	\]
	
	\item These directional derivatives are contracted with $\v{\lambda}$:
	\[
	\v{\lambda}^\top \left( \frac{J_{c_1,c_2} - J_{c_2}}{h} \right)
	= \sum_{r=1}^m \lambda_r \cdot \frac{J_{c_1,c_2}[r] - J_{c_2}[r]}{h}
	\]
	
	\item For each color pair $(c_1, c_2)$, the sparse structure mapping $Q(c_1, c_2)$ provides a list of entries to be updated: each element $(\mathrm{rows}, \mathrm{nz})$ indicates that the flat Hessian index $\mathrm{nz}$ corresponds to a variable pair $(i,j)$, and that this entry is influenced by function components indexed by $r \in \mathrm{rows}$.
	
	\item To reduce redundant multiplications, we precompute:
	\[
	\hat{J}_{c_2}[r] := \lambda_r \cdot J_{c_2}[r]
	\]
	
	\item The final Hessian entry is then assembled as:
	\[
	H^{\mathrm{flat}}_{\mathrm{nz}} \approx \sum_{r \in \mathrm{rows}} \lambda_r \cdot \frac{J_{c_1,c_2}[r] - J_{c_2}[r]}{h}
	= \sum_{r \in \mathrm{rows}} \frac{\lambda_r \cdot J_{c_1,c_2}[r] - \hat{J}_{c_2}[r]}{h}
	\]
	
	\item This approach evaluates all variable combinations defined by $(c_1, c_2)$ simultaneously, leveraging the independence structure imposed by the coloring. Since the variables in the same color do not interfere in any row of $\v{f}$, we can perturb them all at once, which drastically reduces the number of required JVP evaluations.
	
\end{itemize}

This method efficiently recovers the sparse structure and values of the adjoint-weighted Hessian by batching second derivative computations and exploiting sparsity at both the Jacobian and Hessian levels.

	

\end{document}