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


def richardson_extrapolation(f, x, lambd, h, n=2, order=1, base=10):
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
        table[i, 0] = f(x, lambd, h / (base**i))
    
    # Perform Richardson extrapolation (array operations)
    for j in range(1, n):
        factor = base**(order * j)
        for i in range(j, n):
            table[i, j] = (factor * table[i, j-1] - table[i-1, j-1]) / (factor - 1)
    
    return table[n-1, n-1], table


x = np.array([1, 2, 4, 5])
lambd = np.array([1, 2, 5, 1, -5])

"""
A, B = richardson_extrapolation(compute_hessian_entries, x, lambd, 1e-5, n=2, base=2)
print(A)

A, B = richardson_extrapolation(compute_hessian_entries, x, lambd, 1e-3, n=3, base=2)
print(A)
"""
H = compute_hessian_entries(x, lambd, 1e-5)
print(dense_hessian_f_weighted(x, lambd))
print(H)

# grid serach
import numpy as np

h_values = [1, 0.5, 1e-1, 5e-2, 1e-2, 5e-3, 1e-3, 5e-4, 1e-4, 5e-5, 1e-5, 5e-6, 1e-6, 5e-7, 1e-7, 5e-8, 1e-8, 5e-9, 1e-9]
n_values = [1, 2]
base_values = [1.25, 1.5, 1.75, 2, 3, 5, 10, 20, 50, 100]
exact = np.array([10, 580, 0, 1024, 782, -150])
results = []
for h in h_values:
    for n in n_values:
        for base in base_values:
            A, _ = richardson_extrapolation(compute_hessian_entries, x, lambd, h, n=n, order=1, base=base)
            error = np.linalg.norm(A - exact, np.inf)
            results.append({
                'h': h,
                'n': n,
                'base': base,
                'error': error,
                'result': A
            })


results.sort(key=lambda x: x['error'])

print("Top configurations by error:")
for i, res in enumerate(results):
    print(f"{i+1}. h={res['h']:.0e}, n={res['n']}, base={res['base']}, error={res['error']:.10f}")
    print("   Result: [" + ", ".join(f"{x:.10f}" for x in res['result']) + "]")
    print()

## 1-norm                                      | 2-norm                                      | oo-norm
#  n=4: h=1e+00, base=2,    error=0.0000000000 | n=4: h=1e+00, base=2,    error=0.0000000000 | n=4: h=1e+00, base=2,  error=0.0000000000
#  n=3: h=1e-03, base=1.75, error=0.0000000028 | n=3: h=1e-03, base=1.75, error=0.0000000018 | n=3: h=1e-02, base=10, error=0.0000000012
#  n=2: h=5e-05, base=1.75, error=0.0000000998 | n=2: h=5e-05, base=1.75, error=0.0000000775 | n=2: h=5e-05, base=10, error=0.0000000731
#  n=1: h=5e-08, base=x,    error=0.0000198226 | n=1: h=5e-08, base=x,    error=0.0000131057 | n=1: h=1e-08, base=x,  error=0.0000106304

# Re-execute due to state reset

import numpy as np
from sympy import symbols, exp, diff, lambdify, Matrix

# Define symbols
x1, x2, u = symbols('x1 x2 u')

# Constants
R = 1.9872
T = 700 * u

# Reaction rate expressions
k1 = exp(8.86 - 20300 / R / T)
k2 = exp(24.25 - 37400 / R / T)
k3 = exp(23.67 - 33800 / R / T)
k4 = exp(18.75 - 28200 / R / T)
k5 = exp(20.70 - 31000 / R / T)

# Define the function
f = (-k1 * x1) - (k3 + k4 + k5) * x1 * x2

# Variables vector
vars_vec = Matrix([x1, x2, u])

# Compute Hessian
hessian = f.diff(vars_vec).jacobian(vars_vec)

# Lambdify Hessian for numerical evaluation
hessian_func = lambdify((x1, x2, u), hessian, 'numpy')

# Evaluate at given point
x1_val = 0.9528284865048966
x2_val = 0.02791790905016496
u_val = 1.068785721450124
np.set_printoptions(precision=16)

hessian_numeric = np.array(hessian_func(x1_val, x2_val, u_val), dtype=np.float64)

print(hessian_numeric)
