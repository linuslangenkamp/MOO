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
from sympy import symbols, exp, diff, lambdify, Matrix, sin, cos

# Define symbols
x1, x2, x3, x4, u = symbols('x1 x2 x3 x4 u')

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
F = []
dec = 0
if dec == 0:
    F.append((-k1 * x1) - (k3 + k4 + k5) * x1 * x2)
    F.append(k1 * x1 - k2 * x2 + k3 * x1 * x2)
    F.append(k2 * x2 + k4 * x1 * x2)
    F.append(k5 * x1 * x2)
    F.append((-k1 * x1) - (k3 + k4 + k5) * x1 * x2)
    V = [x1, x2, x3, x4, u]
    x1_val = 1.35266244529518659e-01
    x2_val = 3.55066498452544121e-01
    x3_val = 3.34360423673394747e-01
    x4_val = 1.75306833344556989e-01
    u_val = 9.70337348817454814e-01
elif dec == 1:
    F.append(sin(x1 + x2 - 1))
    F.append(cos(x2 + x3 - 2))
    F.append((x3 + u) / (x3**2 + u**2))
    F.append(-1 + x3*x3)
    V = [x1, x2, x3, u]
    x1_val = 1.478975818445809
    x2_val = 1.470723942456369
    x3_val = 1.11596063078672
    u_val = 4.999999970445199

lambd = [0, 0, 0, 0, 1]

#f = sum(lambd[i] * F[i] for i in range(4))
f = F[4]
print(f)
# Variables vector
vars_vec = Matrix(V)

# Compute Hessian
hessian = f.diff(vars_vec).jacobian(vars_vec)

# Lambdify Hessian for numerical evaluation
hessian_func = lambdify((x1, x2, x3, x4, u), hessian, 'numpy')

# Evaluate at given point



np.set_printoptions(precision=5)

hessian_numeric = np.array(hessian_func(x1_val, x2_val, x3_val, x4_val, u_val), dtype=np.float128)
print()
print(hessian_numeric)



#### test simulation

t = [0, 0.08, 0.16, 0.24, 0.32, 0.4, 0.48, 0.56, 0.64, 0.72, 0.8, 0.88, 0.96, 1.04, 1.12, 1.2, 1.28, 1.36, 1.44, 1.52, 1.6, 1.68, 1.76, 1.84, 1.92, 2, 2.08, 2.16, 2.24, 2.32, 2.4, 2.48, 2.56, 2.64, 2.72, 2.8, 2.88, 2.96, 3.04, 3.12, 3.2, 3.28, 3.36, 3.44, 3.52, 3.6, 3.68, 3.76, 3.84, 3.92, 4, 4.08, 4.16, 4.24, 4.32, 4.4, 4.48, 4.56, 4.64, 4.72, 4.8, 4.88, 4.96, 5.04, 5.12, 5.2, 5.28, 5.36, 5.44, 5.52, 5.6, 5.68, 5.76, 5.84, 5.92, 6, 6.08, 6.16, 6.24, 6.32, 6.4, 6.48, 6.56, 6.64, 6.72, 6.8, 6.88, 6.96, 7.04, 7.12, 7.2, 7.28, 7.36, 7.44, 7.52, 7.6, 7.68, 7.76, 7.84, 7.92, 8]
x = [
  [1, 0, 0, 0],
  [0.999731, 0.000263691, 3.03727e-06, 2.1471e-06],
  [0.999442, 0.000537228, 1.22981e-05, 8.69334e-06],
  [0.999131, 0.000820966, 2.80119e-05, 1.98002e-05],
  [0.998799, 0.00111528, 5.04158e-05, 3.56344e-05],
  [0.998443, 0.00142054, 7.97553e-05, 5.63688e-05],
  [0.998064, 0.00173716, 0.000116284, 8.21818e-05],
  [0.997661, 0.00206553, 0.000160266, 0.000113258],
  [0.997232, 0.00240608, 0.000211971, 0.000149787],
  [0.996777, 0.00275924, 0.00027168, 0.000191967],
  [0.996295, 0.00312546, 0.000339684, 0.000240001],
  [0.995784, 0.00350521, 0.000416282, 0.000294098],
  [0.995245, 0.00389897, 0.000501783, 0.000354476],
  [0.994675, 0.00430722, 0.000596509, 0.000421358],
  [0.994074, 0.00473048, 0.000700788, 0.000494975],
  [0.99344, 0.00516926, 0.000814962, 0.000575565],
  [0.992773, 0.00562412, 0.000939383, 0.000663373],
  [0.992071, 0.0060956, 0.00107441, 0.000758653],
  [0.991334, 0.00658428, 0.00122043, 0.000861665],
  [0.990559, 0.00709076, 0.00137782, 0.000972679],
  [0.989745, 0.00761562, 0.00154697, 0.00109197],
  [0.988892, 0.00815951, 0.00172831, 0.00121982],
  [0.987998, 0.00872307, 0.00192225, 0.00135653],
  [0.987061, 0.00930694, 0.00212922, 0.0015024],
  [0.986081, 0.00991182, 0.00234969, 0.00165773],
  [0.985055, 0.0105384, 0.0025841, 0.00182285],
  [0.983982, 0.0111874, 0.00283294, 0.00199809],
  [0.98286, 0.0118595, 0.00309669, 0.00218377],
  [0.981688, 0.0125555, 0.00337585, 0.00238024],
  [0.980465, 0.0132762, 0.00367094, 0.00258787],
  [0.979188, 0.0140224, 0.00398249, 0.002807],
  [0.977856, 0.0147948, 0.00431105, 0.00303802],
  [0.976467, 0.0155942, 0.00465717, 0.0032813],
  [0.97502, 0.0164216, 0.00502142, 0.00353724],
  [0.973512, 0.0172778, 0.00540439, 0.00380622],
  [0.971941, 0.0181636, 0.00580669, 0.00408867],
  [0.970306, 0.01908, 0.00622892, 0.00438499],
  [0.968605, 0.0200279, 0.00667172, 0.00469561],
  [0.966835, 0.0210082, 0.00713574, 0.00502095],
  [0.964995, 0.0220217, 0.00762162, 0.00536148],
  [0.963083, 0.0230696, 0.00813004, 0.00571762],
  [0.961096, 0.0241527, 0.00866169, 0.00608984],
  [0.959032, 0.0252721, 0.00921726, 0.0064786],
  [0.956889, 0.0264287, 0.00979746, 0.00688438],
  [0.954666, 0.0276236, 0.010403, 0.00730764],
  [0.952359, 0.0288577, 0.0110347, 0.00774888],
  [0.949966, 0.0301321, 0.0116932, 0.00820858],
  [0.947486, 0.0314478, 0.0123793, 0.00868723],
  [0.944915, 0.0328059, 0.0130938, 0.00918534],
  [0.942252, 0.0342074, 0.0138375, 0.00970341],
  [0.939494, 0.0356533, 0.0146111, 0.0102419],
  [0.936639, 0.0371446, 0.0154154, 0.0108014],
  [0.933684, 0.0386825, 0.0162513, 0.0113824],
  [0.930627, 0.0402678, 0.0171196, 0.0119854],
  [0.927466, 0.0419018, 0.018021, 0.0126108],
  [0.924199, 0.0435853, 0.0189564, 0.0132593],
  [0.920823, 0.0453194, 0.0199267, 0.0139313],
  [0.917335, 0.047105, 0.0209326, 0.0146272],
  [0.913734, 0.0489431, 0.021975, 0.0153477],
  [0.910018, 0.0508346, 0.0230546, 0.0160932],
  [0.906183, 0.0527805, 0.0241724, 0.0168641],
  [0.902228, 0.0547817, 0.0253291, 0.0176609],
  [0.898151, 0.0568389, 0.0265255, 0.0184842],
  [0.89395, 0.0589529, 0.0277625, 0.0193343],
  [0.889623, 0.0611246, 0.0290408, 0.0202116],
  [0.885167, 0.0633547, 0.0303612, 0.0211167],
  [0.880582, 0.0656438, 0.0317244, 0.0220498],
  [0.875865, 0.0679924, 0.0331312, 0.0230114],
  [0.871015, 0.0704013, 0.0345823, 0.0240018],
  [0.866029, 0.0728708, 0.0360785, 0.0250213],
  [0.860908, 0.0754013, 0.0376203, 0.0260703],
  [0.855649, 0.0779932, 0.0392084, 0.027149],
  [0.850252, 0.0806467, 0.0408435, 0.0282577],
  [0.844715, 0.083362, 0.0425261, 0.0293965],
  [0.839038, 0.0861392, 0.0442567, 0.0305657],
  [0.83322, 0.0889782, 0.0460359, 0.0317654],
  [0.827261, 0.0918789, 0.0478641, 0.0329957],
  [0.82116, 0.0948411, 0.0497418, 0.0342566],
  [0.814918, 0.0978645, 0.0516693, 0.0355483],
  [0.808534, 0.100949, 0.053647, 0.0368705],
  [0.802009, 0.104093, 0.055675, 0.0382234],
  [0.795343, 0.107296, 0.0577537, 0.0396067],
  [0.788538, 0.110559, 0.0598831, 0.0410203],
  [0.781594, 0.113878, 0.0620635, 0.042464],
  [0.774513, 0.117255, 0.0642947, 0.0439375],
  [0.767296, 0.120687, 0.0665769, 0.0454404],
  [0.759945, 0.124173, 0.0689098, 0.0469725],
  [0.752462, 0.127711, 0.0712933, 0.0485332],
  [0.74485, 0.131301, 0.0737273, 0.050122],
  [0.737111, 0.134939, 0.0762114, 0.0517384],
  [0.729248, 0.138625, 0.0787451, 0.0533818],
  [0.721264, 0.142356, 0.0813282, 0.0550515],
  [0.713163, 0.14613, 0.08396, 0.0567469],
  [0.704948, 0.149945, 0.0866399, 0.058467],
  [0.696623, 0.153798, 0.0893674, 0.0602111],
  [0.688193, 0.157687, 0.0921416, 0.0619783],
  [0.679661, 0.161609, 0.0949618, 0.0637676],
  [0.671033, 0.165562, 0.097827, 0.065578],
  [0.662314, 0.169541, 0.100736, 0.0674085],
  [0.653508, 0.173545, 0.103689, 0.069258],
  [0.644622, 0.17757, 0.106683, 0.0711253],
]
u = [
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
  [700],
]


import matplotlib.pyplot as plt
import numpy as np

t = np.array(t)
x = np.array(x)

for i in range(x.shape[1]):
    plt.plot(t, x[:, i], label=f'x[{i}]')

plt.xlabel("Time")
plt.ylabel("State values")
plt.title("States over time")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()