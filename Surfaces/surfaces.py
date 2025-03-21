import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy import symbols, log, diff, lambdify
from tqdm import tqdm
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import Delaunay



# Define symbols
avc, bvc, cvc, vc = sp.symbols(' avc bvc cvc vc')
#ac, cc, vp1, vp2, cp1, cp2, ap1, ap2 = sp.symbols('ac cc vp1 vp2 cp1 cp2 ap1 ap2')
rhot, y, x = sp.symbols('rhot y x')



Pi = np.pi
rc = 1
ac = 4*Pi*rc**2
vc = 4/3*Pi*rc**3
cc = rc
rp1 = 1
cp1 = rp1
ap1 = 4*Pi*rp1**2
vp1 = 4/3*Pi*rp1**3
rp2 = 0.666*rp1
ap2 = 4*Pi*rp2**2
vp2 = 4/3*Pi*rp2**3
cp2 = rp2
sc  = 4*Pi*rc**2
sp1 = 4*Pi*rp1**2
sp2 = 4*Pi*rp2**2



avc1 = ac * cp1 + ap1 * cc + vp1
bvc1 = (ac * cp1) ** 2 / 2 + cc * ac * vp1
cvc1 = ac ** 2 * cc ** 2 * vp1 / 3
avc2 = ac * cp2 + ap2 * cc + vp2
bvc2 = (ac * cp2) ** 2 / 2 + cc * ac * vp2
cvc2 = ac ** 2 * cc ** 2 * vp2 / 3



n0c   = rhot * y
n0p1  = rhot * (1 - y) * (1 - x)
n0p2  = rhot * (1 - y) * x
n1c   = rhot * y * cc
n1p1  = rhot * (1 - y) * (1 - x) * cp1
n1p2  = rhot * (1 - y) * x * cp2
n2c   = rhot * y * (4 * np.pi * cc ** 2)
n2p1  = rhot * (1 - y) * (1 - x) * (4 * np.pi * cp1 ** 2)
n2p2  = rhot * (1 - y) * x * (4 * np.pi * cp2 ** 2)
etac  = rhot * y * vc
etap1 = rhot * (1 - y) * (1 - x) * vp1
etap2 = rhot * (1 - y) * x * vp2



alpha = (1 - etac) * sp.exp(-avc * etac / (vc * (1 - etac)) - bvc * (etac / (1 - etac)) ** 2 / vc ** 2 - cvc * (etac / (1 - etac)) ** 3 / vc ** 3)
alpha1 = alpha.subs({avc: avc1, bvc: bvc1, cvc: cvc1})
alpha2 = alpha.subs({avc: avc2, bvc: bvc2, cvc: cvc2})



fexcess = -n0c * sp.log(1 - rhot * y * vc) + n1c * n2c / (1 - rhot * y * vc) + n2c ** 3 / (24 * np.pi * (1 - rhot * y * vc) ** 2)
fext = fexcess/rhot -  (1 - y) * (1 - x) * sp.log(alpha1) -  (1 - y) * x * sp.log(alpha2)


rjj = 1
rii = 0.666*rjj
epsii = 2.0
epsjj = 2.0
rij   = (0.5 * (rii**2 + rjj**2))**0.5
epsij = 0.944 * epsii
vii   = rii**3 * epsii * np.pi**1.5
print("Vii", vii)

vjj   = rjj**3 * epsjj * np.pi**1.5
print("Vjj", vjj)

vij   = rij**3 * epsij * np.pi**1.5
print("Vij", vij)

vci   =  vii
vcc   = -vii
vcj   = vij
hatVpp = x**2 * vii + (1 - x)**2 * vjj + 2 * x * (1 - x) * vij
hatVcp = x * vci + (1 - x) * vcj
hatVcc = vcc



fid = sp.log(rhot) - 1
fmix = y * sp.log(y) + (1 - y) * x * sp.log((1 - y) * x) + (1 - y) * (1 - x) * sp.log((1 - y) * (1 - x))
fint   = 0.5 * rhot * hatVpp * (1 - y)**2
fintcp = rhot * y * (1 - y) * hatVcp
fintcc = 0.5 * rhot * y**2 * hatVcc




# Total free energy
ftotal = fid + fmix + fint + fintcp + fintcc + fext



mueA = ftotal + rhot * sp.diff(ftotal, rhot) + (1 - x) * sp.diff(ftotal, x) - y * sp.diff(ftotal, y) - y * (1 - x) * sp.diff(ftotal, x, y)
mueB = ftotal + rhot * sp.diff(ftotal, rhot) - sp.diff(ftotal, x) * x - sp.diff(ftotal, y) * y  + sp.diff(ftotal, x, y) * x * y 
mueC = ftotal + rhot * sp.diff(ftotal, rhot) + (1 - y) * sp.diff(ftotal, y)

# Pressure equation
pressure = sp.simplify(rhot ** 2 * sp.diff(ftotal, rhot))

mue_a_func = lambdify((rhot, x, y), mueA, 'numpy')
mue_b_func = lambdify((rhot, x, y), mueB, 'numpy')
mue_c_func = lambdify((rhot, x, y), mueC, 'numpy')




# Solve nonlinear equations

print("now the evaluation starts")

#Parameters for the mesh
rho_t_vals = np.linspace(2, 15, 100)
x_vals = np.linspace(0.01, 1, 100)
y_vals = np.linspace(0.00000000001, 0.01, 100)



specified_mue_c = 430  # Example target value for mu_c

# Create a mesh grid
rho_t_mesh, x_mesh, y_mesh = np.meshgrid(rho_t_vals, x_vals, y_vals, indexing='ij')

# Flatten the mesh grid for vectorized computation
rho_t_flat = rho_t_mesh.ravel()
x_flat = x_mesh.ravel()
y_flat = y_mesh.ravel()

# Compute mu_c values for all points in the mesh

data_points = []
mu_c_values = []
for i in tqdm(range(len(rho_t_flat)), desc="Processing grid points"):
    rho_t_val, x_val, y_val = rho_t_flat[i], x_flat[i], y_flat[i]
    try:
        mu_c_val = mue_c_func(rho_t_val, x_val, y_val)
        mu_c_values.append(mu_c_val)
        if np.isclose(mu_c_val, specified_mue_c, atol=0.1):
            data_points.append((rho_t_val, x_val, y_val))
    except Exception:
        mu_c_values.append(np.nan)
        continue

# Convert the result to a NumPy array
result = np.array(data_points)
mu_c_array = np.array(mu_c_values).reshape(rho_t_mesh.shape)


# Output
print(f"Found {len(result)} points satisfying mu_c = {specified_mue_c}")
if len(result) > 0:
    print("Example points:", result[:5])


if len(result) > 0:
    np.savetxt("result_coordinates_species_C.txt", result, header="", comments="", fmt="%.6f")
    print("Result coordinates exported to 'result_coordinates.txt'")




specified_mue_a = 44.46
data_points = []
mu_a_values = []
for i in tqdm(range(len(rho_t_flat)), desc="Processing grid points"):
    rho_t_val, x_val, y_val = rho_t_flat[i], x_flat[i], y_flat[i]
    try:
        mu_a_val = mue_a_func(rho_t_val, x_val, y_val)
        mu_a_values.append(mu_a_val)
        if np.isclose(mu_a_val, specified_mue_a, atol=0.1):
            data_points.append((rho_t_val, x_val, y_val))
    except Exception:
        mu_a_values.append(np.nan)
        continue

# Convert the result to a NumPy array
result = np.array(data_points)
mu_a_array = np.array(mu_a_values).reshape(rho_t_mesh.shape)


# Output
print(f"Found {len(result)} points satisfying mu_A = {specified_mue_a}")
if len(result) > 0:
    print("Example points:", result[:5])


if len(result) > 0:
    np.savetxt("result_coordinates_species_A.txt", result, header="", comments="", fmt="%.6f")
    print("Result coordinates exported to 'result_coordinates.txt'")






specified_mue_b = 80
data_points = []
mu_b_values = []
for i in tqdm(range(len(rho_t_flat)), desc="Processing grid points"):
    rho_t_val, x_val, y_val = rho_t_flat[i], x_flat[i], y_flat[i]
    try:
        mu_b_val = mue_b_func(rho_t_val, x_val, y_val)
        mu_b_values.append(mu_b_val)
        if np.isclose(mu_b_val, specified_mue_b, atol=0.1):
            data_points.append((rho_t_val, x_val, y_val))
    except Exception:
        mu_b_values.append(np.nan)
        continue

# Convert the result to a NumPy array
result = np.array(data_points)
mu_b_array = np.array(mu_b_values).reshape(rho_t_mesh.shape)


# Output
print(f"Found {len(result)} points satisfying mu_b = {specified_mue_b}")
if len(result) > 0:
    print("Example points:", result[:5])


if len(result) > 0:
    np.savetxt("result_coordinates_species_B.txt", result, header="", comments="", fmt="%.6f")
    print("Result coordinates exported to 'result_coordinates.txt'")



