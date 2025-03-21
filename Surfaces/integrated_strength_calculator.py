import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from scipy.special import j0  # Bessel function of the first kind (order 0)

# Mie Potential Function (Avoiding Hard-Core Repulsion)
def V_mie(r, epsilon):
    if r < sigma:
        return 0  # Ignore r < sigma (hard-core repulsion)
    return 4 * epsilon * ((alpha / r) ** 48.0 - (alpha / r) ** 24.0)

# Gaussian Potential Function
def V_gaussian(r, epsilon):
    return epsilon * np.exp(- (r / sigma) ** 2)

# Hankel Transform (Fourier Transform for Radial Functions)
def hankel_transform(V_func, k, r_min, r_max, epsilon):
    integral, _ = integrate.quad(lambda r: 4 * np.pi * r**2 * V_func(r, epsilon) * j0(k * r), r_min, r_max, limit=10000)
    return integral

# Range of k-values
k_values = np.linspace(0.00001, 5, 100)

# Compute Fourier Transform for Gaussian Potential
epsilon = 2.0  # Energy scale
alpha = 1     # Length scale for Mie potential
sigma = 0.665      # Hard-core size (minimum range)
cutoff = 100.0 # Integration limit

r_min = 0
r_max = cutoff

Vk_gaussian = [hankel_transform(V_gaussian, k, r_min, r_max, epsilon) for k in k_values]

# Compute Fourier Transform for Mie Potential
epsilon = 1.0  
alpha = 1.0   
sigma = 1.0  
cutoff = 5.0  

r_min = sigma
r_max = cutoff

Vk_mie = [hankel_transform(V_mie, k, r_min, r_max, epsilon) for k in k_values]

# Plot results and save as high-resolution PNG
plt.figure(figsize=(8, 5))
plt.plot(k_values, Vk_mie, label="Mie Potential", linestyle="--", color="red")
plt.plot(k_values, Vk_gaussian, label="Gaussian Potential", linestyle="-", color="blue")
plt.xlabel("Wavevector k")
plt.ylabel("Fourier Transform V(k)")
plt.title("Fourier Transform (Hankel Transform) of Potentials")
plt.legend()
plt.grid()
plt.savefig("fourier_transform.png", dpi=300)  # Save high-resolution image

# Save data to a text file
data = np.column_stack((k_values, Vk_mie, Vk_gaussian))
np.savetxt("fourier_transform_data.txt", data, header="k_values  Vk_mie  Vk_gaussian", fmt="%.6e")

# Print values at k → 0
print(f"Mie Potential at k=0: {Vk_mie[0]:.5f}")
print(f"Gaussian Potential at k=0: {Vk_gaussian[0]:.5f}")

