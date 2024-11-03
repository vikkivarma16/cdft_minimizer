import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.special import jv  # Bessel function of the first kind

# Function to generate a Gaussian potential
def gaussian_potential(r, epsilon, sigma):
    return epsilon * np.exp(-(r ** 2) / (2 * sigma ** 2))

# Bessel Fourier Transform (Hankel Transform)
def bessel_fourier_transform(r, V_r):
    k_space = np.linspace(0.1, 10, len(r))  # Choose an appropriate k-space range
    V_k = np.zeros_like(k_space)

    # Calculate the Hankel transform using numerical integration
    for i, k in enumerate(k_space):
        integrand = r * V_r * jv(0, k * r)  # Bessel function of the first kind
        V_k[i] = simps(integrand, r)  # Use Simpson's rule for integration

    return V_k, k_space

# Regular Fourier Transform
def regular_fourier_transform(r, V_r):
    V_k = np.fft.fft(V_r)
    k_space = np.fft.fftfreq(len(r), r[1] - r[0])
    return V_k, k_space

# Main function to generate Gaussian and perform transforms
def main():
    # Parameters for Gaussian
    epsilon = 1.0  # Height of the Gaussian
    sigma = 1.0    # Width of the Gaussian
    r = np.linspace(0, 10, 1000)  # Distance values

    # Generate Gaussian potential
    V_r = gaussian_potential(r, epsilon, sigma)

    # Perform Bessel Fourier Transform
    V_k_bessel, k_space_bessel = bessel_fourier_transform(r, V_r)

    # Perform regular Fourier Transform
    V_k_regular, k_space_regular = regular_fourier_transform(r, V_r)

    # Plotting
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Plot Gaussian potential
    axes[0].plot(r, V_r, label='Gaussian Potential V(r)', color='blue')
    axes[0].set_title('Gaussian Potential V(r)')
    axes[0].set_xlabel('r')
    axes[0].set_ylabel('V(r)')
    axes[0].grid(True)
    axes[0].legend()

    # Plot Bessel Fourier Transform
    axes[1].plot(k_space_bessel, np.real(V_k_bessel), label='Bessel FT V(k)', color='orange')
    axes[1].set_title('Bessel Fourier Transform V(k)')
    axes[1].set_xlabel('k')
    axes[1].set_ylabel('V(k)')
    axes[1].grid(True)
    axes[1].legend()

    # Plot Regular Fourier Transform
    axes[2].plot(k_space_regular, np.real(V_k_regular), label='Regular FT V(k)', color='green')
    axes[2].set_title('Regular Fourier Transform V(k)')
    axes[2].set_xlabel('k')
    axes[2].set_ylabel('V(k)')
    axes[2].grid(True)
    axes[2].legend()

    plt.tight_layout()
    plt.show()

# Execute the main function
main()

