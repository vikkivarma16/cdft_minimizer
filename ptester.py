import numpy as np

def chemical_potential_cs(rho, sigma=1.0, Lambda=1.0):
    """
    Calculate the chemical potential of a hard sphere fluid using the Carnahan-Starling (CS) expression.
    
    Parameters:
    rho (float): Number density of the hard sphere fluid.
    sigma (float): Diameter of the hard spheres. Default is 1.0.
    Lambda (float): Thermal de Broglie wavelength. Default is 1.0, assuming reduced units.
    
    Returns:
    float: Chemical potential of the hard sphere fluid.
    """
    # Calculate the packing fraction, eta
    eta = (np.pi * rho * sigma**3) / 6.0
    
    # Check if eta is close to 1 to avoid division by zero
    if eta >= 1.0:
        raise ValueError("Packing fraction (eta) must be less than 1 for the Carnahan-Starling equation.")
    
    # Carnahan-Starling expression for the chemical potential (beta mu)
    beta_mu = np.log(rho * Lambda**3) + eta* (8 - 9 * eta + 3 * eta**2) / (1 - eta)**3
    e
    print (np.log(rho * Lambda**3))
    
    print (eta*(8 - 9 * eta + 3 * eta**3) / (1 - eta)**3)
    return beta_mu

# Example usage
rho_values = [0.1, 0.2, 0.3, 0.4, 0.5]  # Example densities
chemical_potentials = [chemical_potential_cs(rho) for rho in rho_values]

for rho, mu in zip(rho_values, chemical_potentials):
    print(f"Density: {rho}, Chemical Potential (beta*mu): {mu}")

