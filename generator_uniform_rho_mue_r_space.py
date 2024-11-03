# density functional minimizer/bulk rho and mue generator...

# this part of the code generates the bulk values of the density and the chemical potential in uniform way for each r space point supplied by the r-space files for all kind of the interactions specified in the interaction data json data profile... 

def bulk_rho_mue_generator():
    import numpy as np
    import json
    from scipy.integrate import quad

    
    
    # File paths
    json_file_species = "input_species_properties.json"
    json_file_interaction = "input_interactions_properties.json"
    json_file_thermodynamics = "input_thermodynamic_properties.json"
    r_space_file = "r_space.txt"
    output_file = "rho_mue_r_space.txt"



    # Load thermodynamic properties
    with open(json_file_thermodynamics, "r") as file:
        data_thermodynamic = json.load(file)
    total_rho = data_thermodynamic["thermodynamic_properties"]["rho"]



    # Load species properties
    with open(json_file_species, 'r') as file:
        data_species = json.load(file)
    species = {k: v["rho_frac"] * total_rho for k, v in data_species["species"].items()}  # Calculate rho for each species
    
    
    
    # Load interaction properties
    with open(json_file_interaction, 'r') as file:
        data_interactions = json.load(file)
    interactions = data_interactions["interactions"]
    
    

    # Define interaction potential based on data
    def interaction_potential(r, epsilon, sigma, interaction_type):
        """Calculate the potential based on interaction type."""
        if interaction_type == "lj":
            # Standard Lennard-Jones potential with cutoff at 5 * sigma
            if r < 2**(1/6) * sigma:
                return -epsilon
            elif 2**(1/6) * sigma <= r < 5 * sigma:
                return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
            else:
                return 0
                
                
        elif interaction_type == "wca":
            # Weeks-Chandler-Andersen potential, truncated and shifted Lennard-Jones
            if r < 2**(1/6) * sigma:
                return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6) + epsilon
            else:
                return 0
        
        
        elif interaction_type == "gs":
            # Gaussian interaction potential
            return epsilon * np.exp(-((r / sigma)**2))
        
        
        elif interaction_type == "yk":
            # Yukawa potential with a range controlled by sigma
            kappa = 1.0 / sigma  # Adjust kappa based on desired screening length
            return epsilon * np.exp(-kappa * r) / r if r != 0 else 0
        
        
        elif interaction_type == "hc":
            # Hard-core repulsive potential
            return float('inf') if r < sigma else 0
        
        
        else:
            # Default case if interaction type is not recognized
            return 0

    bulk_mue = {}

    
    for species_type, rho_value in species.items():
        # Initialize total chemical potential for the current species
        mue_total = 0
        for other_species, rho_other in species.items():
            # Determine the interaction data key
            interaction_key1 = f"{species_type}{other_species}"
            interaction_key2 = f"{other_species}{species_type}"

            if interaction_key1 in interactions["primary"]:
                interaction_data = interactions["primary"][interaction_key1]
            elif interaction_key2 in interactions:
                interaction_data = interactions["primary"][interaction_key2]
            else:
                # Skip if no interaction data is found
                continue

            # Unpack interaction parameters
            epsilon = interaction_data["epsilon"]
            sigma_ij = interaction_data["sigma"]
            interaction_type = interaction_data["type"]
            cutoff = interaction_data["cutoff"]

            # Check if it's a hard-core (hc) interaction
            if interaction_type == "hc":
                # Approximate hard-core chemical potential contribution
                try:
                    mue_total += rho_other * np.log(1 - np.pi * rho_value * sigma_ij**3 / 6)
                except ValueError as e:
                    print(f"Error calculating hard-core contribution for {species_type}-{other_species}: {e}")
                continue  # Skip integration for hard-core interactions

            # Set the lower limit of integration based on interaction type
            r_min = sigma_ij if interaction_type in ["lj", "wca", "yk"] else 0

            # Integrate the potential from r_min to cutoff
            try:
                integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
                mue_total += rho_other * integral_result
            except Exception as e:
                print(f"Error integrating potential for {species_type}-{other_species}: {e}")
                continue

        for other_species, rho_other in species.items():
            # Determine the interaction data key
            
            interaction_key1 = f"{species_type}{other_species}"
            interaction_key2 = f"{other_species}{species_type}"

            if interaction_key1 in interactions["secondary"]:
                interaction_data = interactions["secondary"][interaction_key1]
            elif interaction_key2 in interactions:
                interaction_data = interactions["secondary"][interaction_key2]
            else:
                # Skip if no interaction data is found
                continue

            # Unpack interaction parameters
            epsilon = interaction_data["epsilon"]
            sigma_ij = interaction_data["sigma"]
            interaction_type = interaction_data["type"]
            cutoff = interaction_data["cutoff"]

            # Check if it's a hard-core (hc) interaction
            if interaction_type == "hc":
                # Approximate hard-core chemical potential contribution
                try:
                    mue_total += rho_other * np.log(1 - np.pi * rho_value * sigma_ij**3 / 6)
                except ValueError as e:
                    print(f"Error calculating hard-core contribution for {species_type}-{other_species}: {e}")
                continue  # Skip integration for hard-core interactions

            # Set the lower limit of integration based on interaction type
            r_min = sigma_ij if interaction_type in ["lj", "wca", "yk"] else 0

            # Integrate the potential from r_min to cutoff
            try:
                integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
                mue_total += rho_other * integral_result
            except Exception as e:
                print(f"Error integrating potential for {species_type}-{other_species}: {e}")
                continue
                
        for other_species, rho_other in species.items():
            # Determine the interaction data key
            
            interaction_key1 = f"{species_type}{other_species}"
            interaction_key2 = f"{other_species}{species_type}"

            if interaction_key1 in interactions["tertiary"]:
                interaction_data = interactions["tertiary"][interaction_key1]
            elif interaction_key2 in interactions:
                interaction_data = interactions["tertiary"][interaction_key2]
            else:
                # Skip if no interaction data is found
                continue

            # Unpack interaction parameters
            epsilon = interaction_data["epsilon"]
            sigma_ij = interaction_data["sigma"]
            interaction_type = interaction_data["type"]
            cutoff = interaction_data["cutoff"]

            # Check if it's a hard-core (hc) interaction
            if interaction_type == "hc":
                # Approximate hard-core chemical potential contribution
                try:
                    mue_total += rho_other * np.log(1 - np.pi * rho_value * sigma_ij**3 / 6)
                except ValueError as e:
                    print(f"Error calculating hard-core contribution for {species_type}-{other_species}: {e}")
                continue  # Skip integration for hard-core interactions

            # Set the lower limit of integration based on interaction type
            r_min = sigma_ij if interaction_type in ["lj", "wca", "yk"] else 0

            # Integrate the potential from r_min to cutoff
            try:
                integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
                mue_total += rho_other * integral_result
            except Exception as e:
                print(f"Error integrating potential for {species_type}-{other_species}: {e}")
                continue


        
        # Store the total chemical potential for the species
        bulk_mue[species_type] = mue_total

    # Load r-space data
    r_space_data = np.loadtxt(r_space_file)

    # Create output data with r positions, rho, and chemical potential values for each point
    output_data = []
    for r_point in r_space_data:
        row = list(r_point)  # Initial columns are the position (x, y, z) in r-space
        for species_type, rho_value in species.items():
            row.append(rho_value)  # Append rho value for the species
           
            row.append(bulk_mue[species_type])  # Append the calculated chemical potential for the species
        output_data.append(row)

    # Save output to text file
    np.savetxt(output_file, output_data, fmt="%.6f", header="x y z rho mue")

    print(f"...... uniform rho and mue has been assigned to each r space points and exported to the appropriate points ......\n\n\n")

# Run the function
# bulk_rho_mue_generator()
