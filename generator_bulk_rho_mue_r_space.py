# density functional minimizer/bulk rho and mue generator...

# this part of the code generates the bulk values of the density and the chemical potential in uniform way for each r space point supplied by the r-space files for all kind of the interactions specified in the interaction data json data profile... 

def bulk_rho_mue_r_space():
    import numpy as np
    import json
    from scipy.integrate import quad
    import calculator_pair_potential_custom
    
    
    # File paths
    
    json_file_particles_interactions = "input_data_particles_interactions_parameters.json"
    json_file_simulation_thermodynamics = "input_data_simulation_thermodynamic_parameters.json"
    r_space_file = "supplied_data_r_space.txt"
    output_file = "supplied_data_bulk_mue_rho_r_space.txt"



    # Load thermodynamic properties
    with open(json_file_simulation_thermodynamics, "r") as file:
        data_thermodynamic = json.load(file)
    total_rho = data_thermodynamic["simulation_thermodynamic_parameters"]["rho"]
    temperature =  data_thermodynamic["simulation_thermodynamic_parameters"]["temperature"]


    # Load interaction properties
    with open(json_file_particles_interactions, 'r') as file:
        data_interactions = json.load(file)
    interactions = data_interactions["particles_interactions_parameters"]["interactions"]
    data_species = data_interactions["particles_interactions_parameters"]
    species = {k: v["rho_frac"] * total_rho for k, v in data_species["species"].items()}  # Calculate rho for each species



    # Define interaction potential based on data
    def interaction_potential(r, epsilon, sigma, interaction_type):
        """Calculate the potential based on interaction type."""
        if interaction_type == "wca":
            # Standard Lennard-Jones potential with cutoff at 5 * sigma
            if r < 2**(1/6) * sigma:
                return -epsilon
            elif 2**(1/6) * sigma <= r < 5 * sigma:
                return 4 * epsilon * ((sigma / r)**12 - (sigma / r)**6)
            else:
                return 0
                
                
        elif interaction_type == "wca_1":
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
            return 0
            
        elif "custom" in interaction_type:
            dummy_string="pair_potential_integrant_"+interaction_type
            
            func = getattr(calculator_pair_potential_custom, dummy_string)
            v_r = np.zeros_like(r)
            func(r, v_r, epsilon, sigma)
            return v_r
            
            
            '''try:
                func = getattr(calculator_pair_potential_custom, dummy_string)
                v_r = np.zeros_like(r)
                func(r, v_r, epsilon, sigma)
                return v_r
            except AttributeError:
                print(f"       Error:    Unknown interaction type: {interaction_type}")
                exit(0)
            
            '''
            return func(r, epsilon, sigma)
        
        
        else:
            # Default case if interaction type is not recognized
            return 0

    bulk_mue = {}

    # hard core chemical potentials will be calculated correct way only for the additive mixture for the non - additive mixture you should look for some different theories ... and implement in the system.....
    
    
    y1=0
    y2=0
    y3=0
    y0=1 
    eta=0
    total_rho=0
    hc_sigma={}
    
    for species_type, rho_value in species.items():
        # Initialize total chemical potential for the current species
        mue_total = 0
        bulk_mue[species_type] = 0
        # for the hard core kind of interactions
        
        
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
            if interaction_type == "hc" and species_type == other_species:
                # Approximate hard-core chemical potential contribution
                try:
                    y1 += rho_other*sigma_ij
                    y2 += rho_other*sigma_ij*sigma_ij
                    y3 += rho_other*sigma_ij*sigma_ij*sigma_ij
                    total_rho += rho_other
                    eta += (1/6.0)*rho_other*np.pi*sigma_ij**3.0 
                    hc_sigma[species_type] = sigma_ij
                    
                except ValueError as e:
                    print(f"Error calculating hard-core contribution for {species_type}-{other_species}: {e}")
                continue  # Skip integration for hard-core interactions

            # Set the lower limit of integration based on interaction type
            r_min = sigma_ij if interaction_type in ["lj", "yk"] else 0
            
            integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
            mue_total += rho_other * integral_result
            '''
            # Integrate the potential from r_min to cutoff
            try:
                integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
                mue_total += rho_other * integral_result
            except Exception as e:
                print(f"Error integrating potential for {species_type}-{other_species}: {e}")
                continue
            '''

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
                    y1 += rho_other*sigma_ij
                    y2 += rho_other*sigma_ij*sigma_ij
                    y3 += rho_other*sigma_ij*sigma_ij*sigma_ij
                    total_rho += rho_other
                    eta += (1/6.0)*rho_other*np.pi*sigma_ij**3.0 
                    hc_sigma[species_type] = sigma_ij
                except ValueError as e:
                    print(f"Error calculating hard-core contribution for {species_type}-{other_species}: {e}")
                continue  # Skip integration for hard-core interactions

            # Set the lower limit of integration based on interaction type
            r_min = sigma_ij if interaction_type in ["lj", "yk"] else 0

            # Integrate the potential from r_min to cutoff
            integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
            mue_total += rho_other * integral_result
            
            
            '''
            try:
                integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
                mue_total += rho_other * integral_result
            except Exception as e:
                print(f"Error integrating potential for {species_type}-{other_species}: {e}")
                continue
            '''  
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
                    y1 += rho_other*sigma_ij
                    y2 += rho_other*sigma_ij*sigma_ij
                    y3 += rho_other*sigma_ij*sigma_ij*sigma_ij
                    total_rho += rho_other
                    eta += (1/6.0)*rho_other*np.pi*sigma_ij**3.0 
                    hc_sigma[species_type] = sigma_ij
                except ValueError as e:
                    print(f"Error calculating hard-core contribution for {species_type}-{other_species}: {e}")
                continue  # Skip integration for hard-core interactions

            # Set the lower limit of integration based on interaction type
            r_min = sigma_ij if interaction_type in ["lj", "yk"] else 0

            # Integrate the potential from r_min to cutoff
            integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
            mue_total += rho_other * integral_result
            '''
            try:
                integral_result, _ = quad(interaction_potential, r_min, cutoff, args=(epsilon, sigma_ij, interaction_type))
                mue_total += rho_other * integral_result
            except Exception as e:
                print(f"Error integrating potential for {species_type}-{other_species}: {e}")
                continue
            '''
        
        #print("mue before ideal", mue_total/temperature)
        # Store the total chemical potential for the species
        bulk_mue[species_type] = mue_total/temperature + np.log(species[species_type])
        #print ("mue ideal", np.log(species[species_type]))
    try:
        y1=y1/total_rho
        y2=y2/total_rho
        y3=y3/total_rho
        print ("... there are additive hard cores specified in the system and the bulk Carnahan-Sterling approximation have been done ...\n\n")
    except Exception as e:
        print(f"\n")
            
    for key, value in hc_sigma.items():
        sigma = float(hc_sigma[key])
        #print(sigma, "and ", bulk_mue[key])
        #hard_core =    ((sigma / y1) * (3 * eta / (1 - eta) + 3 * eta**2 / (1 - eta)**2)) + (sigma**2 / y2) * (9 * eta**2 / (2 * (1 - eta)**2)) - eta / (1 - eta)
        
        #hard_core = (6 * eta / (np.pi * y3)) * ((y0 * y2) / y3 + (y2**3) / y3**2 + (y0**2) / (3 * y3) )
        
        hard_core = eta* (8 - 9 * eta + 3 * eta**2) / (1 - eta)**3
        
        bulk_mue[key] += hard_core
        #print("hard_core values are given by:", hard_core)
        #print ("total_values are given as :", key, bulk_mue[key])
        
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

    print(f"\n\n... uniform rho and mue has been assigned to each r space points and exported to the supplied data file section ...\n\n\n")

# Run the function
#bulk_rho_mue_r_space()
