import numpy as np
import json

# Constants
EPSILON = 1e-6  # Small value to avoid division by zero
PI = np.pi

# Function to read particle interaction data from the JSON file
def read_particle_data_from_json(json_file):
    with open(json_file, 'r') as file:
        data = json.load(file)
    
    interaction_types = {}
    closest_distances = {}
    
    # Parse interaction data from JSON
    for pair_type, values in data["interactions"]["primary"].items():
        if values["type"] == "hc":
            interaction_types[pair_type] = values["type"]
            closest_distances[pair_type] = values["sigma"]
            
    for pair_type, values in data["interactions"]["secondary"].items():
        if values["type"] == "hc":
            interaction_types[pair_type] = values["type"]
            closest_distances[pair_type] = values["sigma"]
            
    for pair_type, values in data["interactions"]["tertiary"].items():
        if values["type"] == "hc":
            interaction_types[pair_type] = values["type"]
            closest_distances[pair_type] = values["sigma"]
        
    return interaction_types, closest_distances

# Function to identify unique particle types in the system
def identify_particle_types(interaction_types):
    particle_types = []
    
    for pair_type in interaction_types.keys():
        if pair_type[0] not in particle_types:
            particle_types.append(pair_type[0])
        if pair_type[1] not in particle_types:
            particle_types.append(pair_type[1])
    return particle_types

# Function to calculate particle sizes based on closest distances
def calculate_particle_sizes(closest_distances, particle_types):
    particle_sizes = {}
    
    for particle_type in particle_types:
        interaction_key = particle_type * 2  # e.g., 'a' -> 'aa', 'b' -> 'bb'
        # this section can handle the different interactions imposed on the same kind of the particles ...
        for pair_type, sigma in closest_distances.items():
            if interaction_key == pair_type:
                particle_sizes[particle_type] = closest_distances[interaction_key] / 2
            
    return particle_sizes

# Function to calculate weight function in k-space
def calculate_weight_function_k_space(particle_sizes, k_space, dimension):
    weight_functions = {}
    
    if dimension == 1:
        for particle_type, size in particle_sizes.items():
            weight_function = []
            
            for kx, ky, kz in k_space:
                weight_vector = [kx, ky, kz]
                
                
                k_value = np.sqrt(kx**2.0 + ky**2.0 + kz**2.0)
                
                if k_value < EPSILON:
                    weight_vector.extend([
                        1,                       # Weight function at k=0
                        size * 0.5,              # Additional weight terms
                        PI * size**2 * 0.25,
                        PI * 0.125 * size**3,
                        0,                  # n1_x
                        0                  # n2_x
                    ])
                else:
                    weight_vector.extend([
                        np.sin(k_value * PI * size) / (k_value * size * PI),
                        np.sin(k_value * PI * size) / (2.0 * k_value * PI),
                        size * np.sin(k_value * PI * size) / k_value,
                        (np.sin(k_value * PI * size) / (2.0 * k_value**3 * PI**2) - 
                         size * np.cos(k_value * PI * size) / (2.0 * k_value**2 * PI)),
                        (k_value * PI * size * np.cos(k_value * PI * size) - 
                         np.sin(k_value * PI * size)) / (2.0 * size * PI**2 * k_value**2),
                                              # n1_y, n1_z
                        (k_value * PI * size * np.cos(k_value * PI * size) - 
                         np.sin(k_value * PI * size)) / (k_value**2 * PI)                       # n2_y, n2_z
                    ])
                
                weight_function.append(weight_vector)
                
               
                
            weight_functions[particle_type] = np.array(weight_function)
    
    return weight_functions

# Function to export weight functions to files
def export_weight_functions_to_files(weight_functions):
    for particle_type, weight_function in weight_functions.items():
        file_name = f"weight_FMT_k_space_{particle_type}.txt"
        np.savetxt(file_name, weight_function)
        
        

# Main function to execute the calculations
def fmt_weights_1d():
    # File paths
    json_file_path = 'input_interactions_properties.json'
    k_space_file_path = 'k_space.txt'
    
    # Read interaction data
    interaction_types, closest_distances = read_particle_data_from_json(json_file_path)
    
    # Identify particle types and calculate particle sizes
    particle_types = identify_particle_types(interaction_types)
    particle_sizes = calculate_particle_sizes(closest_distances, particle_types)
    
    # Load k-space data from file
    k_space = np.loadtxt(k_space_file_path)
    dimension = 1  # Set dimension for calculation
    
    # Calculate weight functions
    weight_functions = calculate_weight_function_k_space(particle_sizes, k_space, dimension)
    
    # Export weight functions to files
    export_weight_functions_to_files(weight_functions)

# Run the main function

# fmt_weights_1d()
