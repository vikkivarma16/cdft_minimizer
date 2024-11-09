# density functional minimizer/executor

# this is the main code which read and write the data in an executable format and then also run the program for the calculation...






# this part read the input file and print the json file for the further processing of the input file in the selected format...







# Density Functional Minimizer/Executor
# Main code that reads, writes, and processes input data to execute the program

# Density Functional Minimizer/Executor
import numpy as np
import json
# Assuming you have pynufft installed
from pynufft import NUFFT

import pyfftw.interfaces as fftw





from generator_input_data_particles_interactions_parameters import data_exporter_particles_interactions_parameters as interactions_parameters
from generator_input_data_space_confinement_parameters import data_exporter_space_confinement_parameters as confinement_parameters
from generator_input_data_simulation_thermodynamic_parameters import data_exporter_simulation_thermodynamic_parameters as thermodynamic_parameters
from generator_pair_potential_particles_visualization import pair_potential_particles_visualization as visualizer
from generator_k_and_r_space_box import r_k_space as rks
from generator_wall_potential_values_visualization import wall_potential_values_visualization as wp_values
from generator_bulk_rho_mue_r_space import bulk_rho_mue_r_space as brm
from calculator_FMT_weights_1d_cartesian  import fmt_weights_1d as fm_weights






# File paths
input_file_name = 'executor_input.in'

try:
    interactions_parameters(input_file_name)
except Exception as e:
    print("Error exporting input data:", e)
    exit(0)
    
try:
    confinement_parameters(input_file_name)
except Exception as e:
    print("Error exporting input data:", e)
    exit(0)

try:
    thermodynamic_parameters(input_file_name)
except Exception as e:
    print("Error exporting input data:", e)
    exit(0)

try:
    visualizer()
except Exception as e:
    print("Error visualizing particles:", e)
    exit(0)
    
try:
    rks()
except Exception as e:
    print("Error generating k and r space:", e)
    exit(0)

try:
    wp_values()
except Exception as e:
    print("Error generating wall potential:", e)
    exit(0)

try: 
    brm()
except Exception as e:
    print("Error while generating the bulk rho and mue value ... \n")
    exit(0)

try:
    fm_weights()
except Exception as e:
    print ("Error while calculating FMT weights in k space ... \n")
    exit(0)

# Load space properties
json_file_interaction = "input_data_particles_interactions_parameters.json"
json_file_thermodynamics = "input_data_simulation_thermodynamic_parameters.json"

try:
    
    with open(json_file_thermodynamics, "r") as file:
        thermodynamics = json.load(file)["simulation_thermodynamic_parameters"]
    temperature, rho, iteration = thermodynamics["temperature"], thermodynamics["rho"], thermodynamics["iteration_max"]
    print(f"\n\n... thermodynamic properties loaded: Temperature = {temperature}, Rho = {rho}, Max Iterations = {iteration} ...\n\n")
except FileNotFoundError:
    print("Thermodynamic properties file not found.")
    exit()

try:
    with open(json_file_interaction, 'r') as file:
        data_interactions = json.load(file)["particles_interactions_parameters"]["interactions"]
    interaction_types = {}
    closest_distances = {}
    interaction_strength = {}
    cutoff_ranges = {}
    
    interaction_types["primary"] = {k: v["type"] for k, v in data_interactions["primary"].items()}
    closest_distances["primary"] = {k: v["sigma"] for k, v in data_interactions["primary"].items()}
    interaction_strength["primary"] = {k: v["epsilon"] for k, v in data_interactions["primary"].items()}
    cutoff_ranges["primary"] = {k: v["cutoff"] for k, v in data_interactions["primary"].items()}
    
    interaction_types["secondary"] = {k: v["type"] for k, v in data_interactions["secondary"].items()}
    closest_distances["secondary"] = {k: v["sigma"] for k, v in data_interactions["secondary"].items()}
    interaction_strength["secondary"] = {k: v["epsilon"] for k, v in data_interactions["secondary"].items()}
    cutoff_ranges["secondary"] = {k: v["cutoff"] for k, v in data_interactions["secondary"].items()}
    
    interaction_types["tertiary"] = {k: v["type"] for k, v in data_interactions["tertiary"].items()}
    closest_distances["tertiary"] = {k: v["sigma"] for k, v in data_interactions["tertiary"].items()}
    interaction_strength["tertiary"] = {k: v["epsilon"] for k, v in data_interactions["tertiary"].items()}
    cutoff_ranges["tertiary"] = {k: v["cutoff"] for k, v in data_interactions["tertiary"].items()}
    
    with open(json_file_interaction, 'r') as file:
        species = {k: v["rho_frac"] for k, v in json.load(file)["particles_interactions_parameters"]["species"].items()}
    print(f"\n... loaded {len(species)} species for simulation ...\n\n")
        
except FileNotFoundError:
    print("Interaction properties file not found.")
    exit()


print("\n") 
interaction_level = "primary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... primary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}\n")
          
print("\n")          
interaction_level = "secondary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... secondary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}\n")
          
print("\n")       
interaction_level = "tertiary"
# Display interactions for verification
for pair_type, interaction_type in interaction_types[interaction_level].items():
    print(f"... tertiary interaction Pair {pair_type}:  Type = {interaction_type}, Sigma = {closest_distances[interaction_level][pair_type]}, "
          f"Epsilon = {interaction_strength[interaction_level][pair_type]}, Cutoff = {cutoff_ranges[interaction_level][pair_type]}  \n")
          
print("\n\n\n")




import numpy as np

# Initialize empty lists for each column
x, y, z, rho_r, mue_r = [], [], [], [], []
# Read the text file for the bulk rho mue values ... 

for particle, rho in species.items():
    pid=0
    rho_ind= []
    mue_ind= []
    with open("supplied_data_bulk_mue_rho_r_space.txt", "r") as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue
            # Split the line into columns and convert them to floats
            columns = line.strip().split()
            if (pid == 0 ):
                x.append(float(columns[0]))
                y.append(float(columns[1]))
                z.append(float(columns[2]))
            
            li=[]
            chi=[]
            i=3+pid*2
            rho_ind.append(float(columns[i]))
            mue_ind.append (float(columns[i+1]))
    rho_r.append(rho_ind)
    mue_r.append(mue_ind)

    # Convert lists to numpy arrays
x = np.array(x)
y = np.array(y)
z = np.array(z)
rho_r = np.array(rho_r)
mue_r = np.array(mue_r)



# Read the text file for the fmt weights for the k space data 



kx, ky, kz = [], [], []  # Define kx, ky, kz before the loop
fmt_weights = {}  # Define fmt_weights to hold all species weights

for key, rho in species.items():
    fmt_weights_ind = []  # Initialize a list for individual species weights
    
    # Open file for the current species
    with open(f"supplied_data_weight_FMT_k_space_{key}.txt", "r") as file:
        for line in file:
            # Skip comment lines
            if line.startswith("#"):
                continue

            # Split the line into columns and convert them to floats
            columns = line.strip().split()
            kx.append(float(columns[0]))
            ky.append(float(columns[1]))
            kz.append(float(columns[2]))
            
            # Collect rho-related values for this species
            li = []

            li.append(float(columns[3]) + 1j*0 )
            li.append(float(columns[4]) + 1j*0 )
            li.append(float(columns[5]) + 1j*0 )
            li.append(float(columns[6]) + 1j*0 )
            li.append(0 + 1j*float(columns[7]) )
            li.append(0 + 1j*float(columns[8]) )
            
            fmt_weights_ind.append(li)
    
    fmt_weights[key] = fmt_weights_ind # Append the individual weights list to fmt_weights


print ("\n\n... supplied data has been imported successfully ...\n\n\n")

v_ext=[]
with open(f"supplied_data_walls_potential_r_space.txt", "r") as file:
    for line in file:
        # Skip comment lines
        if line.startswith("#"):
            continue

        # Split the line into columns and convert them to floats
        columns = line.strip().split()
        v_ext.append(float(columns[3]))
        

# this is the main regions for the calculation for the 1 d walls confinement DFT simulation ...

nufft = NUFFT()


i = 0
j = 0

iteration = 0
rho_r_initial = np.array(rho_r)
rho_r_current = np.array(rho_r)

piee = np.pi

threshold = 0.001
alpha = 0.4
while (iteration < iteration_max):
    
    rho_r_initial = rho_r_current 
    
    # energy filtration for the hard core fmt terms ...
    rho_alpha_r={}
    fmt_flag = 0
    pid = 0
    for particle1, rho1 in species.items():
        for particle2, rho2 in species.items():
            pair1 = particle1.strip() + particle2.strip()
            pair2 = particle2.strip() + particle1.strip()
            interaction_level = "primary" 
            if pair1 in interactions_types:
                pair = pair1
            elif pair2 in interactions_types:
                pair = pair2
            if interaction_types[interaction_level] == "hc" and particle1 == particle2:
                fmt_flag = 1
                rho_k_ind = np.fft.fft(rho_r_current[pid])
                omega_rho_k = []
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][0] * rho_k_ind[i] 
                    li.append(product)
                omega_rho_k.append(li)
                
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][1] * rho_k_ind[i] 
                    li.append(product)
                omega_rho_k.append(li)
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][2] * rho_k_ind[i]
                    li.append(product)
                omega_rho_k.append(li)   
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][3] * rho_k_ind[i] 
                    li.append(product)
                omega_rho_k.append(li)   
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][4] * rho_k_ind[i]
                    li.append(product)
                omega_rho_k.append(li)   
                li=[]
                for i in range(len(rho_k_ind)):
                    product =  fmt_weights[particle1][i][5] * rho_k_ind[i]
                    li.append(product)
                omega_rho_k.append(li)    
                
                
                omega_rho_k= np.array(omega_rho_k)
                rho_alpha_r_ind=[]
                for i in range(6):
                    rho_alpha= np.fft.ifft(omega_rho_k[i])
                    rho_alpha_r_ind.append(rho_alpha)
                
                rho_alpha_r[particle1] = np.array(rho_alpha_r_ind)
        pid = pid + 1         
    if (fmt_flag == 1):
        rho_alphas = []
        for j in range(6):
            li=[]
            for i in range(len(rho_r)):
                sum_species_rho_alpha = 0 + 0j
                for keys, values in rho_alpha_r.items():
                    sum_species_rho_alpha += float(values[j][i])
                li.append(sum_species_rho_alpha)
            rho_alphas.append(li)
        rho_alpha=np.array(rho_alphas)
        
        
        dphi=[]
        
        li=[]
        for i in range(len(rho_r)):
            dphi0 = -math.log(1.0 - np.real(rho_alpha[3][k]))
            li.append(dphi0)
        dphi.append(li)
        
        li=[]
        for i in range(len(rho_r)):
            dphi1 = np.real(rho_alpha[2][k]) / (1.0 - np.real(rho_alpha[3][k]))
            li.append(dphi1)
        dphi.append(li)
        
        li=[]
        for i in range(len(rho_r)):
            dphi2 = (np.real(rho_alpha[1][k]) / (1.0 - np.real(rho_alpha[3][k])) + (1.0 / (8.0 * piee)) * (np.real(rho_alpha[2][k]) ** 2 - np.real(rho_alpha[4][k]) ** 2) /((1.0 - np.real(rho_alpha[3][k])) ** 2))
            li.append(dphi2)
        dphi.append(li)
        
        li=[]
        for i in range(len(rho_r)):
            dphi3 = (np.real(rho_alpha[0][k]) / (1.0 - np.real(rho_alpha[3][k])) + (np.real(rho_alpha[1][k]) * np.real(rho_alpha[2][k]) - np.real(rho_alpha[4][k]) * np.real(rho_alpha[5][k]))/((1.0 - np.real(rho_alpha[3][k])) ** 2))
            dphi3 += ((1.0 / (12.0 * piee)) * ( np.real(rho_alpha[2][k]) ** 3 - 3.0 * np.real(rho_alpha[5][k])  ** 2)/((1.0 - np.real(rho_alpha[3][k]) ) ** 3))
            li.append(dphi3)
        dphi.append(li)
        
        li=[]
        for i in range(len(rho_r)):
            
            dphi4 = np.real(rho_alpha[5][k])/ (1.0 - np.real(rho_alpha[3][k]) )
            li.append(dphi4)
        dphi.append(li)
        
        li=[]
        for i in range(len(rho_r)):

            dphi5 = (np.real(rho_alpha[4][k])/ (1.0 - np.real(rho_alpha[3][k]) ) + 0.25 * (np.real(rho_alpha[2][k])  * np.real(rho_alpha[5][k]) ) / ((1.0 - np.real(rho_alpha[3][k])) ** 2))
            li.append(dphi5)
        dphi.append(li)
        
        dphi =  np.array(dphi)
                
        dphi_k=[]        
        for i in range(6):
            dphi_k_alpha = np.fft.fft(dphi[i])
            dphi_k.append(dphi_k)
        dphi_k = np.array(dphi_k)
        # dphi/dalpha is being calculated which is same for all the species ... 
        
        
    # energy filtration for the hard core fmt terms had already been performed ... thank you so much ...
    # now the summing of df/drho starts for all the individual component present in the system ...
    
    
    df_external = {} 
    diff_grand = 0.0
    diff_ind_grand=[]
    for particle1, rho1 in species.items():
    
        df_ext_ind = []
        for i in range(len(rho_k_ind)):
            df_ext_ind.append(0.0)
        for particle2, rho2 in species.items():
            pair1 = particle1.strip() + particle2.strip()
            pair2 = particle2.strip() + particle1.strip()
            interaction_level = "primary" 
            if pair1 in interactions_types:
                pair = pair1
            elif pair2 in interactions_types:
                pair = pair2
            if interaction_types[interaction_level] == "hc" and particle1 == particle2:
               
                rho_k_ind = np.fft.fft(rho_r_current[pid])
                omega_dphi_k = []
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][0] * dphi[0][i]
                    li.append(product)
                omega_dphi_k.append(li)
                
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][1] * dphi[1][i]
                    li.append(product)
                omega_dphi_k.append(li)
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][2] * dphi[2][i]
                    li.append(product)
                omega_dphi_k.append(li)   
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][3] * dphi[3][i]
                    li.append(product)
                omega_dphi_k.append(li)   
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][4] * dphi[4][i]
                    li.append(product)
                omega_dphi_k.append(li)   
                li=[]
                for i in range(len(rho_k_ind)):
                    product = fmt_weights[particle1][i][5] * dphi[5][i]
                    li.append(product)
                omega_dphi_k.append(li)    
                
                
                omega_dphi_k= np.array(omega_dphi_k)
                dphi_dalpha_r_ind=[]
                for i in range(6):
                    dphi_rho_alpha= np.fft.ifft(omega_dphi_k[i])
                    dphi_dalpha_r_ind.append(dphi_rho_alpha)
                
                
                for i in range(len(rho_k_ind)):
                    sum_alpha = 0
                    for j in range(6):
                        sum_alpha += dphi_dalpha_r_ind[j][i]
                    df_ext_ind[i] += sum_alpha     
            elif interaction_types[interaction_level] in ["wca", "gs"] :
                print("... it can be implemented easily however had not been done yet ...\n\n")
        df_ext[particle1] = np.array(df_ext_ind)
        
        for i in range(len(rho_k_ind)):
            rho_r_current[pid][i] = alpha * (np.exp( df_ext_ind[i] + mue_r[pid][i] + v_ext[pid][i] )) + (1-alpha) * rho_r_initial[pid][i]
        
        diff = 0.0
        for i in range(len(rho_k_ind)):
            diff += (rho_r_currnt[pid][i] - rho_r_initial[pid][i])**2.0
        diff = diff/len(rho_k_ind)
        diff_ind_grand.append(diff)
        
        diff_grand= diff_grand + diff;
        pid = pid + 1
                

    diff_grand = diff_grand/pid

    if (diff_grand < threshold):
        print ("\n\n... density has been minimized thank you so much for your patience ...\n\n\n") 
    
















    
    
# testing region ... 
'''
import numpy as np
import matplotlib.pyplot as plt

# Example r_space values from 0 to a certain cutoff (e.g., 10)
r_space = np.linspace(0, 10, 1000)  # 1000 points in r-space
r_space_data = np.exp(-r_space)  # Sample data, e.g., a decaying exponential

# Perform FFT
k_space_data = np.fft.fft(r_space_data)

# Get the k-space values, assuming spacing between r-space points is uniform
dr = r_space[1] - r_space[0]  # calculate spacing in r-space
k_values = np.fft.fftfreq(len(r_space), d=dr) * 2 * np.pi  # Convert to angular frequency



plt.close('all')
plt.clf()  # Clears current figure for a new plot
# Plot r-space data
plt.figure(figsize=(12, 5))
plt.subplot(1, 2, 1)
plt.plot(r_space, r_space_data, label='r-space data')
plt.xlabel('r')
plt.ylabel('Amplitude')
plt.title('Real-space Data')
plt.legend()

# Plot k-space data (magnitude)
plt.subplot(1, 2, 2)
plt.plot(np.fft.fftshift(k_values), np.fft.fftshift(np.abs(k_space_data)), label='k-space data')
plt.xlabel('k')
plt.ylabel('Magnitude')
plt.title('k-space Data')
plt.legend()

plt.tight_layout()
plt.show()
'''
