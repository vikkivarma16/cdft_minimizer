# density functional minimizer/input data generator

# this is part of the code prints all the input data in the JSON format so that it could be read and executed by the different modules


# this code is also a gateway to block the code and stop it from running for a given configuration which is not developed yet ...





def data_exporter(input_file):

    import json
    # Initialize dictionaries to hold the data
    # number of dimension and confinement types
    
    
    
    # a submodule to check whether a value is convertible or not 
    
    def is_convertible_to_float(value):
        try:
            float(value)  # Try to convert the string to a float
            return True
        except ValueError:
            return False
    
    
    
    
    
    space_properties={}
    
    # confinement type specific position space simulation properties
    
    box_properties = {}
    sphere_properties={}
    cylinder_properties={}
    
    
    # type of particles in the system
    species_data = {}
    
    
    # interaction among the particles in the system
    interactions_data = {}
    interactions_data ["primary"] = {}
    interactions_data ["secondary"] = {}
    interactions_data ["tertiary"] = {}
    
    
    # thermodynamic parameters defined for the system
    thermodynamic_properties = {}

    
    with open(input_file, 'r') as file:
        lines = file.readlines()

    # Parsing the input file

    # Setting up the default box properties
    space_properties["dimension"]=1
    space_properties["confinement"]="pbox"
   
    box_properties["box_length"] = [10.0, 10.0, 10.0]
    box_properties["box_points"] = [1000, 1000, 1000]
    
    # Default species data
    species = "a"
    species_data[species] = {"rho_frac": 1.0}
    
    
    for key in ["temperature", "rho", "iteration_max"]:
                thermodynamic_properties[key] = "NA"
    
    
    for line in lines:
        line = line.strip()  # Remove leading and trailing whitespace
        if not line or line.startswith("#"):
            continue  # Skip empty lines and comments

        # split the line into key and value
        if "=" in line and "interaction" not in line:
            key, value = line.split("=", 1)
            key = key.strip()
            value = value.strip()  # Remove quotes around strings

            if key == "space_dimension":
                space_properties["dimension"]= int(value)
                
        # setting up default confinement properties depending upon the dimension of the system
                
                space_properties["confinement"]="pbox"
                box_properties["box_length"]=[10.0, 0.0, 0.0]
                box_properties["box_points"]=[0.0, 0.0, 0.0]
                for i in range(int(value)):
                    box_properties["box_length"][i]=10.0
            
            elif key == "space_confinement":
                space_properties["confinement"] = value
                
            elif key == "box_extension":
                box_properties["box_length"]=[10, 10, 10]
                box_properties["box_length"] = [float(x) for x in value.split(",")]
            elif key == "box_points":
                box_properties["box_points"]=[100, 100, 100]
                
                box_properties["box_points"]= [float(x) for x in value.split(",")]
                
            elif key ==  "sphere_extension":  
                sphere_properties["radius_of_spherical_confinement"]= [5.0]
                sphere_properties["radius_of_spherical_confinement"]= value
            elif key ==  "sphere_points":
                sphere_properties["sphere_points"]= [100]
                sphere_properties["sphere_points"]= value
                

            elif key == "species":
                species_names = value.split(",")
                
                rho_frac = 1.0 / len(species_names)
                for name in species_names:
                    species_data[name.strip()] = {"rho_frac": rho_frac}  # Initialize species with rho_frac
                
                # Initialize default interaction data for species
                i=0
                j=0
                flag={}
                for name1 in species_names:
                    j=0
                    for name2 in species_names:
                        if j>=i :
                            interaction = name1.strip() + name2.strip()
                            flag[interaction]=[0,  0,  0]
                            interactions_data["primary"][interaction] = {
                                "type": "gs",
                                "sigma": 1.1,
                                "cutoff": 3.4,
                                "epsilon": 2.0
                            }
                        j = j+1
                    i = i+1
                interactions_data["secondary"] = {}
                interactions_data["tertiary"] = {}
                
            elif key == "species_fraction":
                rho_values = value.split(",")
                
                total_rho=0.0
                for i, name in enumerate(species_data.keys()):
                    total_rho = total_rho + float(rho_values[i].strip())
    
                for i, name in enumerate(species_data.keys()):
                    species_data[name]["rho_frac"] = float(rho_values[i].strip()) / total_rho

            elif key in ["temperature", "rho", "iteration_max"]:
                thermodynamic_properties[key] = float(value)

        # Handle interactions
        elif "interaction" in line:
            
            interaction, dummy = line.split("=", 1)
            
            
            dummy, properties_str = line.split(":", 1)
            
            dummy, pair= interaction.split(":", 1) 
            
            # Remove leading/trailing spaces from interaction type
            pair = pair.strip()

            # Step 2: Split the properties part (after the '=') by commas to get each key-value pair
            properties_list = properties_str.split(",")

            # Initialize a dictionary to store the property values
            properties_dict = {}

            # Step 3: Process each key-value pair and store them in a dictionary
            
            
    
            if flag[pair][0]==0:   
                flag[pair][0]=1
                for prop in properties_list:
                
                   
                    param_key, param_value = prop.split("=", 1)
                    param_key = param_key.strip()
                    param_value = param_value.strip()
                    

                    # Only replace default values with input values
                
                    for interaction in interactions_data["primary"]:
                        if pair in interaction:
                            
                            if param_key != pair:
                                if is_convertible_to_float(param_value)==True:
                                    interactions_data["primary"][pair][param_key] = float(param_value)
                                    
                                else:
                                    interactions_data["primary"][pair][param_key] = param_value
                                    
                                if "type" == "hc":
                                    interactions_data["primary"][pair]["cutoff"]= interaction_data[pair]["sigma"]
                                    
                            elif param_key == pair :
                            
                                interactions_data["primary"][pair]["type"] = param_value
            elif flag[pair][1]==0:
                
                temp={}
                flag[pair][1]=1
                
                
                for prop in properties_list:
                
                    param_key, param_value = prop.split("=", 1)
                    param_key = param_key.strip()
                    param_value = param_value.strip()
                    
                    # Only replace default values with input values
                
                    if param_key != pair:
                        if is_convertible_to_float(param_value)==True:
                            temp[param_key] = float(param_value)
                            
                        else:
                            temp[param_key] = param_value
                            
                        if "type" == "hc":
                            temp["cutoff"]= temp["sigma"]
                            
                    elif param_key == pair :

                        temp["type"] = param_value
                
                interactions_data["secondary"][pair]=temp
                
                
            elif flag[pair][2]==0:
                temp={}
                flag[pair][2]=1
                for prop in properties_list:
                
                   
                    param_key, param_value = prop.split("=", 1)
                    param_key = param_key.strip()
                    param_value = param_value.strip()
                    

                    # Only replace default values with input values
                
                    
                    if param_key != pair:
                        if is_convertible_to_float(param_value)==True:
                            temp[param_key] = float(param_value)
                            
                        else:
                            temp[param_key] = param_value
                            
                        if "type" == "hc":
                            temp["cutoff"]= temp["sigma"]
                            
                    elif param_key == pair :

                        temp["type"] = param_value
                
                interactions_data["tertiary"][pair]=temp
                
            else:
            
                print("Error: too many parameters have been defined for the interaction between the pair", pair, "which is not appropriate... please reduce and adjust the maximum number of interactions to 3 ")
                exit(0)
                
                
                
                        
            
            # Save to interactions_data

    # Save the parsed data to JSON files
    
    if space_properties["confinement"]=="pbox":
        with open('input_box_properties.json', 'w') as f:
            json.dump({"box_properties": box_properties}, f, indent=4)
    elif space_properties["confinement"]=="abox":
        print("\n\n properties has not been specified for the abox.. the code is still in development... Thank you so much...\n\n\n")
        exit(0)
    elif space_properties["confinement"]=="asphere":
        print("\n\n properties has not been specified for the abox.. the code is still in development... Thank you so much...\n\n\n")
        exit(0)
    elif space_properties["confinement"]=="psphere":
        print("\n\n periodicity is not applicable in the spherical simulation confinement... please chose the physical options... Thank you so much...\n\n\n")
        exit(0)

    with open('input_species_properties.json', 'w') as f:
        json.dump({"species": species_data}, f, indent=4)

    with open('input_interactions_properties.json', 'w') as f:
        json.dump({"interactions": interactions_data}, f, indent=4)

    if thermodynamic_properties["temperature"]== "NA" or thermodynamic_properties["rho"]== "NA" or thermodynamic_properties["iteration_max"]== "NA":
        print ("error: specify thermodynamic properties like rho=value, temperature= value, iteration_max= value, in row wise stacking format otherwise simulation will not run..")
        exit(0)
    else:
        with open('input_thermodynamic_properties.json', 'w') as f:
            json.dump({"thermodynamic_properties": thermodynamic_properties}, f, indent=4)
            
    with open('input_space_properties.json', 'w') as f:
        json.dump({"space_properties": space_properties}, f , indent=4)

    print("\n\n\n  .............. the data has been successfully formatted and exported to the JSON files ............  \n\n\n")
    
    return 0
# Specify the input file path
#input_file_path = 'executor_input.in'  # Replace with your actual input file name
#dummy = data_exporter(input_file_path)




# this part of the code run the various segments of the code before running the actual minimizer codes.....



