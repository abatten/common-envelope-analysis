from yt.mods import *
import numpy as np
import math
import os
import time as timeclass
import ConfigParser

import cemodules.cefunctions as cef
# Time Check
initial_time = timeclass.time()

# Disable YT Outputs on screen (True = Yes, False = No)
yt_outputs = "True"
mylog.disabled = yt_outputs

def read_inlist(ipath):
    """
    Reads the inlist and then returns variables based
    on the input file. 

    Return order:
            root_dir, exclude_dir, plot_dir,
            initial_path, final_path_plus_one,
            particle_number, smoothing_length,
            output_file_name, output_file_append,
            use_smoothed_potential
    """
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split("/")[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, "r"))

    # Read in the config file
    root_dir = config.get("Common", "root_dir")
    exclude_dir = config.get("Commonn", "exclude_dir")
    plot_dir = config.get("Common", "plot_dir")
    initial_path = config.getint("Common", "initial_path")
    final_path_plus_one = config.getint("Common", "final_path_plus_one")
    particle_number = config.getint("Common", "particle_number")
    smoothing_length = config.getfloat("Common", "smoothing_length")
    
    output_file_name = config.get("Mass Loss", "output_file_name")
    output_file_append = config.getboolean("Mass Loss", "output_file_append")
    use_smoothed_potential = config.getboolean("Mass Loss", "use_smoothed_potential")

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("PARTICLE NUMBER: " + str(particle_number))
    print("SMOOTHING LENGTH: " + str(smoothing_length))
    print("OUTPUT FILE NAME: " + str(output_file_name))
    print("OUTPUT FILE APPEND: " + str(output_file_append))
    print("USE SMOOTHED POTENTIAL: " + str(use_smoothed_potential))

    return (root_dir, exclude_dir, plot_dir, 
            initial_path, final_path_plus_one,
            particle_number, smoothing_length, 
            output_file_name, output_file_append, 
            use_smoothed_potential)


# New yt field containing the indexes of the cells of the primary, 
# selected via a boundness criterion
def _BoundCells(field, data):            
    """
    New YT field containing the indexes of the cells of the primary
    selected via a boundness criterion.        
    """
    # Gets the length, time and mass units used in the current simulation    
    length_unit1 = data.pf.parameters["LengthUnits"]
    time_unit1 = data.pf.parameters["TimeUnits"]
    mass_unit1 = data.pf.parameters["MassUnits"]
        
    # Total Energy = Thermal Energy + Kinetic Energy + Potential Energy Gas to Gas + Potential Energy Particle to Gas
    # Step 1: Thermal Energy of Gas 
    current_Etherm_gas = data["ThermalEnergy"] * data["CellMass"]

    # Step 2: Kinetic Energy of Gas
    current_Ekin_gas = data["KineticEnergy"] * data["CellVolume"]
        
    # Step 3: Potential Energy Gas to Gas
    # The potential energy is computed only if the potential field is available
    if (data.pf.parameters["SelfGravity"] == 1):
            current_Epot_gas = data["Grav_Potential"]  * (length_unit1/time_unit1)**2.0 * data["CellMass"]
    else:
            current_Epot_gas = 0.0
        
    # Step 4: Potential Energy Particle to Gas
    # Particles (working only for the common envelope problem type)
    if (data.pf.parameters["ProblemType"]) == 41:
            # Whole box needed to find particles
            ce = data.pf.h.all_data()
             
            # Determine the size of the smallest cell for the smoothig length                
            smallest_cell_length = data.pf.h.get_smallest_dx() * length_unit1
            current_Epot_particle = dict()

            # Coordinates of the gas
            gas_coords = (data["x"], data["y"], data["z"])
            for i in range(particle_number):
                particle_coords = (ce["particle_position_x"][i],
                                   ce["particle_position_y"][i],
                                   ce["particle_position_z"][i])
                # Calculate the distance between teh gas and the particle
                particle_gas_distance = cef.distance(gas_coords, particle_coords, length_unit1) 
                # Smoothed Gravitational Potential
                # M Ruffert 1993
                current_Epot_particle[i] = cef.grav_pot(ce["ParticleMass"][i], data["CellMass"], 
                                                        particle_gas_distance, use_smoothed_potential, 
                                                        smoothing_length, smallest_cell_length)

            # Initialise the Potential Energy as 0 before the loop
            current_Epot_particle_to_gas = 0
            #for j in range(len(ce["particle_index"])):
            for j in range(particle_number):   
                current_Epot_particle_to_gas += current_Epot_particle[j]

    else:
        current_Epot_particle_to_gas = 0.0

    # Computes the total energy of the gas in each cell
    # Step 1 + Step 2 + Step 3 + Step 4
    Etotal = (current_Etherm_gas + current_Ekin_gas + 
             current_Epot_gas + current_Epot_particle_to_gas)

    # Locations of bound and unbound cells
    bound = Etotal < 0.0
    unbound = Etotal >= 0.0

    return(bound)

def ce_mass_loss(directory, index, outfile):
    pf = load(directory)
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->")
    print("READING ROOT DIRECTORY " + str_index + ": " + directory)
    print("<-------------->")

    # Gets the length, time and mass units used in the current simulation
    length_unit1 = pf.parameters["LengthUnits"]
    time_unit1 = pf.parameters["TimeUnits"]
    mass_unit1 = pf.parameters["MassUnits"]

    # Adds the whole data as an object called common_envelope
    # It is an array that can be access the data by knowing their 
    # name through: common_envelope["Dataname"]:
    ce = pf.h.all_data()
    
    # What is the current time?
    current_time = pf.current_time

    # Boolean Arrays of the bound and unbound cells
    bound_cells = ce["BoundCells"]
    unbound_cells = np.invert(bound_cells)

    # Mass Components
    current_mass_gas_total = np.sum(ce["CellMassMsun"])
    current_mass_gas_bound = np.sum(ce["CellMassMsun"][bound_cells])
    current_mass_gas_unbound = np.sum(ce["CellMassMsun"][unbound_cells])

    # Write to file
    outfile.write(str(current_time) + " " + str(current_mass_gas_total) + " " + str(current_mass_gas_bound) + " " + str(current_mass_gas_unbound))
    outfile.write("\n")

    return pf

def open_file(file_name, append):
    if append == False:  # Overwrite existing file
        output_file = open(file_name, "w")

        # Write the header information in the file
        output_file.write("Time (yr), Gas Mass Total (Msun), Gas Mass Bound (Msun), Gas Mass Unbound (Msun)")
        output_file.write("\n")
    elif append:  # Append to existing file
        output_file = open(output_file_name, "a")
    return output_file

if __name__ == "__main__":
    # Read the inlist file and return the values of the variables.
    if len(sys.argv) >= 2:  # Check for supplied inlist
        inlist_path = sys.argv[1]
        (root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one,
        particle_number, smoothing_length, output_file_name, 
        output_file_append, use_smoothed_potential) = read_inlist(inlist_path)
    else:
        # Bring up error if there is no inlist
        print("Inlist File Not Supplied!!!")
        sys.exit(0)

    # Adds the field to the available yt fields
    add_field("BoundCells", function=_BoundCells, units=r"Boolean array")

    # Sort the Root Directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, output_file_append)
    # Calculate 
    for index in range(initial_path, final_path_plus_one):
        pf = ce_mass_loss(root_dir_list[index], index, output_file)
  
    output_file.close()

    print(" ")
    print("<-------------->")
    print("FINISHED: CE_MASS_LOSS")
    print("<-------------->")
