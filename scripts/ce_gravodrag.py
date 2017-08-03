from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

import cemodules.cefunctions as cef

mylog.disabled = True  

##############################################################################
##############################################################################
##############################################################################

def read_inlist(ipath):
    """
    Reads the contents of the inlist file and returns the values of the
    appropriate variables. Returns the values in the following order:

    root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one,
    particle_number, output_file_name, output_file_append
    """

    print(" ")
    print("<-------------->")  
    print("READING INLIST FILE") 
    print("<-------------->") 
    inlist_name = ipath.split("/")[-1] 
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, "r"))
 
    #  The directory of all the enzo outputs
    root_dir = config.get("Common Section", "root_dir") 
    #  Directories to exclude from reading
    exclude_dir = config.get("Common Section", "exclude_dir")
    #  Directory to save the output file
    plot_dir = config.get("Common Section", "plot_dir") 
    #  Where to begin in the directory list: Usually 0
    initial_path = config.getint("Common Section", "initial_path")
    #  Where to end in the directory list + 1 (To work with python arrays)
    final_path_plus_one = config.getint("Common Section", 
                                        "final_path_plus_one")
    #  The number of particles in the simulation
    particle_number =config.getint("Common Section", "particle_number")
    #  The name of the output text file
    output_file_name = config.get("Gravodrag Section", "output_file_name")
    #  Append to the bottom of the text file or create new?
    output_file_append = config.getboolean("Gravodrag Section", 
                                           "output_file_append") 
    
    print("INLIST FILE: " + inlist_name) 
    print("ROOT DIRECTORY: " + str(root_dir)) 
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("PARTICLE NUMBER: " + str(particle_number)) 
    print("OUTPUT FILE NAME: " + str(output_file_name))
    print("OUTPUT FILE APPEND: " + str(output_file_append))

    return(root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one,           particle_number, output_file_name, output_file_append)


def open_file(file_name, append, particle_number):
    """
    Opens file and writes the header information.

    Returns:
        output_file
    """

    if (append == False):  #  Overwrite
        output_file = open(file_name, "w")
        #  Create first part of header
        header = "Time (yr), Cycle(#)"

        density_dict = {}
        rel_vel_dict = {}
        gravodrag_dict = {}
        
        #  Adjust the number of columns for number of particles
        for i in range(particle_number - 1):  # One less! Don't want the core
            density_dict[str(i+1)] = ("%s_%s_%s" % 
                                     ("Particle", str(i+1),"Density"))
            rel_vel_dict[str(i+1)] = ("%s_%s_%s" %
                                     ("Particle", str(i+1),"Rel_Vel"))
            gravodrag_dict[str(i+1)] = ("%s_%s_%s" %
                                       ("Particle", str(i+1),"Gravodrag"))

            header = ", ".join([header, density_dict[str(i+1)], 
                                rel_vel_dict[str(i+1)], 
                                gravodrag_dict[str(i+1)]])

        #  Write the header line of the file
        output_file.write(header + "\n")
        return output_file


def ce_gravodrag(directory, index, outfile):
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->")
    print("READ ROOTDIRECTORY " + str_index + ":", directory)
    print("<-------------->")
    
    #  Load the profile pf
    pf = load(directory) 

    #  Get the length, time and mass units used
    length_unit1 = pf.parameters["LengthUnits"]
    time_unit1 = pf.parameters["TimeUnits"]
    mass_unit1 = pf.parameters["MassUnits"]
    
    #  Current cycle and current time will be written in the out_file
    current_time = pf.current_time                                    
    if (index == initial_path):
        current_cycle = 0
    else:
        current_cycle = pf.parameters["InitialCycleNumber"]

    #  Adds the whole data as an object called common_envelope.
    #  It is an array and you can acces the data by
    #  knowing their name through: ce["Dataname"]:
    ce = pf.h.all_data()

    #  List of particle masses
    particle_masses = ce["ParticleMassMsun"]

    #  The core of the primary has the largest mass
    prim_mass = np.max(particle_masses)

    #  Find the index of the primary
    for i in range(len(particle_masses)):
        if particle_masses[i] == prim_mass:
            prim_index = i
            break
        else:
            pass
    print("Primary Index: ", prim_index)

    #  Get the indicies of the particles. i.e track specific particles
    pdex = ce["particle_index"]
    comp_coords = {}
    comp_velocity = {}
    
    if len(pdex) != particle_number:
        print("ERROR")
        print("More particles in simulation than expecting!")
        print("Expecting: ", particle_number)
        print("Found: ", len(pdex))
        sys.exit(0)
    else:
        pass

    #  Get positions of companions
    for i in range(particle_number):
        if i != prim_index:  #  Ignore Primary as already done
            coord = [ce['particle_position_x'][i] * length_unit1,
                           ce['particle_position_y'][i] * length_unit1,
                           ce['particle_position_z'][i] * length_unit1]

            #  Assign values to dictionary
            comp_coords[pdex[i]] = coord

    # Gravodrag = xi * pi * accretion_rad**2 * density * relative_vel**3

    density_in_cell = {}
    velocity_in_cell = {}
    relative_velocity = {}
    accretion_radius = {}
    sound_speed_in_cell = {}
    comp_mass = {}

    for i in range(particle_number):
        if i != prim_index:
            #  Convert to code coordinates
            code_coords = [x/length_unit1 for x in comp_coords[pdex[i]]]

            #  PART 1: AVERAGE DENSITY
            #  Find the density at the position of the companion 
            density = pf.h.find_field_value_at_point(['Density'], code_coords)
            density_in_cell[pdex[i]] = density

            #  PART 2: AVERAGE RELATIVE VELOCITY
            #  Find the velocity of the gas in the companion cell
            gas_velocity = pf.h.find_field_value_at_point(["x-velocity",
                                                           "y-velocity",
                                                           "z-velocity"],
                                                           code_coords)
            part_velocity = [ce["particle_velocity_x"][pdex[i]],
                        ce["particle_velocity_y"][pdex[i]],
                        ce["particle_velocity_z"][pdex[i]]]

            velocity_in_cell[pdex[i]] = gas_velocity
            comp_velocity[pdex[i]] = part_velocity
            relative_velocity[pdex[i]] = [x-y for x,y in 
                                          zip(part_velocity, gas_velocity)]


            #  PART 3: ACCRETION RADIUS
            #  accretion_radius = 2 * G * Mass / rel_vel**2 + sound_speed**2
            
            #  Find the sound speed in the cell
            sound_speed = pf.h.find_field_value_at_point(["SoundSpeed"], 
                                                          code_coords)
            sound_speed_in_cell[pdex[i]] = sound_speed
            rel_vel_mag = np.linalg.norm(relative_velocity[pdex[i]])
            print(rel_vel_mag)            


    print("Density In Cell: ", density_in_cell)
    print("Gas Velocity: ", velocity_in_cell) 
    print("Comp Velocity: ", comp_velocity)
    print("Relative Velocity: ", relative_velocity)


    #  PART 4: XI
    
    data = str(current_time) + " " + str(current_cycle) + " "
    outfile.write(data)
    outfile.write("\n")

    return pf


if __name__ == "__main__":
    #  Read the inlist file and return the values of the variables
    if len(sys.argv) >= 2:
        inlist_path = sys.argv[1]
        (root_dir, exclude_dir, plot_dir, initial_path, 
         final_path_plus_one, particle_number, output_file_name, 
         output_file_append) = read_inlist(inlist_path)

    else:
        print("Inlist File Not Supplied!")
        sys.exit(0)

    #  Sort the root directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    #  Set the output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, output_file_append, 
                            particle_number)

    #  Calculate the gravodrag for every directory and write them
    for index in range(initial_path, final_path_plus_one):
        pf = ce_gravodrag(root_dir_list[index], index, output_file)

    output_file.close()

    print(" ")
    print("<-------------->")  
    print("FINISHED: CE_GRAVODRAG") 
    print("<-------------->") 
 
