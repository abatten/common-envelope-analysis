from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import numpy as np
import matplotlib.pyplot as plt
import ConfigParser

import cemodules.cefunctions as cef
import cemodules.constants as const

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
    root_dir = config.get("Common", "root_dir") 
    #  Directories to exclude from reading
    exclude_dir = config.get("Common", "exclude_dir")
    #  Directory to save the output file
    plot_dir = config.get("Common", "plot_dir") 
    #  Where to begin in the directory list: Usually 0
    initial_path = config.getint("Common", "initial_path")
    #  Where to end in the directory list + 1 (To work with python arrays)
    final_path_plus_one = config.getint("Common",
                                        "final_path_plus_one")
    #  The number of particles in the simulation
    particle_number = config.getint("Common", "particle_number")
    #  The name of the output text file
    output_file_name = config.get("Gravodrag", "output_file_name")
    #  Append to the bottom of the text file or create new?
    output_file_append = config.getboolean("Gravodrag", 
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

    return(root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one, 
           particle_number, output_file_name, output_file_append)


def open_file(file_name, append, particle_number):
    """
    Opens file and writes the header information.

    Returns:
        output_file
    """

    if (append == False):  #  Overwrite
        output_file = open(file_name, "w")
        #  Create first part of header
        header = "Time (yr) Cycle(#)"

        density_dict = {}
        rel_vel_dict = {}
        gravodrag_dict = {}
        accr_rad_dict = {}
        sound_speed_dict = {}
        velocity_dict = {}
        comp_vel_dict = {}

        #  Adjust the number of columns for number of particles
        for i in range(particle_number - 1):  # One less! Don't want the core
            density_dict[str(i+1)] = ("%s_%s_%s" % 
                                     ("Particle", str(i+1),"Density (g/cm^3)"))
            rel_vel_dict[str(i+1)] = ("%s_%s_%s" %
                                     ("Particle", str(i+1),"Rel_Vel"))
            gravodrag_dict[str(i+1)] = ("%s_%s_%s" %
                                       ("Particle", str(i+1),"Gravodrag (g cm/s^2)"))
            accr_rad_dict[str(i+1)] = ("%s_%s_%s" %
                                      ("Particle", str(i+1),"Accretion_Rad (cm)"))
            sound_speed_dict[str(i+1)] = ("%s_%s_%s" %
                                         ("Particle", str(i+1),"Sound_Speed"))
            velocity_dict[str(i+1)] = ("%s_%s_%s" %
                                      ("Particle", str(i+1),"Gas_Vel"))
            comp_vel_dict[str(i+1)] = ("%s_%s_%s" %
                                      ("Particle", str(i+1),"Vel"))

            header = ",".join([header, gravodrag_dict[str(i+1)], 
                                density_dict[str(i+1)], 
                                accr_rad_dict[str(i+1)],
                                str(sound_speed_dict[str(i+1)]) + " (cm/s)", 
                                str(rel_vel_dict[str(i+1)]) + "_X (cm/s)",
                                str(rel_vel_dict[str(i+1)]) + "_Y (cm/s)",
                                str(rel_vel_dict[str(i+1)]) + "_Z (cm/s)",
                                str(velocity_dict[str(i+1)]) + "_X (cm/s)",
                                str(velocity_dict[str(i+1)]) + "_Y (cm/s)",
                                str(velocity_dict[str(i+1)]) + "_Z (cm/s)",
                                str(comp_vel_dict[str(i+1)]) + "_X (cm/s)",
                                str(comp_vel_dict[str(i+1)]) + "_Y (cm/s)",
                                str(comp_vel_dict[str(i+1)]) + "_Z (cm/s)"])

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
    lu = pf.parameters["LengthUnits"]
    tu = pf.parameters["TimeUnits"]
    mu = pf.parameters["MassUnits"]
    
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

    #  List of unordered particle masses
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

    #  Check that there is the correct number of particles
    comp_coords = {}
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
            coord = [ce['particle_position_x'][i] * lu,
                           ce['particle_position_y'][i] * lu,
                           ce['particle_position_z'][i] * lu]
            #  Assign values to dictionary
            comp_coords[pdex[i]] = coord

    #  Gravodrag = xi * pi * accretion_rad**2 * density * relative_vel**3
    density_in_cell = {}
    velocity_in_cell = {}
    comp_velocity = {}
    relative_velocity = {}
    accretion_radius = {}
    sound_speed_in_cell = {}
    comp_mass = {}
    gravodrag = {}

    for i in range(particle_number):
        if i != prim_index:
            #  Convert to code coordinates
            code_coords = [x / lu for x in comp_coords[pdex[i]]]

            #  PART 1: AVERAGE DENSITY
            #  Find the density at the position of the companion 
            #  Returns an array so take the [0] element
            density = pf.h.find_field_value_at_point(['Density'], 
                                                     code_coords)[0]
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
            relative_velocity[pdex[i]] = [x - y for x,y in 
                                          zip(part_velocity, gas_velocity)]

            #  PART 3: ACCRETION RADIUS
            #  accretion_radius = 2 * G * Mass / rel_vel**2 + sound_speed**2
            
            #  Find companion mass            
            #  const.YMSUN converts mass from solar masses to g
            comp_mass[pdex[i]] = particle_masses[i] * const.YMSUN

            #  Find the sound speed in the cell
            #  Returns an array so take the [0] element
            sound_speed = pf.h.find_field_value_at_point("SoundSpeed", 
                                                          code_coords)[0]
            sound_speed_in_cell[pdex[i]] = sound_speed

            #  Calculate magnitude of the relative velocity
            rel_vel_mag = np.linalg.norm(relative_velocity[pdex[i]])
            #  Calculate the accretion radius
            accretion_radius[pdex[i]] = (2.0 * const.YG * comp_mass[pdex[i]] /
                                        (rel_vel_mag**2.0 
                                        + sound_speed_in_cell[pdex[i]]))

            #  PART 4: XI
            #  Assume that xi = 1
            xi = 1

            #  PART 5: GRAVODRAG
            #  Gravodrag = xi * pi * accretion_rad**2 * density * rel_vel**3
            gravodrag[pdex[i]] = (xi * np.pi * accretion_radius[pdex[i]]**2.0 *
                                  density_in_cell[pdex[i]] * 
                                  rel_vel_mag**3.0)

    print("Accretion Radius (cm):", accretion_radius)
    print("Comp Mass (g)" , comp_mass)
    print("Density in Cell (g/cm^3): ", density_in_cell)
    print("Gas Velocity (cm/s): ", velocity_in_cell) 
    print("Comp Velocity (cm/s): ", comp_velocity)
    print("Relative Velocity (cm/s): ", relative_velocity)
    print("Sound Speed in Cell (cm/s): ", sound_speed_in_cell)
    print("Gravodrag (g cm/s^2)", gravodrag)

    data = str(current_time) + " " + str(current_cycle) + " "
   
    for i in range(particle_number):
        if i != prim_index:
            data += (str(gravodrag[pdex[i]]) + " " 
                    + str(density_in_cell[pdex[i]]) + " "
                    + str(accretion_radius[pdex[i]]) + " "
                    + str(sound_speed_in_cell[pdex[i]]) + " "
                    + str(relative_velocity[pdex[i]][0]) + " "
                    + str(relative_velocity[pdex[i]][1]) + " "
                    + str(relative_velocity[pdex[i]][2]) + " "
                    + str(velocity_in_cell[pdex[i]][0]) + " "
                    + str(velocity_in_cell[pdex[i]][1]) + " "
                    + str(velocity_in_cell[pdex[i]][2]) + " "
                    + str(comp_velocity[pdex[i]][0]) + " "
                    + str(comp_velocity[pdex[i]][1]) + " "
                    + str(comp_velocity[pdex[i]][2]) + " ")

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
