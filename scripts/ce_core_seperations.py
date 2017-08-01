from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import ConfigParser

import cemodules.cefunctions as cef

mylog.disabled = True

#################################################################################
def read_inlist(ipath):
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split('/')[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, 'r'))

    # Read in the config file
    root_dir = config.get('Common Section', 'root_dir')
    exclude_dir = config.get('Common Section', 'exclude_dir')
    plot_dir = config.get('Common Section', 'plot_dir')
    initial_path = config.getint('Common Section', 'initial_path')
    final_path_plus_one = config.getint("Common Section", "final_path_plus_one")
    output_file_name = config.get("Seperation Section", "output_file_name")
    particle_number = config.getint("Common Section", "particle_number")

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("OUTPUT FILE NAME: " + str(output_file_name))
    print("PARTICLE NUMBER: " + str(particle_number))

    return (root_dir, exclude_dir, plot_dir, initial_path,
            final_path_plus_one, output_file_name, particle_number)


def open_file(file_name, num_particles):
    output_file = open(file_name, 'w' )

    header = "Time (yr), Cycle(#)"

    # Create header based on the number of particles
    dict = {}
    for i in range(num_particles-1):
       dict[str(i+1)] = "Seperation_%s_%s" % ("Companion", i+1) + " (cm)"
       header = ", ".join([header, dict[str(i+1)]])
    
    # Write the first line of information in the file
    output_file.write(header + "\n")

    return output_file


def seperations(directory, index, output_file, particle_number):

    str_index = cef.index2str(index)       

    print(" ")
    print("<-------------->")
    print("READ ROOTDIRECTORY " + str_index + ":", directory)
    print("<-------------->")
    pf = load(directory)
      
    #gets the length, time and mass units used in the current simulation    
    length_unit1 = pf.parameters['LengthUnits']
    time_unit1 = pf.parameters['TimeUnits']
    mass_unit1 = pf.parameters['MassUnits']

    current_cycle = pf.parameters['InitialCycleNumber']
    current_time = pf.current_time

    #adds the whole data as an object called common_envelope. It is an array and you can acces the data by knowing their name through: common_envelope["Dataname"]:
    common_envelope = pf.h.all_data()

    if (index == initial_path):
        current_cycle = 0
    else:
        current_cycle = pf.parameters['InitialCycleNumber']

    yr = 365.35 * 24 * 60 * 60
    current_time = pf.current_time / yr


    # Find which index corresponds to the primary star (the largest mass).
    particle_masses = common_envelope["ParticleMassMsun"]
    primary_mass = np.max(particle_masses)
    for i in range(len(particle_masses)):
        if particle_masses[i] == primary_mass:
            primary_index = i
            break
        else:
            pass    

    print("Primary Index: ", primary_index)

    # Coordinates of the Primary Star
    primary_coords = [common_envelope['particle_position_x'][primary_index] * length_unit1,
                      common_envelope['particle_position_y'][primary_index] * length_unit1,
                      common_envelope['particle_position_z'][primary_index] * length_unit1]

    sep = {}

    for i in range(particle_number):
        particle_indicies = common_envelope['particle_index']

        if i != primary_index:

            companion_coords = [common_envelope['particle_position_x'][i] * length_unit1,
                                common_envelope['particle_position_y'][i] * length_unit1,
                                common_envelope['particle_position_z'][i] * length_unit1]

            seperation = cef.distance(primary_coords, companion_coords)
            sep[particle_indicies[i]] = seperation
            print("Companion " + str(particle_indicies[i]), seperation)

    row = str(current_time * time_unit1) + " " +str (current_cycle)
    for i in range(len(sep.items())):
        row = " ".join([row, str(sep.items()[i][1])])


    output_file.write(row + "\n")


if __name__ == "__main__":

    (root_dir, exclude_dir, plot_dir, initial_path,
     final_path_plus_one, output_file_name, particle_number) = read_inlist(sys.argv[1])


    # Sort the root directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, particle_number)


    for index in range(initial_path, final_path_plus_one):
        seperations(root_dir_list[index], index, output_file, particle_number)

