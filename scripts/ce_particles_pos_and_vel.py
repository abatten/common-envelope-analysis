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
    output_file_name = config.get("Position-Velocity Section", "output_file_name")
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
    pos_dict = {}
    posvec = ["x", "y", "z"]
    vel_dict = {}
    velvec = ["vx", "vy", "vz"]
    for i in range(num_particles):
        for j in range(3):
            pos_dict[str(j+1)] = "%s_%s_%s" % ("Particle", str(i+1), posvec[j])
            vel_dict[str(j+1)] = "%s_%s_%s" % ("Particle", str(i+1), velvec[j])

            header = ", ".join([header, pos_dict[str(j+1)], vel_dict[str(j+1)]])

    # Write the first line of information in the file
    output_file.write(header + "\n")
    return output_file

def positions_velocities(directory, index, outfile):
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

    # Adds the whole data as an object called common_envelope. 
    # It is an array and you can acces the data by 
    # knowing their name through: ce["Dataname"]:
    ce = pf.h.all_data()

    if (index == initial_path):
        current_cycle = 0
    else:
        current_cycle = pf.parameters['InitialCycleNumber']

    particle_masses = ce['ParticleMassMsun']
    
    # The core of the primary has the largest mass
    primary_mass = np.max(particle_masses)

    for i in range(len(particle_masses)):
        if particle_masses[i] == primary_mass:
            primary_index = i
            break
        else:
            pass
    print("Primary Index: ", primary_index)

    primary_coords = [ce['particle_position_x'][primary_index] * length_unit1,
                      ce['particle_position_y'][primary_index] * length_unit1,
                      ce['particle_position_z'][primary_index] * length_unit1]

    primary_vel = [ce['particle_velocity_x'][primary_index] * length_unit1,
                   ce['particle_velocity_y'][primary_index] * length_unit1,
                   ce['particle_velocity_z'][primary_index] * length_unit1]

    particle_indicies = ce['particle_index']
    print(particle_indicies)

    coords = {}
    vels = {}
    for i in range(particle_number):
        if i != primary_index:
            companion_coords = [str(ce['particle_position_x'][i] * length_unit1),
                                str(ce['particle_position_y'][i] * length_unit1),
                                str(ce['particle_position_z'][i] * length_unit1)]

            companion_vels = [str(ce['particle_velocity_x'][i] * length_unit1),
                              str(ce['particle_velocity_y'][i] * length_unit1),
                              str(ce['particle_velocity_z'][i] * length_unit1)]
       

            coords[particle_indicies[i]] = companion_coords
            vels[particle_indicies[i]] = companion_vels

    data = str(current_time) + " " + str(current_cycle) + " "
    for i in range(len(primary_coords)):
        data += str(primary_coords[i]) + " " + str(primary_vel[i]) + " "

    for i in range(len(coords.items())):
        for j in range(len(coords.items()[i][1])):
            data += str(coords.items()[i][1][j]) + " " + str(vels.items()[i][1][j]) + " "

#    output_file.write(str(current_time)+" "+str(current_cycle)+" "+str(primary_coords[0])+" "
#                     +str(primary_vel[0])+" "+str(primary_coords[1])+" "+str(primary_vel[1])+" "
#                     +str(primary_coords[2])+" "+str(primary_vel[2])+" ")


   
#    for i in range(len(coords.items())):
#        print(coords.items()[i], vels.items()[i])
#        print(coords.items()[i][1][1:-1], vels.items()[i][1])
#        output_file.write(str(coords.items()[i][1]) + " " + str(vels.items()[i][1]) + " ")

    
#    for i in range(len(particle_indicies)):
#        if i != primary_index:
#            output_file.write(str(common_envelope["particle_position_x"][particle_indicies[i]]*length_unit1)+" "
#                              +str(common_envelope["particle_velocity_x"][particle_indicies[i]])+" " 
#                              +str(common_envelope["particle_position_y"][particle_indicies[i]]*length_unit1)+" "
#                              +str(common_envelope["particle_velocity_y"][particle_indicies[i]])+" "
#                              +str(common_envelope["particle_position_z"][particle_indicies[i]]*length_unit1)+" "
#                              +str(common_envelope["particle_velocity_z"][particle_indicies[i]])+" ")

    outfile.write(data)

    outfile.write("\n")
         

if __name__ == "__main__":

    (root_dir, exclude_dir, plot_dir, initial_path,
     final_path_plus_one, output_file_name, particle_number) = read_inlist(sys.argv[1])


    # Sort the root directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, particle_number)
    

    for index in range(initial_path, final_path_plus_one):
        positions_velocities(root_dir_list[index], index, output_file)
    output_file.close()

