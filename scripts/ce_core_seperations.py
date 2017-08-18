from __future__ import absolute_import, division, print_function

import yt.mods as yt

import numpy as np
import sys
import ConfigParser

import cemodules.cefunctions as cef


def read_inlist(ipath):
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split('/')[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, 'r'))

    # Read in the config file
    root_dir = config.get("Common", "root_dir")
    exclude_dir = config.get("Common", "exclude_dir")
    plot_dir = config.get("Common", "plot_dir")
    initial_path = config.getint("Common", "initial_path")
    final_path_plus_one = config.getint("Common", "final_path_plus_one")
    output_file_name = config.get("Seperation", "output_file_name")
    particle_number = config.getint("Common", "particle_number")

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
    output_file = open(file_name, 'w')

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
    pf = yt.load(directory)

    #  Gets the length, time and mass units used in the current simulation
    lu = pf.parameters['LengthUnits']
    tu = pf.parameters['TimeUnits']

    current_cycle = pf.parameters['InitialCycleNumber']
    current_time = pf.current_time

    #  Adds the whole data as an object called ce.
    #  It is an array and you can acces the data by knowing their
    #  name through: ce["Dataname"]:
    ce = pf.h.all_data()

    if (index == initial_path):
        current_cycle = 0
    else:
        current_cycle = pf.parameters['InitialCycleNumber']

    yr = 365.35 * 24 * 60 * 60
    current_time = pf.current_time / yr

    # Find which index corresponds to the primary star (the largest mass).
    particle_masses = ce["ParticleMassMsun"]
    prim_mass = np.max(particle_masses)
    for i in range(len(particle_masses)):
        if particle_masses[i] == prim_mass:
            prim_index = i
            break
        else:
            pass

    print("Primary Index: ", prim_index)

    # Coordinates of the Primary Star
    primary_coords = [ce['particle_position_x'][prim_index] * lu,
                      ce['particle_position_y'][prim_index] * lu,
                      ce['particle_position_z'][prim_index] * lu]

    sep = {}
    for i in range(particle_number):
        #  Get list of particle indicies to track individual particles
        pdex = ce['particle_index']

        if i != prim_index:
            companion_coords = [ce['particle_position_x'][i] * lu,
                                ce['particle_position_y'][i] * lu,
                                ce['particle_position_z'][i] * lu]

            seperation = cef.distance(primary_coords, companion_coords)
            sep[pdex[i]] = seperation
            print("Companion " + str(pdex[i]), seperation)

    row = str(current_time * tu) + " " + str(current_cycle)
    for i in range(len(sep.items())):
        row = " ".join([row, str(sep.items()[i][1])])

    output_file.write(row + "\n")


if __name__ == "__main__":
    yt.mylog.disabled = True

    (root_dir, exclude_dir,
     plot_dir, initial_path,
     final_path_plus_one, output_file_name,
     particle_number) = read_inlist(sys.argv[1])

    # Sort the root directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, particle_number)

    for index in range(initial_path, final_path_plus_one):
        seperations(root_dir_list[index], index, output_file, particle_number)
