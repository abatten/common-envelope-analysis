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
    inlist_name = ipath.split("/")[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, "r"))

    # Read in Config File
    root_dir = config.get("Common", "root_dir")
    exclude_dir = config.get("Common", "exclude_dir")
    plot_dir = config.get("Common", "plot_dir")
    initial_path = config.getint("Common", "initial_path")
    final_path_plus_one = config.getint("Common", "final_path_plus_one")
    output_file_name = config.get("KeplerVel", "output_file_name")
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
    output_file = open(file_name, "w")

    header = "Time (yr),  Cycle(#)"

    # Create header based on the number of particles
    vel_dict = {}
    kep_dict = {}
    for i in range(num_particle-1):
        vel_dict[str(i+1)] = "Vel_%s" % (i+1) + " (cm/s)"
        kep_dict[str(i+1)] = "KepVel_%s" % (i+1) + " (cm/s)"

        header = ",  ".join([header, vel_dict[str(i+1)], kep_dict[str(i+1)]])
    output_file.write(header + "\n")

    return output_file

def velocities(directory, index, output_file, particle_number):
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->") 
    print("READ ROOT DIRECTORY " + str_index + ":", directory)
    print("<-------------->")
    pf = yt.load(directory)

    lu = pf.parameters["LengthUnits"]
    tu = pf.parameters["TimeUnits"]

    current_cycle = pf.parameters["InitialCycleNumber"]
    year = 365.25 * 24 * 60 * 60
    current_time = pf.current_time / year

    ce = pf.h.all_data()

    #particle_masses = ce["ParticleMassSun"]
    #prim_mass = np.max(particle_masses)
    #for i in range(len(particle_masses)):
    #    if particle_masses[i] == prim_mass:
    #        prim_index = i
    #        break
    #    else:
    #        pass
    prim_index = cef.primary_index(ce)  # Finds the index of the primary core

    primary_coords = 

    
