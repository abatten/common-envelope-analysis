from __future__ import absolute_import, division, print_function

import yt.mods as yt

import numpy as np
import sys
import ConfigParser

import cemodules.cefunctions as cef
import cemodules.constants as CONST

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
    for i in range(num_particles-1):
        vel_dict[str(i+1)] = "Vel_%s" % (i+1) + " (cm/s)"
        kep_dict[str(i+1)] = "KepVel_%s" % (i+1) + " (cm/s)"

        header = ",  ".join([header, vel_dict[str(i+1)], kep_dict[str(i+1)]])
    output_file.write(header + "\n")

    return output_file


def kepler_velocity(seperation, mass):
    GRAVCONST = CONST.YG
    prim_mass = mass * CONST.YMSUN
    
    velocity = np.sqrt((GRAVCONST * prim_mass) / seperation)

    return velocity


def velocities(directory, index, output_file, particle_number):
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->") 
    print("READ ROOT DIRECTORY " + str_index + ":", directory)
    print("<-------------->")

    pf = yt.load(directory)

    lu = pf.parameters["LengthUnits"]  # Length Units
    tu = pf.parameters["TimeUnits"]    # Time Units

    current_cycle = pf.parameters["InitialCycleNumber"]
    year = 365.25 * 24 * 60 * 60
    current_time = pf.current_time * tu / year

    ce = pf.h.all_data()

    particle_indicies = ce["particle_index"]
    particle_masses = ce["ParticleMass"]
    prim_mass = np.max(particle_masses)
    prim_index = cef.primary_index(ce)  # Finds the index of the primary core
    prim_coords = cef.primary_coords(ce, lu)  # Finds the corrds of the primary core

    comp_coords = {}
    comp_vels = {}
    comp_seps = {}

    data = str(current_time) + " "

    for i in range(particle_number):
        if i != prim_index:
            coords = [ce["particle_position_x"][i] * lu,
                      ce["particle_position_y"][i] * lu,
                      ce["particle_position_z"][i] * lu]
       
            seps = [coords[0] - prim_coords[0],
                    coords[1] - prim_coords[1],
                    coords[2] - prim_coords[2]]
        
            vels = [ce["particle_velocity_x"][i] * lu,
                    ce["particle_velocity_y"][i] * lu,
                    ce["particle_velocity_z"][i] * lu]

            comp_coords[particle_indicies[i]] = coords
            comp_seps[particle_indicies[i]] = seps
            comp_vels[particle_indicies[i]] = vels
           
            seperation = np.sqrt(seps[0]**2 + seps[1]**2 + seps[2]**2)


            kep_velocity = kepler_velocity(seperation, particle_masses[particle_indicies[i]])
            
            comp_velocity = np.sqrt(comp_vels[particle_indicies[i]][0]**2 + 
                                    comp_vels[particle_indicies[i]][1]**2 + 
                                    comp_vels[particle_indicies[i]][2]**2)

            print(kep_velocity)
            print(comp_velocity)

            data += str(comp_velocity) + " " + str(kep_velocity)
    output_file.write(data)
    output_file.write("\n")  




if __name__ == "__main__":
    # Read the inlist file and return the values of the variables. 
    if len(sys.argv) >= 2:  # Check for supplied inlist
        inlist_path = sys.argv[1]
        (root_dir, exclude_dir, plot_dir, initial_path,
        final_path_plus_one, output_file_name, particle_number) = read_inlist(inlist_path)
    else: 
        # Bring up error if there is no inlist
        print("Inlist File Not Supplied!!!")
        sys.exit(0)

    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    output_file_name = plot_dir + output_file_name

    output_file = open_file(output_file_name, particle_number)

    for index in range(initial_path, final_path_plus_one):
        velocities(root_dir_list[index], index, output_file, particle_number)

    output_file.close()






