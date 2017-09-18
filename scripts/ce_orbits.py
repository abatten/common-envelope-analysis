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

    # Read in the config file
    root_dir = config.get("Common", "root_dir")
    exclude_dir = config.get("Common", "exclude_dir")
    plot_dir = config.get("Common", "plot_dir")
    initial_path = config.getint("Common", "initial_path")
    final_path_plus_one = config.getint("Common", "final_path_plus_one")
    output_file_name = config.get("Orbits", "output_file_name")
    particle_number = config.getint("Common", "particle_number")

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))   
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir)) 
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("OUTPUT FILE NAME: " + str(output_file_name))

    return (root_dir, exclude_dir, plot_dir, initial_path, 
            final_path_plus_one, output_file_name, particle_number)

def open_file(file_name, num_particles):
    """
    Creates file and writes the header information

    """
    output_file = open(file_name, "w")

    header = "Time (yr), Cycle(#)"
    
    sepx_dict = {}
    sepy_dict = {}
    sepz_dict = {}
    vx_dict = {}
    vy_dict = {}
    vz_dict = {}
    energy_dict = {}
    semi_major_dict = {}
    eccen_dict = {}
    for i in range(particle_number - 1):
        sepx_dict[str(i+1)] = "Part_%s_rx" % (i+1) + " (cm)"
        sepy_dict[str(i+1)] = "Part_%s_ry" % (i+1) + " (cm)"
        sepz_dict[str(i+1)] = "Part_%s_rz" % (i+1) + " (cm)"
        vx_dict[str(i+1)] = "Part_%s_vx" % (i+1) + " (cm/s)"
        vy_dict[str(i+1)] = "Part_%s_vy" % (i+1) + " (cm/s)"
        vz_dict[str(i+1)] = "Part_%s_vz" % (i+1) + " (cm/s)"
        energy_dict[str(i+1)] = "Part_%s_E" % (i+1) + " (dynes)"
        semi_major_dict[str(i+1)] = "Part_%s_a" % (i+1) + " (cm)"
        eccen_dict[str(i+1)] = "Part_%s_e" % (i+1)

        header = ", ".join([header, 
                            semi_major_dict[str(i+1)],
                            energy_dict[str(i+1)],
                            eccen_dict[str(i+1)],
                            sepx_dict[str(i+1)],
                            sepy_dict[str(i+1)],
                            sepz_dict[str(i+1)],
                            vx_dict[str(i+1)],
                            vy_dict[str(i+1)],
                            vz_dict[str(i+1)]])

    # Write the header of the file
    output_file.write(header + "\n")
    return output_file


def orbits(directory, index, output_file):
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->")
    print("READ ROOTDIRECTORY " + str_index + ":", directory)
    print("<-------------->")

    pf = yt.load(directory)

    #  Gets the length, time and mass units used in the current simulation
    lu = pf.parameters["LengthUnits"]
    tu = pf.parameters["TimeUnits"]
    mu = pf.parameters["MassUnits"]
    current_cycle = pf.parameters["InitialCycleNumber"]
    year = 365.25 * 24 * 60 * 60
    current_time = pf.current_time / year

    #  Adds the whole data as an object called ce. 
    #  It is an array and you can acces the data by knowing their
    #  name through: ce["Dataname"]
    ce = pf.h.all_data()
    
    particle_masses = ce["ParticleMassMsun"]

    number_of_particles = len(particle_masses)

    prim_mass = np.max(particle_masses)
    for i in range(number_of_particles):
        if particle_masses[i] == prim_mass:
            prim_index = i
            break
        else:
            pass

    prim_coords = np.array([ce["particle_position_x"][prim_index] * lu,
                            ce["particle_position_y"][prim_index] * lu,
                            ce["particle_position_z"][prim_index] * lu])


if __name__ == "__main__":
    yt.mylog.disabled = True

    # Check the user supplied the inlist
    if len(sys.argv) >= 2:
        inlist_path = sys.argv[1]
    else:
        print("INLIST NOT SUPPLIED!!!")
        sys.exit(0)

    (root_dir, exclude_dir,
    plot_dir, initial_path,
    final_path_plus_one, output_file_name,
    particle_number) = read_inlist(inlist_path)

    # Sort the root directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    #Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, 2)

    for index in range(initial_path, final_path_plus_one):
        orbits(root_dir_list[index], index, output_file)
