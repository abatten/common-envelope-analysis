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
    for i in range(num_particles - 1):
        sepx_dict[str(i+1)] = "Part_%s_rx" % (i+1) + "(cm)"
        sepy_dict[str(i+1)] = "Part_%s_ry" % (i+1) + "(cm)"
        sepz_dict[str(i+1)] = "Part_%s_rz" % (i+1) + "(cm)"
        vx_dict[str(i+1)] = "Part_%s_vx" % (i+1) + "(cm/s)"
        vy_dict[str(i+1)] = "Part_%s_vy" % (i+1) + "(cm/s)"
        vz_dict[str(i+1)] = "Part_%s_vz" % (i+1) + "(cm/s)"
        energy_dict[str(i+1)] = "Part_%s_E" % (i+1) + "(dynes)"
        semi_major_dict[str(i+1)] = "Part_%s_a" % (i+1) + "(cm)"
        eccen_dict[str(i+1)] = "Part_%s_e" % (i+1)

        header = " ".join([header, 
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


def energy_calc(velocity, seperation, prim_mass):
    v_mag = np.linalg.norm(velocity)
    s_mag = np.linalg.norm(seperation)
    GM = CONST.YG * prim_mass * CONST.YMSUN

    #print("GM:", GM)
    #print("v_mag:", v_mag)
    #print("s_mag:", s_mag)
    #print("First Term:",  (v_mag**2 / 2.0))
    #print("Second Term:", (GM / s_mag))
    #print(" ")
    energy = (v_mag**2 / 2.0) - (GM / s_mag)

    return energy

def semimajor_calc(energy, prim_mass):
    GM = CONST.YG * prim_mass * CONST.YMSUN
    semimajor = - GM / (2 * energy)

    return semimajor

def eccen_calc(velocity, seperation, prim_mass):
    GM = CONST.YG * prim_mass * CONST.YMSUN
    v_mag = np.linalg.norm(velocity)
    s_mag = np.linalg.norm(seperation)

    term1 = (v_mag**2 * seperation) / GM
    term2 = ((np.dot(seperation, velocity)) * velocity) / GM
    term3 = seperation / s_mag

    eccen = term1 - term2 - term3

    return eccen


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
    current_time = pf.current_time * tu / year

    #  Adds the whole data as an object called ce. 
    #  It is an array and you can acces the data by knowing their
    #  name through: ce["Dataname"]
    ce = pf.h.all_data()
    
    particle_masses = ce["ParticleMassMsun"]

    number_of_particles = len(particle_masses)

    pdex = ce["particle_index"]

    prim_mass = np.max(particle_masses)
    prim_index = cef.primary_index(ce)
    prim_coords = cef.primary_coords(ce ,lu)


    part_pos_dict = {}
    part_vel_dict = {}
    semimajor_dict = {}
    energy_dict = {}
    eccen_dict = {}
    for i in range(number_of_particles):
        if i != prim_index:
            comp_coords = np.array([ce["particle_position_x"][i] * lu,
                                    ce["particle_position_y"][i] * lu,
                                    ce["particle_position_z"][i] * lu])
           
            comp_vel = np.array([ce["particle_velocity_x"][i],
                                 ce["particle_velocity_y"][i],
                                 ce["particle_velocity_z"][i]])
          
            seperation = cef.distance(prim_coords, comp_coords)
            part_pos_dict[pdex[i]] = comp_coords
            part_vel_dict[pdex[i]] = comp_vel

            #print("Seperation:  ", seperation)
            #print("comp vel:    ", comp_vel)

            energy = energy_calc(comp_vel, seperation, prim_mass)
            semimajor = semimajor_calc(energy, prim_mass)
            eccen = eccen_calc(comp_vel, seperation, prim_mass)
            eccen_mag = np.linalg.norm(eccen)

            energy_dict[pdex[i]] = energy
            semimajor_dict[pdex[i]] = semimajor
            eccen_dict[pdex[i]] = eccen_mag


    data = str(current_time) + " " + str(current_cycle)
    for i in range(len(energy_dict.items())):
        data = " ".join([data, 
                         str(semimajor_dict.items()[i][1]),
                         str(energy_dict.items()[i][1]),
                         str(eccen_dict.items()[i][1]),
                         str(part_pos_dict.items()[i][1][0]),
                         str(part_pos_dict.items()[i][1][1]),
                         str(part_pos_dict.items()[i][1][2]),
                         str(part_vel_dict.items()[i][1][0]),
                         str(part_vel_dict.items()[i][1][1]),
                         str(part_vel_dict.items()[i][1][2])])

    output_file.write(data)
    output_file.write("\n")

                         
            


            #print(eccen_mag)


            #print("Energy:", energy)
            #print("semimajor:", semimajor)
            #print("eccen:", eccen)


    


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
    output_file = open_file(output_file_name, particle_number)

    for index in range(initial_path, final_path_plus_one):
        orbits(root_dir_list[index], index, output_file)
