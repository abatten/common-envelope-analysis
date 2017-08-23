from __future__ import absolute_import, division, print_function, unicode_literals

import yt.mods as yt

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import ConfigParser

import cemodules.cefunctions as cef
import cemodules.constants as CONST

yt.mylog.disabled = True

##############################################################################
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
    output_file_name = config.get("COM", "output_file_name")

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("OUTPUT FILE NAME: " + str(output_file_name))

    return (root_dir, exclude_dir, plot_dir, initial_path,
            final_path_plus_one, output_file_name)

def open_file(file_name):
    output_file = open(file_name, "w")
    header = "Time (yr), Cycle(#), CoMx, CoMy, CoMz, TotalMass(Msun)"

    output_file.write(header + "\n")

    return output_file


def centre_of_mass(directory, index, outfile):
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->")
    print("READ ROOTDIRECTOY" + str_index + ":", directory)
    print("<-------------->")

    pf = yt.load(directory)

    lu = pf.parameters["LengthUnits"]
    tu = pf.parameters["TimeUnits"]

    current_time = pf.current_time * tu / CONST.XYR
    
    current_cycle = pf.parameters["InitialCycleNumber"]

    ce = pf.h.all_data()
    com = ce.quantities.functions['CenterOfMass'][1](ce, use_cells=True, use_particles=True)

    data = str(current_time) + " " + str(current_cycle) + " "

    for i in range(len(com)):
        data += str(com[i]) + " "
    
    outfile.write(data)
    outfile.write("\n")

    return

if __name__ == "__main__":
    yt.mylog.disabled = True
    (root_dir, exclude_dir, plot_dir, initial_path,
    final_path_plus_one, output_file_name) = read_inlist(sys.argv[1])

    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name)

    for index in range(initial_path, final_path_plus_one):
        centre_of_mass(root_dir_list[index], index, output_file)


