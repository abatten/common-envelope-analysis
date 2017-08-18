from __future__ import absolute_import, division, print_function, unicode_literals

import yt.mods as yt

import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import ConfigParser

import cemodules.cefunctions as cef

yt.mylog.disabled = True

##############################################################################
def read_inlist(ipath):
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split('/')[-1]
    config =ConfigParser()
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
    header = "Time (yr), Cycle(#), CoMx, CoMy, CoMz"

    output_file.write(header + "\n")

    return output_file


def centre_of_mass(directory, index, outfile):
    


    
    outfile.write(data)
    outfile.write("\n")

    return


