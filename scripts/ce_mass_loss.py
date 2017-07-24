from yt.mods import *
import numpy as np
import math
import os
import time as timeclass

import modules.cefunctions as cef
# Time Check
initial_time = timeclass.time()

# Disable YT Outputs on screen (True = Yes, False = No)
yt_outputs = "True"

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
    output_file_name = config.get("Mass Loss Section", "output_file_name")
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
