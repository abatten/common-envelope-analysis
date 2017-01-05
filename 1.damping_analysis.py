#from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np
import math
import os


global_storage = False # If the data from all the directories is to be stored 
global_plots = False # If on the plots of all outputs are read and done in a single figure
single_plots = False # If all the plots single outputs are done


# The path where the data directories are stored. NOTE: root_dir must only contain output directories
root_dir = "/path/to/directory"

# The path where the plots will be written
plot_dir = "/path/to/plot/directory"



initial_path = 0
final_path_plus_one = 1 #Must be plus one to match with python range functions


### Constants for conversions ###

length_unit = 6.955*10**10 # 1Rsun in cm
length_unit_km = 1.0*10**5 # 1km in cm



system_radius = 6.0e13 # 4AU in cm
# system_radius = 3.0e13 # 2AU in cm
# system_radius = 1.2e14 # 8AU in cm


root_dir_list = []
index = 0


def read_root():
	root_dir = input("Path to root directory: \n")
	if root_dir == 0:
		print('Yes')
	else:
		print("|--------|")
		print("READ ROOT DIRECTORY " + str(index) + ":", root_dir)
		print("|--------|")


read_root()