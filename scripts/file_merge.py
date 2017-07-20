from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import os
import sys
import ConfigParser
import glob

directory = str(sys.argv[1])
file_name = str(sys.argv[2])

#directory = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/'
#file_name = 'energy_componets'

output_file_name = 'combined_' + file_name + '_file.txt'

print(" ")
print("<---------->")
print("COMBINING OUTPUT FILES")
print("<---------->")
print(" ")

# List all files in the directory with the given name
file_list = glob.glob(directory + "/" + file_name + "*")
file_numbers = []

# Find the suffix number for each file
for i in range(len(file_list)):
    file_numbers.append(file_list[i].split('_'))  #
    file_numbers[i][2] = file_numbers[i][2][:-4]  # Remove .txt from the files

# Sort the files into the correct order
dict = {}
for i in range(len(file_list)):
    dict[int(file_numbers[i][2])] = file_list[i]  # Associate the suffix number with the file

# Write the combined output file
out_file = open(directory + output_file_name, 'w')
for i in range(len(file_list)):
    print(dict[i])
    f = open(dict[i], 'r')
    lines = f.readlines()
    if i == 0:  # If it is the first file, write the header
        out_file.write(lines[0])
    out_file.write(lines[1])
    

print(" ")
print("<---------->")
print("COMBINE COMPLETE")
print("<---------->")
print(" ")

