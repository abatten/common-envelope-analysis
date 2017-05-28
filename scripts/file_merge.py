from __future__ import absolute_import, division, print_function, unicode_literals

import numpy as np
import os
import sys
import ConfigParser
import glob

#directory = str(sys.argv[1])
#file_name = str(sys.argv[2])

directory = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/'
file_name = 'energy_componets'

file_list = glob.glob(directory + "/" + file_name + "*")
print(len(file_list))
file_numbers = []

for i in range(len(file_list)):
    file_numbers.append(file_list[i].split('_'))  # 
    file_numbers[i][2] = file_numbers[i][2][:-4]  # Remove .txt from the files


dict = {}
for i in range(len(file_list)):
    dict[int(file_numbers[i][2])] = file_list[i]

for i in range(len(file_list)):
    print(dict[i])


#sorted_file_numbers = file_numbers[:][2].sort()
#print(sorted_file_numbers)
