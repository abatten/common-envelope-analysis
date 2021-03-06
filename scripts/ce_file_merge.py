from __future__ import absolute_import, print_function

import glob

# directory = str(sys.argv[1])
# file_name = str(sys.argv[2])
directory = "/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2_mass_added/"
file_name = "ce_mass_loss_three"
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
    file_numbers.append(file_list[i].split('_'))
    print(file_numbers[i][-1], file_list[i])
    file_numbers[i][-1] = file_numbers[i][-1][:-4]  # Remove .txt from files

# Sort the files into the correct order
dict = {}
for i in range(len(file_list)):
    print(file_numbers[i][-1])
    #  Associate the suffix number with the file
    dict[int(file_numbers[i][-1])] = file_list[i]

# Create output file
out_file = open(directory + output_file_name, 'w')
# Write the combined output file
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
print("CREATED: " + directory + output_file_name)
print("<---------->")
print(" ")
