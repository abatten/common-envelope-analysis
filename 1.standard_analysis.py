from __future__ import absolute_import, division, print_function, unicode_literals
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np
import math
import os, sys

####################################################################################################################################
####################################################### User Controls ##############################################################
####################################################################################################################################


global_storage = False # If the data from all the directories is to be stored 
global_plots = False # If on the plots of all outputs are read and done in a single figure
single_plots = False # If all the plots single outputs are done

# The path where the data directories are stored. NOTE: root_dir should only contain output directories
#root_dir = "/path/to/directory"
root_dir = "/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj"

# The directories in root_dir to exclude from the analysis.
exclude_dir = ["1part"]

# The path where the plots will be written. Ensure it ends with a /
#plot_dir = "/path/to/plot/directory/"
plot_dir  = "/disks/ceres/makemake/acomp/abatten/jstaff/sim2/testplots/"

initial_path = 0 # Usually start with zero.
final_path_plus_one = 1 # The number of directories to be analysed. Must be plus one to match with python range functions

# The array of plots to be made. Example: [2,6,9] will produce 3 plots, with indicies [2], [6], [9] in plot_function_list
plot_vector = []

'''
# Single Plots
0: Density vs Radius
1: Radial Velocity vs Radius
2: Kinetic Energy vs Radius
3: Thermal Energy vs Radius
4: Total Energy vs Radius
5: Grav Potential vs Radius
6: X-Velocity vs Radius
7: Y-Velocity vs Radius
8: Z-Velocity vs Radius
9: Density Slice along Z-Axis
10: Density Projection along Z-Axis
11: Pressure Slice along Z-Axis
12: Thermal Energy SLice along Z-Axis
13: Density Gradient Modulus Slice along Z-Axis
14: Velocity Modulus Slice along Z-Axis
15: Temperature Slice along Z-Axis
16: Mach Number Slice along Z-Axis
17: Entropy Slice along Z-Axis

# Global Plots
101: Global Density vs Radius
102: Global Radial Velocity vs Radius
'''

###################################################################################################################################
###################################################################################################################################
###################################################################################################################################

############################################### CONSTANTS FOR CONVERSIONS #########################################################

length_unit = 6.955*10**10 # 1 Rsun in cm
length_unit_km = 1.0*10**5 # 1 km in cm


# system_radius = 3.0e13 # 2 AU in cm
system_radius = 6.0e13 # 4 AU in cm
# system_radius = 1.2e14 # 8 AU in cm


#################################################### TYPES OF PLOTS  #############################################################

plot_type_list = ["Density vs Radius", "Radial Velocity vs Radius", "Density Slice Along Z"]


def density_vs_radius(profile, radial_profile, str_index):
        print("Density vs Radius")
	plt.figure(1) #working on display 1
        plt.clf()
#       plt.ylim(ymin = 1.0e-12, ymax = 1.0e-3) #sets the limits of the y axis
#       plt.xlim(xmin = 1.0e0, xmax = 1.0e3) #sets the limits of the x axis
        plt.xscale("log")
        plt.yscale("log")
        plt.title("Density vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
	plt.ylabel("Density ($g \cdot cm^{-3}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["Density"], '-b')
        plt.savefig(plot_dir + "Density_Radius_" + str_index + ".png")
	print("Created Density vs Radius plot at: " + plot_dir + "Density_Radius_" + str_index + ".png")

def radial_velocity_vs_radius(profile, radial_profile, str_index):
	print("Radial Velocity vs Radius")
        plt.figure(2) #working on display 2
        plt.clf()
        plt.subplots_adjust(left=0.15) #makes room for the labes on the y axis
        #plt.ylim(ymin=-0.0001,ymax=0.0002) #sets the limits of the y axis
        #plt.ylim(ymin=-0.0005,ymax=0.0025) #sets the limits of the y axis
        plt.ylim(ymin=-2.0, ymax=7.0) #sets the limits of the y axis
        plt.xlim(xmin=1.0e0, xmax=1.0e3) #sets the limits of the x axis
        plt.xscale("log")
       #plt.yscale("log")
        plt.title("Radial velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Radial velocity ($km \cdot s^{-1}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["RadialVelocity"]/length_unit_km, '-b') 
        plt.savefig(plot_dir + "RadialVelocity_Radius_" + str_index + ".png")
        print("Created Radial Velocity vs Radius plot at: " + plot_dir + "RadialVelocity_Radius_" + str_index + ".png")

def kinetic_energy_vs_radius(profile, radial_profile, str_index):
        plt.figure(3) #working on display 3
        plt.clf()
       #plt.ylim(ymin=-0.001,ymax=0.002) #sets the limits of the y axis
       #plt.xlim(xmin=1.0e0,xmax=1.0e3) #sets the limits of the x axis
        plt.xscale("log")
       #plt.yscale("log")
        plt.title("Kinetic energy vs Radius (t="+str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Kinetic energy ($erg$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["KineticEnergy"], '-b') 
        plt.savefig(plot_dir + "KineticEnergy_Radius_" + str_index + ".png")
        print("Created Kinetic Energy vs Radius plot at: " + plot_dir + "KineticEnergy_Radius_" + str_index + ".png")

def thermal_energy_vs_radius():
        plt.figure(4) #working on display 4
        plt.clf()
       #plt.ylim(ymin=-0.001,ymax=0.002) #sets the limits of the y axis
       #plt.xlim(xmin=1.0e0,xmax=1.0e3) #sets the limits of the x axis
        plt.xscale("log")
       #plt.yscale("log")
        plt.title("Thermal energy vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Thermal energy ($erg$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["ThermalEnergy"], '-b') 
        plt.savefig(plot_dir + "ThermalEnergy_Radius_" + str_index + ".png")































def density_slice_along_z(profile, radial_profile, str_index):
	print ("Plotting Density Slice Along Z")
	sp = SlicePlot(profile, 'z', "Density", width = 1.0) #centered in the center of the axis
        #sp = SlicePlot(pf, 'y', "Density",'m', width = 1.0) #centered at the max density point
        #sp.set_cmap("Density","Hue Sat Value 2") #sets the color map for the desired field
#       sp.annotate_velocity(factor=8, normalize=False) #overplots velocity field vectors  
        sp.annotate_grids() #overplots the enzo grid
       #sp.annotate_contour('Density', ncont=1, clim=(1.e-8,1.e-8),plot_args={"colors": "red"}) 
       #sp.annotate_contour('Density', ncont=1, clim=(1.0e-10,1.0e-10),plot_args={"colors": "red"})
       #sp.annotate_contour('Density', ncont=1, clim=(6.9314488311118865e-12,6.9314488311118865e-12),plot_args={"colors": "red"}) 
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
        sp.annotate_text((0.7,1.05), "time = "+str(profile.current_time)+"yr") #annotates the current simulation time on the plot
       #sp.annotate_sphere((common_envelope['particle_position_x'][1],common_envelope['particle_position_y'][1]),0.09621083333333333)
#       sp.annotate_image_line((0.0, 0.25), (1.0, 0.25), plot_args={'linewidth':5,'color':'blue'})
#       sp.annotate_image_line((0.0, 0.375), (1.0, 0.375), plot_args={'linewidth':5,'color':'red'})
#       sp.annotate_image_line((0.0, 0.5), (1.0, 0.5), plot_args={'linewidth':5,'color':'yellow'})
#       sp.annotate_image_line((0.0, 0.625), (1.0, 0.625), plot_args={'linewidth':5,'color':'grey'})
#       sp.annotate_image_line((0.0, 0.75), (1.0, 0.75), plot_args={'linewidth':5,'color':'black'})
        sp.save(plot_dir + "density_slice_z_" + str_index + ".png")
        sp.save(plot_dir + "density_slice_z_" + str_index + ".ps")



# List of the plot function names so they can be called in a loop.
plot_function_list = [density_vs_radius, radial_velocity_vs_radius, kinetic_energy_vs_radius, density_slice_along_z]

##################################################################################################################################

def calc_str_index(input_index):
        '''
        calc_str_index converts an integer to a 4 character string to represent that integer.
        This is used so that the output files will be sorted in the correct order.
        i.e. 0000, 0001, 0002 etc rather than 1,10,11,12, etc

        Example:
        calc_str_index(4) returns "0004"
        calc_str_index(105) returns "0105"
        '''

        str_index = str(input_index)
        num_zeros = 4 - len(str_index) # How many zeroes to add.
        new_str_index = [] # Stores the digits of the new index
        for i in range(num_zeros):
                new_str_index.append('0')
        new_str_index.append(str_index)
        str_index = ''.join(new_str_index)
        return(str_index)

def read_root():
        global root_dir
        if root_dir == "0":
                print("No supplied Root Directory")
		sys.exit(0)
        else:
             	if len(sys.argv) >= 2: # If the user specifies an aditional argument
			global plot_vector
			plot_vector = map(int,sys.argv[1].split(",")) # Map that argument to a vector to specify which plots will be created.
	                print("<-------------->")
			print("PLOTS TO BE CREATED ")
	                print("<-------------->")
			for i in range (len(plot_vector)):
				try:
					print(str(plot_vector[i]) + " : " + plot_type_list[plot_vector[i]])
				except:	
					print("ERROR:")
					print(str(plot_vector[i]) + " IS NOT IN THE LIST OF AVALIABLE PLOTS.")
					sys.exit(0)
		print("<-------------->")
                print("READ ROOT DIRECTORY " + " : " + root_dir)
                print("<-------------->")

def read_dir(input_directory, index):
	str_index = calc_str_index(index)

	print("<-------------->")
	print("READ DIRECTORY " + str_index + " : " + input_directory)
	print("<-------------->")

	pf = load(input_directory)

	common_envelope = pf.h.all_data()
	
	# Create binned profiles to evaluate velocity, density, pressure, temperature as a function of radius
	radial_profile = BinnedProfile1D(common_envelope, 200, 'Radius', 0.0, system_radius, log_space = False)

        radial_profile.add_fields("Density")
        radial_profile.add_fields("x-velocity")
	radial_profile.add_fields("y-velocity")
	radial_profile.add_fields("z-velocity")
	radial_profile.add_fields("RadialVelocity")
	radial_profile.add_fields("KineticEnergy")
	radial_profile.add_fields("ThermalEnergy")
	radial_profile.add_fields("TotalEnergy")
	radial_profile.add_fields("Grav_Potential")

	# Store the data if needed.
	if global_storage:
		radial_profile_list.append(radial_profile)

	for i in range(len(plot_vector)):
		plot_function_list[plot_vector[i]](pf, radial_profile,  str_index)		

###################################################### SORT ROOT DIRECTORY FILES ##########################################################

root_dir_list = [] #empty list where to append the values of the root directories to read        
def sort_root():
	print("[SORTING ROOT DIRECTORY FILES]")

	#cycles subfolders and files in the root folder
	for root, dirs, files in os.walk(root_dir):
        	# The os.walk includes the root directory to the list of browsed folders. We don't want this!!
		for i in range(len(exclude_dir)):
	        	if  exclude_dir[i] not in root.split("/"):
				current_folder = root
                		if (current_folder != root_dir):
                        	# cycles subfolders and files in the current sub-folder
	                       		for sub_root, sub_dirs, sub_files in os.walk(root):
	                               		# sorts the files in the subfolder to have the file to pass to yt in posi$
                                		sub_files.sort()
	                               		# Appends the current browsed subdfolder + enzo target file to the list o
                                		root_dir_list.append(os.path.join(root,sub_files[0]))

	# sorts list by directory name
	return root_dir_list.sort()


############################################################################################################################################


if __name__ == "__main__":
	if 'root_dir' not in globals(): # If the user did not define a root directory, ask them for the location
        	root_dir = raw_input("Path to root directory: \n")

        if 'plot_dir' not in globals(): # If the user did not define a plot directory, ask them for the location
                plot_dir = raw_input("Path to plot directory: \n")
	
	sort_root()	
	read_root()

	# Read all the directories in root and produce the plots.
	for index in range(initial_path, final_path_plus_one):
        	read_dir(root_dir_list[index], index)

