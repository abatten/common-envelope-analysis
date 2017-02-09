"""
This script is designed to analyse the data of a pre-common
envelope simulation. WIth a single star that is inserted in enzo
from a 1D stellar evolution code output, via the YT package.

The script can generate a number of plots as listed below. 
To run this script you can either specify the plots directly 
by adding them to the plot_vector (seperated by commas) or 
by adding them as arguments after you 

--> python <file_name> 0,9,10

Set root_dir to the directory where the data files are stored.
Set plot_sir to where you want the plots to be saved.
Set final_path_plus_one to the number of directories + 1.

       Radial Plots                       Slice Plots
------------------------------------------------------------------------
0: Density vs Radius         |  9: Density Slice along Z-Axis   
1: Radial Velocity vs Radius | 10: Density Projection along Z-Axis
2: Kinetic Energy vs Radius  | 11: Pressure Slice along Z-Axis
3: Thermal Energy vs Radius  | 12: Thermal Energy SLice along Z-Axis
4: Total Energy vs Radius    | 13: Density Gradient Modulus Slice along Z-Axis
5: Grav Potential vs Radius  | 14: Velocity Modulus Slice along Z-Axis
6: X-Velocity vs Radius      | 15: Temperature Slice along Z-Axis
7: Y-Velocity vs Radius      | 16: Mach Number Slice along Z-Axis
8: Z-Velocity vs Radius      | 17: Entropy Slice along Z-Axis
"""

from __future__ import absolute_import, division, print_function, unicode_literals
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys

################################################################################
################################# User Controls ################################
################################################################################


global_storage = False # If the data from all the directories is to be stored 
global_plots = False # If the plots of all outputs are combined into one figure
single_plots = False # If all the plots single outputs are done

# The path where the data directories are stored. 
# NOTE: root_dir should only contain output directories
root_dir = "/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj/1part/"

# The directories in root_dir to exclude from the analysis.
exclude_dir = []

# The path where the plots will be written. Ensure it ends with a /
plot_dir  = "/disks/ceres/makemake/acomp/abatten/jstaff/sim2_part2/plots/"

# Number to start counding from in path.
initial_path = 0 # Usually start with zero.

# The number of directories to be analysed. 
# Plus one to match with python range functions
final_path_plus_one = 177


# The array of plots to be made. 
# Example: [2,6,9] will produce 3 plots, corresponding to 2,6,9 in plot_function_list
plot_vector = []

# Global Plots # Not made yet
#101: Global Density vs Radius
#102: Global Radial Velocity vs Radius

################################################################################
################################################################################
################################################################################

########################### CONSTANTS FOR CONVERSIONS ##########################

length_unit = 6.955*10**10 # 1 Rsun in cm
length_unit_km = 1.0*10**5 # 1 km in cm


# system_radius = 3.0e13 # 2 AU in cm
system_radius = 6.0e13 # 4 AU in cm
# system_radius = 1.2e14 # 8 AU in cm


################################ TYPES OF PLOTS  ###############################

def density_vs_radius(profile, radial_profile, ce, str_index):
	plt.figure(1)
        plt.clf()
#       plt.ylim(ymin = 1.0e-12, ymax = 1.0e-3)
#       plt.xlim(xmin = 1.0e0, xmax = 1.0e3)
        plt.xscale("log")
        plt.yscale("log")
        plt.title("Density vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
	plt.ylabel("Density ($g \cdot cm^{-3}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["Density"], '-b')
        plt.savefig(plot_dir + "Density_Radius_" + str_index + ".png")
	print("Created Density vs Radius plot at: " + plot_dir + "Density_Radius_" + str_index + ".png")

def radial_velocity_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(2)
        plt.clf()
        plt.subplots_adjust(left=0.15) # makes room for the labels on the y axis
#       plt.ylim(ymin=-0.0001,ymax=0.0002) 
#       plt.ylim(ymin=-0.0005,ymax=0.0025) 
        plt.ylim(ymin=-2.0, ymax=7.0) 
        plt.xlim(xmin=1.0e0, xmax=1.0e3) 
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("Radial velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Radial velocity ($km \cdot s^{-1}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["RadialVelocity"]/length_unit_km, '-b') 
        plt.savefig(plot_dir + "RadialVelocity_Radius_" + str_index + ".png")
        print("Created Radial Velocity vs Radius plot at: " + plot_dir + "RadialVelocity_Radius_" + str_index + ".png")

def kinetic_energy_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(3)
        plt.clf()
#       plt.ylim(ymin=-0.001,ymax=0.002) 
#       plt.xlim(xmin=1.0e0,xmax=1.0e3)
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("Kinetic energy vs Radius (t="+str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Kinetic energy ($erg$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["KineticEnergy"], '-b') 
        plt.savefig(plot_dir + "KineticEnergy_Radius_" + str_index + ".png")
        print("Created Kinetic Energy vs Radius plot at: " + plot_dir + "KineticEnergy_Radius_" + str_index + ".png")

def thermal_energy_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(4) 
        plt.clf()
#       plt.ylim(ymin=-0.001,ymax=0.002)
#       plt.xlim(xmin=1.0e0,xmax=1.0e3) 
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("Thermal energy vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Thermal energy ($erg$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["ThermalEnergy"], '-b') 
        plt.savefig(plot_dir + "ThermalEnergy_Radius_" + str_index + ".png")
        print("Created Thermal Energy vs Radius plot at: " + plot_dir + "ThermalEnergy_Radius_" + str_index + ".png")

def total_energy_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(4)
        plt.clf()
#       plt.ylim(ymin=-0.001,ymax=0.002) 
#       plt.xlim(xmin=1.0e0,xmax=1.0e3) 
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("Total energy vs Radius (t="+str(profile.current_time)+" $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Total energy ($erg$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["TotalEnergy"], '-b') 
        plt.savefig(plot_dir + "TotalEnergy_Radius_" + str_index + ".png")
        print("Created Total Energy vs Radius plot at: " + plot_dir + "TotalEnergy_Radius_" + str_index + ".png")

def grav_potential_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(4)
        plt.clf()
#       plt.ylim(ymin=-0.001,ymax=0.002)
#       plt.xlim(xmin=1.0e0,xmax=1.0e3) 
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("Gravitational Potential vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("Gravitational Potential ($pd$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["Grav_Potential"], '-b') 
        plt.savefig(plot_dir + "GravPotential_Radius_" + str_index + ".png")
        print("Created Grav Potential vs Radius plot at: " + plot_dir + "GravPotential_Radius_" + str_index + ".png")


def x_velocity_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(3)
        plt.clf()
#       plt.ylim(ymin=-0.05,ymax=0.01) 
#       plt.xlim(xmin=1.0e0,xmax=1.0e3)
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("x-velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("x-velocity ($km \cdot s^{-1}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["x-velocity"]/length_unit_km, '-b') 
        plt.savefig(plot_dir + "XVelocity_Radius_t" + str(profile.current_time) + ".png")
        print("Created X-velocity vs Radius plot at: " + plot_dir + "XVelocity_Radius_t" + str(profile.current_time) + ".png")

def y_velocity_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(3)
        plt.clf()
#       plt.ylim(ymin=-0.05,ymax=0.01) 
#       plt.xlim(xmin=1.0e0,xmax=1.0e3) 
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("y-velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("y-velocity ($km \cdot s^{-1}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["y-velocity"]/length_unit_km, '-b')
        plt.savefig(plot_dir + "YVelocity_Radius_t" + str(profile.current_time) + ".png")
        print("Created Y-velocity vs Radius plot at: " + plot_dir + "YVelocity_Radius_t" + str(profile.current_time) + ".png")

def z_velocity_vs_radius(profile, radial_profile, ce, str_index):
        plt.figure(3) 
        plt.clf()
#       plt.ylim(ymin=-0.05,ymax=0.01) 
#       plt.xlim(xmin=1.0e0,xmax=1.0e3)
        plt.xscale("log")
#       plt.yscale("log")
        plt.title("z-velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
        plt.xlabel("Radius ($R_{\odot}$)")
        plt.ylabel("z-velocity ($km \cdot s^{-1}$)")
        plt.plot(radial_profile["Radius"]/length_unit, radial_profile["z-velocity"]/length_unit_km, '-b')
        plt.savefig(plot_dir + "ZVelocity_Radius_t" + str(profile.current_time) + ".png")
        print("Created Z-velocity vs Radius plot at: " + plot_dir + "ZVelocity_Radius_t" + str(profile.current_time) + ".png")

def density_slice_along_z(profile, radial_profile, ce, str_index):
	print ("Plotting Density Slice Along Z")
	sp = SlicePlot(profile, 'z', "Density", width = 1.0) #centered in the center of the axis
#       sp = SlicePlot(pf, 'y', "Density",'m', width = 1.0) #centered at the max density point
#       sp.set_cmap("Density","Hue Sat Value 2") 
#       sp.annotate_velocity(factor=8, normalize=False) #overplots velocity field vectors  
#       sp.annotate_grids() #overplots the enzo grid
#       sp.annotate_contour('Density', ncont=1, clim=(1.e-8,1.e-8),plot_args={"colors": "red"}) 
#       sp.annotate_contour('Density', ncont=1, clim=(1.0e-10,1.0e-10),plot_args={"colors": "red"})
#       sp.annotate_contour('Density', ncont=1, clim=(6.9314488311118865e-12,6.9314488311118865e-12),plot_args={"colors": "red"}) 
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
        sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") #annotates the current simulation time on the plot
#       sp.annotate_sphere((common_envelope['particle_position_x'][1],ce['particle_position_y'][1]),0.09621083333333333)
#       sp.annotate_image_line((0.0, 0.25), (1.0, 0.25), plot_args={'linewidth':5,'color':'blue'})
#       sp.annotate_image_line((0.0, 0.375), (1.0, 0.375), plot_args={'linewidth':5,'color':'red'})
#       sp.annotate_image_line((0.0, 0.5), (1.0, 0.5), plot_args={'linewidth':5,'color':'yellow'})
#       sp.annotate_image_line((0.0, 0.625), (1.0, 0.625), plot_args={'linewidth':5,'color':'grey'})
#       sp.annotate_image_line((0.0, 0.75), (1.0, 0.75), plot_args={'linewidth':5,'color':'black'})
        sp.save(plot_dir + "density_slice_z_" + str_index + ".png")
#       sp.save(plot_dir + "density_slice_z_" + str_index + ".ps")

def density_projection_along_z(profile, radial_profile, ce, str_index):
        sp = ProjectionPlot(profile, 'z', "Density", width = 1.0) #centered in the center of the axis 
#       sp.annotate_grids() #overplots the enzo grid 
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles 
        sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") 
        sp.save(plot_dir + "density_proj_z_" + str_index + ".png")
#       sp.save(plot_dir + "density_proj_z_" + str_index + ".ps")

def pressure_slice_along_z(profile, radial_profile, ce, str_index):
        sp = SlicePlot(profile, 'z', "Pressure", width = 1.0)
#       sp.set_cmap("Density","RdBu") #sets the color map for the desired field
        sp.annotate_velocity(factor=10, normalize=False) #overplots velocity field vectors  
#       sp.annotate_grids() #overplots the enzo grid
#       sp.annotate_contour('Density', ncont=1,clim=(6.9314488311118865e-12,6.9314488311118865e-12),plot_args={"colors": "red"}) 
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
        sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time)+"yr") 
#       sp.annotate_sphere((common_envelope['particle_position_x'][0],common_envelope['particle_position_y'][0]),0.1167)
        sp.save(plot_dir + "pressure_slice_z_" + str_index + ".png")

def thermal_energy_slice_along_z(profile, radial_profile, ce, str_index):
        sp = SlicePlot(profile, 'z', "ThermalEnergy", width = 1.0)
        sp.annotate_velocity(factor=10, normalize=False) #overplots velocity field vectors
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particLes
        sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") 
        sp.save(plot_dir + "thermalenergy_slice_z_" + str_index + ".png")

def density_gradient_mod_slice_along_z(profile, radial_profile, ce, str_index):
        sp = SlicePlot(profile, 'z', "gradDensityMagnitude", width = 1.0) 
#       sp.set_cmap("Density","RdBu") #sets the color map for the desired field
#       sp.set_zlim('Density', np.min(common_envelope['Density']), np.max(common_envelope['Density']))
#       sp.set_log('Density',True)
        sp.annotate_velocity(factor=16, normalize=False) #overplots velocity field vectors  
#       sp.annotate_grids() #overplots the enzo grid
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='white') #overplots the projection on the axis of the particles
        sp.annotate_text((0.7, 1.05),"time = " + str(profile.current_time) + "yr") 
#       sp.annotate_sphere((0.5,0.5),0.1167)
        sp.save(plot_dir + "grad_densitymodulus_slice_z_" + str_index + ".png")

def velocity_mod_slice_along_z(profile, radial_profile, ce, str_index):
        sp = SlicePlot(profile, 'z', "VelocityMagnitude", width = 1.0) 
#       sp.set_cmap("Density","RdBu") #sets the color map for the desired field
#       sp.set_zlim('Density', np.min(ce['Density']), np.max(ce['Density']))
#       sp.set_log('Density',True)
        sp.annotate_velocity(factor=16, normalize=False) #overplots velocity field vectors  
#       sp.annotate_grids() #overplots the enzo grid
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
        sp.annotate_text((0.7,1.05),"time = "+str(profile.current_time)+"yr")  
        sp.annotate_sphere((ce['particle_position_x'][1], ce['particle_position_y'][1]), 0.05795833333333333)
        sp.save(plot_dir + "velocitymodulus_slice_z_" + str_index + ".png")

def temperature_slice_along_z(profile, radial_profile, ce, str_index):
        sp = SlicePlot(profile, 'z', "Temperature", width = 1.0)
#       sp.set_cmap("Density", "RdBu") #sets the color map for the desired field
#       sp.set_zlim('Density', np.min(ce['Density']), np.max(ce['Density']))
#       sp.set_log('Density', True)
#       sp.annotate_velocity(factor=16, normalize=False) #overplots velocity field vectors  
#       sp.annotate_grids() #overplots the enzo grid
        sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
        sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr")  
#       sp.annotate_sphere((0.5,0.5),0.1167)
        sp.save(plot_dir + "temperature_slice_z_" + str_index + ".png") 

def mach_number_slice_along_z(profile, radial_profile, ce, str_index):
         sp = SlicePlot(profile, 'z', "MachNumber", width = 1.0) #centered in the center of the axis
#        sp.annotate_velocity(factor=8, normalize=False) #overplots velocity field vectors  
         sp.annotate_contour('MachNumber', ncont=1, clim=(1.,1.), plot_args={"colors": "black"}) 
         sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
         sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr")
         sp.save(plot_dir + "machnumber_slice_z_" + str_index + ".png")
#        sp.save(plot_dir+"machnumber_slice_z_"+str(index)+".ps") 

def entropy_slice_along_z(profile, radial_profile, ce, str_index):
          sp = SlicePlot(profile, 'z', "Entropy", width = 1.0) #centered in the center of the axis
#         sp.annotate_velocity(factor=8, normalize=False) #overplots velocity field vectors  
          sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
          sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") 
          sp.save(plot_dir + "entropy_slice_z_" + str_index + ".png")
#         sp.save(plot_dir+"entropy_slice_z_"+str(index)+".ps")

# List of the plot function names so they can be called in a loop.
plot_functions = [ ["Density vs Radius", density_vs_radius], 
                   ["Radial Velocity vs Radius", radial_velocity_vs_radius],
                   ["Kinetic Energy vs Radius", kinetic_energy_vs_radius],
                   ["Thermal Energy vs Radius", thermal_energy_vs_radius],
                   ["Total Energy vs Radius", total_energy_vs_radius], 
                   ["Grav Potential vs Radius", grav_potential_vs_radius],
                   ["X-Velocity vs Radius", x_velocity_vs_radius], 
                   ["Y-Velocity vs Radius", y_velocity_vs_radius],
                   ["Z-Velocity vs Radius", z_velocity_vs_radius], 
                   ["Density Slice along Z-axis", density_slice_along_z], 
                   ["Density Projection along Z-axis", density_projection_along_z], 
                   ["Pressure Slice aling Z-axis", pressure_slice_along_z],
                   ["Thermal Energy Slice along Z-axis", thermal_energy_slice_along_z], 
                   ["Density Gradient Modulus Slice along Z-axis", density_gradient_mod_slice_along_z], 
                   ["Velocity Modulus Slice along Z-Axis", velocity_mod_slice_along_z],
                   ["Temperature Slice along Z-Axis", temperature_slice_along_z],
                   ["Mach Number Slice along Z-Axis", mach_number_slice_along_z],
                   ["Entropy Slice along Z-Axis", entropy_slice_along_z] ]


################################################################################

def calc_str_index(input_index):
        '''
        calc_str_index converts an integer to a 4 character string to
	 represent that integer.

        This is used so that the output files will be sorted in the 
	correct order.
	
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
			# Map that argument to the plot_vector.
			plot_vector = map(int, sys.argv[1].split(","))
	                print("<-------------->")
			print("PLOTS TO BE CREATED ")
			for i in range (len(plot_vector)):
				try:
					print(str(plot_vector[i]) + " : " + plot_functions[plot_vector[i]][0])
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
#       radial_profile.add_fields("x-velocity")
#	radial_profile.add_fields("y-velocity")
#	radial_profile.add_fields("z-velocity")
#	radial_profile.add_fields("RadialVelocity")
#	radial_profile.add_fields("KineticEnergy")
#	radial_profile.add_fields("ThermalEnergy")
#	radial_profile.add_fields("TotalEnergy")
#	radial_profile.add_fields("Grav_Potential")

	# Store the data if needed.
	if global_storage:
		radial_profile_list.append(radial_profile)

	for i in range(len(plot_vector)):
		plot_functions[plot_vector[i]][1](pf, radial_profile, common_envelope, str_index)		

################################################################################

#empty list where to append the values of the root directories to read        
root_dir_list = []

def sort_root(root_dir_list):
	"""
	Finds the paths to all the enzo data, sorts them based on 
	directory name and adds them to root_dir_list. 

	"""

	print("[SORTING ROOT DIRECTORY FILES]")

	#cycles subfolders and files in the root folder
	for root, dirs, files in os.walk(root_dir):
        	# Skip the direcories that are listed in exclude_dir
		for i in range(len(exclude_dir)):
	        	if  exclude_dir[i] not in root.split("/"):
				current_folder = root

                		# os.walk includes the root directory.
                		# We don't want the root directory!!
                		if (current_folder != root_dir):
                        	# cycles subfolders and files in the current sub-folder
	                       		
					for sub_root, sub_dirs, sub_files in os.walk(root):
	                               		# sorts the files in the subfolder to have the file to pass to yt in position [0]
                                		sub_files.sort()
	                               		# Appends the path of the enzo target file to root_dir_list 
                                		root_dir_list.append(os.path.join(root,sub_files[0]))

	# sorts list by directory name
	root_dir_list.sort()

################################################################################


if __name__ == "__main__":
	# If a root directory isn't defined, ask user for the location
	if 'root_dir' not in globals():
        	root_dir = raw_input("Path to root directory: \n")
        
	# If a plot directory isn't defined, ask user for the location
        if 'plot_dir' not in globals():
                plot_dir = raw_input("Path to plot directory: \n")
	
	sort_root(root_dir_list)
        read_root()

        # If a the final path length isn't defined, ask user for the number
	if 'final_path_plus_one' not in globals():
		final_path_plus_one = input("Number of data directories in root: \n")

	# Read all the directories in root and produce the plots.
	for index in range(initial_path, final_path_plus_one):
		read_dir(root_dir_list[index], index)

        if global_storage:
                print("<--- DATA STORED IN radial_profile --->")
        else:
                print("<--- DATA NOT STORED --->")

