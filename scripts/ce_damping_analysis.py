#!/usr/bin/env python
"""
This script is designed to analyse the data of a pre-common
envelope simulation. With a single star that is inserted in enzo
from a 1D stellar evolution code output, via the YT package.

The script can generate a number of plots as listed below.
To run this script you specify the plots by adding them as
arguments when calling the python script.

>>> python <file_name> <path/to/inlist> <>

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
                             | 18: Gravitational Potential Slice aling Z-axis

"""

from __future__ import absolute_import, division, print_function

import yt.mods as yt

import matplotlib.pyplot as plt
import sys
import ConfigParser

import cemodules.cefunctions as cef
import cemodules.constants as CONST
########################## CONSTANTS FOR CONVERSIONS #########################

length_unit = 6.955*10**10  # 1 Rsun in cm
length_unit_km = 1.0*10**5  # 1 km in cm

# system_radius = 3.0e13  # 2 AU in cm
system_radius = 6.0e13  # 4 AU in cm
# system_radius = 1.2e14  # 8 AU in cm

############################### TYPES OF PLOTS  ##############################


def density_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(1)
    plt.clf()
    plt.xscale("log")
    plt.yscale("log")
    plt.title("Density vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("Density ($g \cdot cm^{-3}$)")
    plt.plot(radial_profile["Radius"]/CONST.YRSUN, radial_profile["Density"], '-b')
    plt.savefig(plot_dir + "Density_Radius_" + str_index + ".png")
    print("Created Density vs Radius plot at: " + plot_dir + "Density_Radius_" + str_index + ".png")


def radial_velocity_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(2)
    plt.clf()
    plt.subplots_adjust(left=0.15)  # makes room for the labels on the y axis
    plt.ylim(ymin=-2.0, ymax=7.0)
    plt.xlim(xmin=1.0e0, xmax=1.0e3)
    plt.xscale("log")
    plt.title("Radial velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("Radial velocity ($km \cdot s^{-1}$)")
    plt.plot(radial_profile["Radius"]/CONST.YRSUN, radial_profile["RadialVelocity"]/length_unit_km, '-b')
    plt.savefig(plot_dir + "RadialVelocity_Radius_" + str_index + ".png")
    print("Created Radial Velocity vs Radius plot at: " + plot_dir + "RadialVelocity_Radius_" + str_index + ".png")


def kinetic_energy_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(3)
    plt.clf()
    plt.xscale("log")
    plt.title("Kinetic energy vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("Kinetic energy ($erg$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["KineticEnergy"], '-b')
    plt.savefig(plot_dir + "KineticEnergy_Radius_" + str_index + ".png")
    print("Created Kinetic Energy vs Radius plot at: " + plot_dir + "KineticEnergy_Radius_" + str_index + ".png")


def thermal_energy_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(4)
    plt.clf()
    plt.xscale("log")
    plt.title("Thermal energy vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("Thermal energy ($erg$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["ThermalEnergy"], '-b')
    plt.savefig(plot_dir + "ThermalEnergy_Radius_" + str_index + ".png")
    print("Created Thermal Energy vs Radius plot at: " + plot_dir + "ThermalEnergy_Radius_" + str_index + ".png")


def total_energy_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(4)
    plt.clf()
    plt.xscale("log")
    plt.title("Total energy vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("Total energy ($erg$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["TotalEnergy"], '-b')
    plt.savefig(plot_dir + "TotalEnergy_Radius_" + str_index + ".png")
    print("Created Total Energy vs Radius plot at: " + plot_dir + "TotalEnergy_Radius_" + str_index + ".png")


def grav_potential_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(4)
    plt.clf()
    plt.xscale("log")
    plt.title("Gravitational Potential vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("Gravitational Potential ($pd$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["Grav_Potential"], '-b')
    plt.savefig(plot_dir + "GravPotential_Radius_" + str_index + ".png")
    print("Created Grav Potential vs Radius plot at: " + plot_dir + "GravPotential_Radius_" + str_index + ".png")


def x_velocity_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(3)
    plt.clf()
    plt.xscale("log")
    plt.title("x-velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("x-velocity ($km \cdot s^{-1}$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["x-velocity"]/length_unit_km, '-b')
    plt.savefig(plot_dir + "XVelocity_Radius_t" + str(profile.current_time) + ".png")
    print("Created X-velocity vs Radius plot at: " + plot_dir + "XVelocity_Radius_t" + str(profile.current_time) + ".png")


def y_velocity_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(3)
    plt.clf()
    plt.xscale("log")
    plt.title("y-velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("y-velocity ($km \cdot s^{-1}$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["y-velocity"]/length_unit_km, '-b')
    plt.savefig(plot_dir + "YVelocity_Radius_t" + str(profile.current_time) + ".png")
    print("Created Y-velocity vs Radius plot at: " + plot_dir + "YVelocity_Radius_t" + str(profile.current_time) + ".png")


def z_velocity_vs_radius(profile, radial_profile, ce, str_index):
    plt.figure(3)
    plt.clf()
    plt.xscale("log")
    plt.title("z-velocity vs Radius (t=" + str(profile.current_time) + " $yr$)")
    plt.xlabel("Radius ($R_{\odot}$)")
    plt.ylabel("z-velocity ($km \cdot s^{-1}$)")
    plt.plot(radial_profile["Radius"]/length_unit, radial_profile["z-velocity"]/length_unit_km, '-b')
    plt.savefig(plot_dir + "ZVelocity_Radius_t" + str(profile.current_time) + ".png")
    print("Created Z-velocity vs Radius plot at: " + plot_dir + "ZVelocity_Radius_t" + str(profile.current_time) + ".png")


def density_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "Density", width=1.0) #centered in the center of the axis
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
    sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") #annotates the current simulation time on the plot
#   sp.annotate_image_line((0.5, 0.0), (0.5, 1.0), plot_args={'linewidth':1,'color':'black'})
#   sp.annotate_image_line((0.0, 0.5), (1.0, 0.5), plot_args={'linewidth':1,'color':'black'})
    sp.save(plot_dir + "density_slice_z_" + str_index + ".png")


def density_projection_along_z(profile, radial_profile, ce, str_index):
    sp = yt.ProjectionPlot(profile, 'z', "Density", width = 1.0) #centered in the center of the axis 
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles 
    sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") 
    sp.save(plot_dir + "density_proj_z_" + str_index + ".png")


def pressure_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "Pressure", width = 1.0)
    sp.annotate_velocity(factor=10, normalize=False) #overplots velocity field vectors  
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
    sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time)+"yr") 
    sp.save(plot_dir + "pressure_slice_z_" + str_index + ".png")


def thermal_energy_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "ThermalEnergy", width = 1.0)
    #sp.annotate_velocity(factor=10, normalize=False) #overplots velocity field vectors
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particLes
    sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") 
    sp.save(plot_dir + "thermalenergy_slice_z_" + str_index + ".png")


def density_gradient_mod_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "gradDensityMagnitude", width = 1.0) 
    sp.annotate_velocity(factor=16, normalize=False) #overplots velocity field vectors  
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='white') #overplots the projection on the axis of the particles
    sp.annotate_text((0.7, 1.05),"time = " + str(profile.current_time) + "yr") 
    sp.save(plot_dir + "grad_densitymodulus_slice_z_" + str_index + ".png")


def velocity_mod_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "VelocityMagnitude", width = 1.0) 
    sp.annotate_velocity(factor=16, normalize=False)  # Overplots velocity field vectors  
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black')  # Overplots the projection on the axis of the particles
    sp.annotate_text((0.7, 1.05),"time = "+str(profile.current_time)+"yr")  
    sp.annotate_sphere((ce['particle_position_x'][1], ce['particle_position_y'][1]), 0.05795833333333333)
    sp.save(plot_dir + "velocitymodulus_slice_z_" + str_index + ".png")


def temperature_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "Temperature", width = 1.0)
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black')  # Overplots the projection on the axis of the particles
    sp.annotate_text((0.7, 1.05), "time = " + str(profile.current_time) + "yr")  
    sp.save(plot_dir + "temperature_slice_z_" + str_index + ".png") 


def mach_number_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "MachNumber", width = 1.0)  # Centered in the center of the axis
    sp.annotate_contour('MachNumber', ncont=1, clim=(1., 1.), plot_args={"colors": "black"}) 
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black')  # Overplots the projection on the axis of the particles
    sp.annotate_text((0.7, 1.05), "time = " + str(profile.current_time) + "yr")
    sp.save(plot_dir + "machnumber_slice_z_" + str_index + ".png")


def entropy_slice_along_z(profile, radial_profile, ce, str_index):
    sp = yt.SlicePlot(profile, 'z', "Entropy", width = 1.0) #centered in the center of the axis
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
    sp.annotate_text((0.7,1.05), "time = " + str(profile.current_time) + "yr") 
    sp.save(plot_dir + "entropy_slice_z_" + str_index + ".png")


def grav_potential_slice_along_z(profile, radial_profile, ce, str_index):
    print("Plotting Gravitational Potential Slice Along Z")
    sp = yt.SlicePlot(profile, 'z', "Grav_Potential", width=1.0)  # Centered in the center of the axis
    sp.annotate_grids()  # Overplots the enzo grid
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black')  # Overplots the projection on the axis of the particles
    sp.annotate_text((0.7, 1.05), "time = " + str(profile.current_time) + "yr")  # Annotates the current simulation time on the plot
    sp.save(plot_dir + "grav_pot_slice_z_" + str_index + ".png")


# List of the plot function names so they can be called in a loop.
plot_functions = [["Density vs Radius", density_vs_radius],
                  ["Radial Velocity vs Radius", radial_velocity_vs_radius],
                  ["Kinetic Energy vs Radius", kinetic_energy_vs_radius],
                  ["Thermal Energy vs Radius", thermal_energy_vs_radius],
                  ["Total Energy vs Radius", total_energy_vs_radius],
                  ["Grav Potential vs Radius", grav_potential_vs_radius],
                  ["X-Velocity vs Radius", x_velocity_vs_radius],
                  ["Y-Velocity vs Radius", y_velocity_vs_radius],
                  ["Z-Velocity vs Radius", z_velocity_vs_radius],
                  ["Density Slice Z-axis", density_slice_along_z],
                  ["Density Projection Z-axis", density_projection_along_z],
                  ["Pressure Slice Z-axis", pressure_slice_along_z],
                  ["Thermal Energy Slice Z-axis", thermal_energy_slice_along_z],
                  ["Density Grad Mod Slice Z-axis", density_gradient_mod_slice_along_z],
                  ["Velocity Mod Slice Z-Axis", velocity_mod_slice_along_z],
                  ["Temperature Slice Z-Axis", temperature_slice_along_z],
                  ["Mach Number Slice Z-Axis", mach_number_slice_along_z],
                  ["Entropy Slice Z-Axis", entropy_slice_along_z],
                  ["Grav. Pot. Slice Z-axis", grav_potential_slice_along_z]]


def read_inlist(ipath):
    """
    Reads the contents of the inlist file and returns the values of
    the appropriate variables. Returns the values in the following order:
    root_dir, exclude_dir, plot_dir, initial_path,
    final_path_plus_one, global_storage
    """
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split('/')[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, 'r'))

    global_storage = config.getboolean('Damping Analysis', 'global_storage')
    root_dir = config.get('Common', 'root_dir')
    exclude_dir = config.get('Common', 'exclude_dir').split(',')
    plot_dir = config.get('Common', 'plot_dir')
    initial_path = config.getint('Common', 'initial_path')
    final_path_plus_one = config.getint('Common', 'final_path_plus_one')

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("GLOBAL STORAGE: " + str(global_storage))

    return (root_dir, exclude_dir, plot_dir,
            initial_path, final_path_plus_one,
            global_storage)


def args2plots():
    if len(sys.argv) >= 3:  # If the user specifies an aditional argument
        # Map that argument to the plot_vector.
        plot_vector = map(int, sys.argv[2].split(","))
        print(" ")
        print("<-------------->")
        print("PLOTS TO BE CREATED ")
        print("<-------------->")
        for i in range(len(plot_vector)):
            try:
                print(str(plot_vector[i]) + " : " +
                      plot_functions[plot_vector[i]][0])
            except:
                print("ERROR:")
                print(str(plot_vector[i]) + "NOT ONE OF AVALIABLE PLOTS.")
                sys.exit(0)
    else:
        print(" ")
        print("<-------------->")
        print("NO PLOTS TO BE CREATED")
        print("<-------------->")
        plot_vector = None

    return plot_vector


def read_data(input_directory, index):
    str_index = cef.index2str(index)
    print(" ")
    print("<-------------->")
    print("READ DIRECTORY " + str_index + " : " + input_directory)
    print("<-------------->")

    pf = yt.load(input_directory)
    common_envelope = pf.h.all_data()

    # Create binned profiles to evaluate:
    # velocity, density, pressure and temperature as a function of radius
    radial_profile = yt.BinnedProfile1D(common_envelope, 200, "Radius", 0.0,
                                        system_radius, log_space=False)

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

    # Make the specified plots
    if plot_vector is not None:
        for i in range(len(plot_vector)):
            plot_functions[plot_vector[i]][1](pf, radial_profile,
                                              common_envelope, str_index)


if __name__ == "__main__":

    # Read the inlist file and return the values of the variables.
    if len(sys.argv) >= 2:
        (root_dir, exclude_dir, plot_dir,
         initial_path, final_path_plus_one,
         global_storage) = read_inlist(sys.argv[1])

        root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)
        plot_vector = args2plots()
        # Read all the directories in root and produce the plots.
        for index in range(initial_path, final_path_plus_one):
            read_data(root_dir_list[index], index)

        if global_storage:
            print(" ")
            print("<-------------->")
            print("DATA STORED IN radial_profile")
            print("<-------------->")
            print(" ")
        else:
            print(" ")
            print("<-------------->")
            print("DATA NOT STORED")
            print("<-------------->")
            print(" ")
