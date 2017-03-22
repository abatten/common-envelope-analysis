from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import ConfigParser

import cefunctions as cef

mylog.disabled = True

######################################################################

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
    final_path_plus_one = config.getint('Common Section', 'final_path_plus_one')
    output_file_name = config.get('Angular Momentum Section', 'output_file_name')
    output_file_append = config.getboolean('Angular Momentum Section', 'output_file_append')
    angular_momentum_wrt_com = config.getboolean("Angular Momentum Section", "angular_momentum_wrt_com")
    smoothing_length = config.getfloat('Angular Momentum Section', 'smoothing_length')

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("OUTPUT FILE NAME: " + str(output_file_name))
    print("OUTPUT FILE APPEND: " + str(output_file_append))
    print("ANGULAR MOMENTUM WRT COM: " + str(angular_momentum_wrt_com))
    print("SMOOTHING LENGTH: " + str(smoothing_length))
    return (root_dir, exclude_dir, plot_dir, initial_path, 
            final_path_plus_one, output_file_name, 
            output_file_append,angular_momentum_wrt_com, 
            smoothing_length)

def BoundMinDensity(field, data):
    # Get the length, time and mass units used in the current simulation  
    length_unit1 = data.pf.parameters["LengthUnits"]
    time_unit1 = data.pf.parameters["TimeUnits"]
    mass_unit1 = data.pf.parameters["MassUnits"]

    # Step 1: Thermal Energy of Gas
    current_Etherm_gas = data['ThermalEnergy'] * data['CellMass']
    # Step 2: Kinetic Energy of the Gas
    current_Ekin_gas = data['KineticEnergy'] * data['CellVolume']
    # Step 3: Potential Energy of the Gas
    # The potential energy is computed only if the potential field is available
    if (data.pf.parameters['SelfGravity'] == 1):
        current_Epot_gas = data['Grav_Potential']  * (length_unit1/time_unit1)**2.0 * data['CellMass']
    else:
        current_Epot_gas = 0.0

    # Step 4: Potential Energy of the particle to the gas
    # Working only for the common envelope problem type
    if (data.pf.parameters['ProblemType']) == 41:
        # Whole box needed to find the particles
        box = data.pf.h.all_data()

        # Determine the size of the smallest cell for the smoothing length
        smallest_cell_length =  data.pf.h.get_smallest_dx() * length_unit1
        radius_particle = dict()
        current_Epot_particle = dict()

        for i in range(len(box['particle_index'])):

            data_coords = (data['x'], data['y'], data['z'])
            particle_coords = (box['particle_position_x'][i],
                               box['particle_position_y'][i],
                               box['particle_position_z'][i])

            # Calculate the distance the gas is from the particle.
            radius_part = cef.distance(data_coords, particle_coords, length_unit1)

            # Smoothed Gravitational Potential
            # See M. Ruffert 1993
            current_Epot_particle[i] = cef.grav_pot(box['ParticleMass'][i], data['CellMass'], radius_part, use_smoothed_potential, smoothing_length, smallest_cell_length)

        for i in range(len(box['particle_index'])):

            data_coords = (data['x'], data['y'], data['z'])
            particle_coords = (box['particle_position_x'][i],
                               box['particle_position_y'][i],
                               box['particle_position_z'][i])

            # Calculate the distance the gas is from the particle.
            radius_part = cef.distance(data_coords, particle_coords, length_unit1)

            # Smoothed Gravitational Potential
            # See M. Ruffert 1993
            current_Epot_particle[i] = cef.grav_pot(box['ParticleMass'][i], data['CellMass'], radius_part, use_smoothed_potential, smoothing_length, smallest_cell_length)

        current_Epot_part_to_gas = 0
        for i in range(len(box['particle_index'])):
            current_Epot_part_to_gas = current_Epot_part_to_gas + current_Epot_particle[i]

    else:
        current_Epot_part_to_gas = 0.0

    # Computes the total energy of the gas in each cell
    # Step 1 + Step 2 + Step 3 + Step 4
    Etot = current_Etherm_gas + current_Epot_gas + current_Ekin_gas + current_Epot_part_to_gas

    # Locations of bound and unbound cells
    bound = Etot < 0.0
    unbound = Etot > 0.0
    # Creates a density array that is 0 where the gas is unbound and not background, the actual density where not
    density_bound = data['Zeros']
    density_bound[np.where(bound)] = data['Density'][np.where(bound)]
    
    # Bound minimum density
    bound_min_density = np.min(density_bound)

def initialise_arrays():
    current_time = np.ndarray([0],dtype=float) #times

    cycle = np.ndarray([0],dtype=float) #top grid cycles

    Lp_x = np.ndarray([0],dtype=float) #particles x angular momentum
    Lp_y= np.ndarray([0],dtype=float) #particles y angular momentum
    Lp_z = np.ndarray([0],dtype=float) #particles z angular momentum

    Lg_x = np.ndarray([0],dtype=float) #gas x angular momentum
    Lg_y = np.ndarray([0],dtype=float) #gas y angular momentum
    Lg_z = np.ndarray([0],dtype=float) #gas z angular momentum

    Ltot_x = np.ndarray([0],dtype=float) #tot x angular momentum
    Ltot_y = np.ndarray([0],dtype=float) #tot y angular momentum
    Ltot_z = np.ndarray([0],dtype=float) #tot z angular momentum

    return current_time, cycle, Lp_x, Lp_y, Lp_z, Lg_x, Lg_y, Lg_z, Ltot_x, Ltot_y, Ltot_z 

def open_file(file_name, append):
    if (append == False): # Overwrite
        output_file = open(file_name, 'w' )

    # Write the first line of information in the file
        output_file.write("Time(yr), Cycle(#), Lp_x(gcm/s), Lp_y(gcm/s), Lp_z(gcm/s), Lg_x(gcm/s),"
                          "Lg_y(gcm/s), Lg_z(gcm/s), Ltot_x(gcm/s), Ltot_y(gcm/s), Ltot_z(gcm/s)"+"\n")

    elif (append == True): # Append
        output_file = open(file_name, 'a' )

    return output_file

def ce_angular_momentum(directory, index, file):
    return None


if __name__ == "__main__":
    # Read the inlist file and return the values of the variables.
    if len(sys.argv) >= 2:
        inlist_path = sys.argv[1]

        (root_dir, exclude_dir, plot_dir, initial_path, 
         final_path_plus_one, output_file_name, 
         output_file_append, angular_momentum_wrt_com, 
         smoothing_length) = read_inlist(inlist_path)
    
    else:
        print("Inlist File not supplied!")
        sys.exit(0)
    
    # Adds the field to the avaliable yt fields
    add_field("BoundMinDensity", function=BoundMinDensity, units=r"\rm{g}/\rm{cm}^{3}")
    
    # Initialise arrays    
    (current_time, cycle, Lp_x, Lp_y, Lp_z, 
     Lg_x, Lg_y, Lg_z, Ltot_x, Ltot_y, Ltot_z) = initialise_arrays()

    # Sort the root directory
    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, output_file_append)
