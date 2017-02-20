from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import ConfigParser

from index2str import index2str
from rootsort import root_sort
################################################################################
################################################################################
################################################################################

def read_inlist(ipath):
    """
    Reads the contents of the inlist file and returns the values of
    the appropriate variables. Returns the values in the following order:
	root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one, 
	output_file_name, output_file_append, use_smoothed_potential, smoothing_length
    """
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split('/')[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, 'r'))

    root_dir = config.get('Common Section', 'root_dir')
    exclude_dir = config.get('Common Section', 'exclude_dir')
    plot_dir = config.get('Common Section', 'plot_dir')
    initial_path = config.get('Common Section', 'initial_path')
    final_path_plus_one = config.get('Common Section', 'final_path_plus_one')
    output_file_name = config.get('Energy Section', 'output_file_name')
    output_file_append = config.get('Energy Section', 'output_file_append')
    use_smoothed_potential = config.get('Energy Section', 'use_smoothed_potential')
    smoothing_length = config.get('Energy Section', 'smoothing_length')

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir))
    print("OUTPUT DIRECTORY: " + str(plot_dir))
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("OUTPUT FILE NAME: " + str(output_file_name))
    print("OUTPUT FILE APPEND: " + str(output_file_append))
    print("USE SMOOTHED POTENTIAL: " + str(use_smoothed_potential))
    print("SMOOTHING LENGTH: " + str(smoothing_length))

    return root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one, output_file_name, output_file_append, use_smoothed_potential, smoothing_length

def calc_radius(point_1, point_2):
	"""
	Calculates the distance between two 3D vectors using pythagoras.
	
	point_1: A list of 3 floats or ints.
	point_2: A list of 3 floats or ints.
	"""
	
	radius = ((point_2[0]-point_1[0])**2.0 
			+ (point_2[1]-point_1[1])**2.0 
			+ (point_2[2]-point_1[2])**2.0)**0.5
	
	return radius


def calc_smoothed_grav_pot(Pmass, Cmass, radius, smoothing_length, smallest_cell):
    grav_constant = 6.67e-8
    top = grav_constant * Pmass * Cmass
    bottom = (radius_part**2.0 + (smoothing_length**2.0 * smallest_cell**2.0 * np.exp((-radius_part**2.0)/(smoothing_length**2.0 * smallest_cell**2.0))))**0.5
    potential = - top / bottom
	
    return potential

	
def SelectPrimary(field, data):
    """
    Defines a new YT field containing the indexes of the primary. 
    Is calculated via a boundess criterion.

    """

    # Gets the length, time and mass units used in the current simulation    
    length_unit1 = data.pf.parameters['LengthUnits']
    time_unit1 = data.pf.parameters['TimeUnits']
    mass_unit1 = data.pf.parameters['MassUnits']

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

    	#whole box needed to find the particles
    	box = data.pf.h.all_data()
		
		# Determine the size of the smallest cell for the smoothing length
        smallest_cell_length =  data.pf.h.get_smallest_dx() * length_unit1
		
        radius_particle = dict()
        current_Epot_particle = dict()

        for i in range(len(box['particle_index'])):
            data_coords = (data['x'], data['y'], data['z'])
            particle_coords = (box['particle_position_x'][i], box['particle_position_y'][i], box['particle_position_z'][i])
			
			# Calculate the distance the gas is from the particle.
            radius_part = calc_radius(data_coords, particle_coords)
		
			# Smoothed Gravitational Potential
			# See M. Ruffert 1993
            current_Epot_particle[i] = -(6.67e-8 * box['ParticleMass'][i] * data['CellMass']) / (radius_part**2.0 + (smoothing_length**2.0 * smallest_cell_length**2.0 * np.exp((-radius_part**2.0)/(smoothing_length**2.0 * smallest_cell_length**2.0))))**0.5
		

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

        # Background check
        background = data['Density'] < np.min(density_bound[np.nonzero(density_bound)]) # Boolean array

        primary = np.invert(background) # Boolean array

        return (primary)















if __name__ == "__main__":
	# Read the inlist file and return the values of the variables.
    if len(sys.argv) >= 2:
    	inlist_path = sys.argv[1]
        root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one, output_file_name, output_file_append, use_smoothed_potential, smoothing_length = read_inlist(inlist_path)
    else:
        print("Inlist File not supplied!")
        sys.exit(0)

    # Sort the root directory
    root_dir_list = root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name

    if (output_file_append == 0): # Overwrite
        output_file = open(output_file_name, 'w' )

    # Write the first line of information in the file
        output_file.write("Time (yr), Total, Total kinetic, Total potential, Total thermal, Gas kinetic, Gas potential, Gas thermal, Particle kinetic, Part-part potential, Part-gas potential, Top-grid cycle (#), Thermal primary, Thermal vacuum (energies in erg)"+"\n")

    elif (output_file_append == 1): # Append
        output_file = open(output_file_name, 'a' )
	
    test = calc_smoothed_grav_pot(10, 5, 300, smoothing_length, 20)
    print(test)
