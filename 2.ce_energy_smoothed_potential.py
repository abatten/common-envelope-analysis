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

mylog.disabled = True

################################################################################
################################################################################
################################################################################

def read_inlist(ipath):
    """
    Reads the contents of the inlist file and returns the values of the 
    appropriate variables. Returns the values in the following order:
    
    root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one, 
    output_file_name, output_file_append, use_smoothed_potential, 
    smoothing_length
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
    initial_path = config.getint('Common Section', 'initial_path')
    final_path_plus_one = config.getint('Common Section', 'final_path_plus_one')
    output_file_name = config.get('Energy Section', 'output_file_name')
    output_file_append = config.getboolean('Energy Section', 'output_file_append')
    use_smoothed_potential = config.getboolean('Energy Section', 'use_smoothed_potential')
    smoothing_length = config.getfloat('Energy Section', 'smoothing_length')
    select_primary = config.getboolean('Energy Section', 'select_primary')

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
    print("SELECT PRIMARY: " + str(select_primary))
    return (root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one,
           output_file_name, output_file_append, use_smoothed_potential, 
           smoothing_length, select_primary)

def distance(point_1, point_2, units=1):
	"""
	Calculates the distance between two 3D vectors using pythagoras.
	
	point_1: A list of 3 floats or ints.
	point_2: A list of 3 floats or ints.
	"""
	
	distance = (((point_2[0]-point_1[0])*units)**2.0 
			+ ((point_2[1]-point_1[1])*units)**2.0 
			+ ((point_2[2]-point_1[2])*units)**2.0)**0.5
	return distance

def grav_pot(Pmass, Cmass, rad, smoothed, smoothing_length=3.0, smallest_cell=0):
    """
    Calculates the gravitational potential between two objects.
    The potential can be smoothed or not. For the smoothed potential see M. Ruffert 1993
    """
    if smoothed:
        grav_constant = 6.67e-8
        top = grav_constant * Pmass * Cmass
        slsc = smoothing_length * smallest_cell
        bottom = (rad**2.0 + (slsc**2.0 * np.exp((-rad**2.0)/(slsc**2.0))))**0.5
        potential = - top / bottom

    elif not smoothed:
        grav_constant = 6.67e-8
        top = grav_constant * Pmass * Cmass
        bottom = rad 
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
            radius_part = distance(data_coords, particle_coords, length_unit1)
		
	    # Smoothed Gravitational Potential
	    # See M. Ruffert 1993
	    current_Epot_particle[i] = grav_pot(box['ParticleMass'][i], data['CellMass'], radius_part, use_smoothed_potential, smoothing_length, smallest_cell_length)	
        
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

def del_used_variables(delete_variables):
    if delete_variables == True:
        del pf, length_unit1, time_unit1, mass_unit1, common_envelope, current_cycle, current_time, primary_boolean, vacuum_boolean,\
            current_Etherm_primary, current_Etherm_vacuum, current_Etherm_gas, current_Epot_gas, current_Ekin_gas, current_Ekin_part,\
            radius_part_0, current_Epot_part_1, radius_part_to_part, current_Epot_part_to_part, cell_length, smoothing_length,\
            current_Epot_part_to_gas, current_Epot_part, current_Etherm_tot, current_Epot_tot, current_Ekin_tot, current_Etot

def open_file(file_name, append):
    if (append == False): # Overwrite
        output_file = open(file_name, 'w' )

    # Write the first line of information in the file
        output_file.write("Time (yr), Total, Total Kinetic, Total Potential, Total Thermal, " 
                          "Gas Kinetic, Gas potential, Gas Thermal, Particle Kinetic, "
                          "Part-part Potential, Part-gas Potential, Top-grid Cycle (#), "
                          "Thermal Primary, Thermal Vacuum (energies in erg)" + "\n")


    elif (append == True): # Append
        output_file = open(file_name, 'a' )

    return output_file

def ce_energies(directory, index, file):
    pf = load(directory)
    str_index = index2str(index)

    print(" ")
    print("<-------------->")
    print("READING ROOT DIRECTORY " + str_index + ": " + directory)
    print("<-------------->")

    #gets the length, time and mass units used in the current simulation    
    length_unit1 = pf.parameters['LengthUnits']
    time_unit1 = pf.parameters['TimeUnits']
    mass_unit1 = pf.parameters['MassUnits']

    # Adds the whole data as an object called common_envelope. 
    # It is an array that can be access the data by knowing their name through: common_envelope["Dataname"]:
    common_envelope = pf.h.all_data()

    # Current output values
    current_cycle = pf.parameters['InitialCycleNumber']
    current_time = pf.current_time

    #### GAS ####

    # GAS THERMAL ENERGY   
    # The thermal energy is computed separately for primary and vacuum if required
    if select_primary == True:
        primary = common_envelope['SelectPrimary'] # Boolean array of locations of the primary
        vacuum = np.invert(primary) # Boolean array of locations of the vacuum

        current_Etherm_primary = np.sum(common_envelope['ThermalEnergy'][primary] * common_envelope['CellMass'][primary])
        current_Etherm_vacuum = np.sum(common_envelope['ThermalEnergy'][vacuum] * common_envelope['CellMass'][vacuum])

        current_Etherm_gas = current_Etherm_primary + current_Etherm_vacuum

    else: # Primary and vacuum togeter
        current_Etherm_primary = 0.0
        current_Etherm_vacuum = 0.0

        current_Etherm_gas = np.sum(common_envelope['ThermalEnergy'] * common_envelope['CellMass'])
 
    # GAS POTENTIAL ENERGY
    if (pf.parameters['SelfGravity'] == 1):

        current_Epot_gas = np.sum(0.5 * common_envelope['Grav_Potential']  * (length_unit1/time_unit1)**2.0 * common_envelope['CellMass'])
    else:
        current_Epot_gas = 0.0

    # GAS KINETIC ENERGY
    current_Ekin_gas = np.sum(common_envelope['KineticEnergy'] * common_envelope['CellVolume'])

    #### PARTICLES #### 
  
    if pf.parameters['ProblemType'] == 41: # Only works for CE problem type. (i.e.41)
        # KINETIC ENERGY
        current_Ekin_particles = np.sum(0.5 * common_envelope['ParticleMass'] * (common_envelope['ParticleVelocityMagnitude'])**2.0)

        # POTENTIAL ENERGY
        smallest_cell_length =  pf.h.get_smallest_dx() * length_unit1
        # Loop over all the particle and find their potential energy and add them together
        current_Epot_particles_to_gas = 0
        current_Epot_particles = dict()
        for i in range(len(common_envelope['particle_index'])):

            data_coords = (common_envelope['x'], common_envelope['y'], common_envelope['z'])

            particle_coords = (common_envelope['particle_position_x'][i], 
                               common_envelope['particle_position_y'][i], 
                               common_envelope['particle_position_z'][i])

            # Calculate the distance the gas is from the particle.
            radius_particle = distance(data_coords, particle_coords, length_unit1)
            # Calculate Gravitational Potential
            current_Epot_particles[i] = grav_pot(common_envelope['ParticleMass'][i], common_envelope['CellMass'], radius_particle, use_smoothed_potential, smoothing_length, smallest_cell_length)
            # Add the current potential to the running total potential.
            current_Epot_particles_to_gas = current_Epot_particles_to_gas + np.sum(current_Epot_particles[i])

        grav_pot_energy = [0]*len(common_envelope['particle_index'])
        current_Epot_particles = 0
        for i in range(len(common_envelope['particle_index'])):
            for j in range(len(common_envelope['particle_index'])):
                if i != j: # No self potential.
                    particle_i_coords = (common_envelope['particle_position_x'][i],
                                         common_envelope['particle_position_y'][i],
                                         common_envelope['particle_position_z'][i])

                    particle_j_coords = (common_envelope['particle_position_x'][j],
                                         common_envelope['particle_position_y'][j],
                                         common_envelope['particle_position_z'][j])

                    dist = distance(particle_i_coords, particle_j_coords, length_unit1)
                    grav_pot_energy[i] += grav_pot(common_envelope['ParticleMass'][i], common_envelope['ParticleMass'][j], dist, use_smoothed_potential, smoothing_length, smallest_cell_length)
            current_Epot_particles += grav_pot_energy[i]

    else: # If not type 41 set everything to 0.
        current_Ekin_part = 0.0
        current_Epot_particles = 0.0
        current_Epot_particles_to_gas = 0.0

    # Calculate totals
    current_Etherm_tot = current_Etherm_gas
    current_Epot_tot = current_Epot_particles + current_Epot_gas
    current_Ekin_tot = current_Ekin_gas + current_Ekin_particles
    current_Etot = np.sum(current_Etherm_tot + current_Epot_tot + current_Ekin_tot)

    # Write energies to file
    file.write(str(current_time) + " " + str(current_Etot) + " " + str(current_Ekin_tot) + " " + str(current_Epot_tot) + "  " \
             + str(current_Etherm_tot) + " " + str(current_Ekin_gas) + " " + str(current_Epot_gas) + " " + str(current_Etherm_gas) + " "\
             + str(current_Ekin_particles) + " " + str(current_Epot_particles) + " " + str(current_Epot_particles_to_gas) + " " \
             + str(current_cycle) + " " + str(current_Etherm_primary) + " " + str(current_Etherm_vacuum) + "\n")


if __name__ == "__main__":
    # Read the inlist file and return the values of the variables.
    if len(sys.argv) >= 2:
        inlist_path = sys.argv[1]
        root_dir, exclude_dir, plot_dir, initial_path, final_path_plus_one, output_file_name, output_file_append, use_smoothed_potential, smoothing_length, select_primary = read_inlist(inlist_path)
    else:
        print("Inlist File not supplied!")
        sys.exit(0)

    # Add the SelectPrimary field to the avaliable YT fields.
    add_field("SelectPrimary", function=SelectPrimary, units=r"Boolean array")

    # Sort the root directory
    root_dir_list = root_sort(root_dir, exclude=exclude_dir)

    # Set output file name and open it to write
    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, output_file_append)

    # Calculate Energies for every directory and write them
    for index in range(initial_path, final_path_plus_one):
        ce_energies(root_dir_list[index], index, output_file)

