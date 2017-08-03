from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import ConfigParser

import cemodules.cefunctions as cef

mylog.disabled = True

######################################################################

def read_inlist(ipath):
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split('/')[-1]
    cfig = ConfigParser.ConfigParser()
    cfig.readfp(open(inlist_name, 'r'))
   
    # Read in the config file
    root_dir = cfig.get('Common Section', 'root_dir')
    exclude_dir = cfig.get('Common Section', 'exclude_dir')
    plot_dir = cfig.get('Common Section', 'plot_dir')
    initial_path = cfig.getint('Common Section', 'initial_path')
    final_path_plus_one = cfig.getint('Common Section', 'final_path_plus_one')
    output_file_name = cfig.get('Angular Momentum Section', 'output_file_name')
    output_file_append = cfig.getboolean('Angular Momentum Section', 'output_file_append')
    angular_momentum_wrt_com = cfig.getboolean("Angular Momentum Section", "angular_momentum_wrt_com")
    smoothing_length = cfig.getfloat('Common Section', 'smoothing_length')
    particle_number = cfig.getint('Common Section', 'particle_number')

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
    print("PARTICLE NUMBER: " + str(particle_number))

    return (root_dir, exclude_dir, plot_dir, initial_path, 
            final_path_plus_one, output_file_name, 
            output_file_append,angular_momentum_wrt_com, 
            smoothing_length, particle_number)

def BoundMinDensity(field, data):
    #  Get the length, time and mass units used in the current simulation  
    lu = data.pf.parameters["LengthUnits"]
    tu = data.pf.parameters["TimeUnits"]
    mu = data.pf.parameters["MassUnits"]

    #  Step 1: Thermal Energy of Gas
    Etherm_gas = data['ThermalEnergy'] * data['CellMass']

    #  Step 2: Kinetic Energy of the Gas
    Ekin_gas = data['KineticEnergy'] * data['CellVolume']

    #  Step 3: Potential Energy of the Gas
    #  Potential energy is computed only if the potential field is available
    if (data.pf.parameters['SelfGravity'] == 1):
        Epot_gas = data['Grav_Potential'] * (lu / tu)**2.0 * data['CellMass']
    else:
        Epot_gas = 0.0

    # Step 4: Potential Energy of the particle to the gas
    # Working only for the common envelope problem type
    if (data.pf.parameters['ProblemType']) == 41:
        # Whole ce needed to find the particles
        ce = data.pf.h.all_data()

        # Determine the size of the smallest cell for the smoothing length
        smallest_cell_length =  data.pf.h.get_smallest_dx() * lu
 
        radius_particle = dict()
        Epot_particle = dict()

        for i in range(len(ce['particle_index'])):
            data_coords = (data['x'], data['y'], data['z'])
            particle_coords = (ce['particle_position_x'][i],
                               ce['particle_position_y'][i],
                               ce['particle_position_z'][i])

            # Calculate the distance the gas is from the particle/.
            radius_part = cef.distance(data_coords, particle_coords, lu)

            # Smoothed Gravitational Potential
            # See M. Ruffert 1993
            Epot_particle[i] = cef.grav_pot(ce['ParticleMass'][i], 
                                                    data['CellMass'], 
                                                    radius_part, False, 
                                                    smoothing_length, 
                                                    smallest_cell_length)

        Epot_part_to_gas = 0
        for i in range(len(ce['particle_index']):
            Epot_part_to_gas = Epot_part_to_gas + Epot_particle[i]

    else:
        Epot_part_to_gas = 0.0
    
    # Computes the total energy of the gas in each cell
    # Step 1 + Step 2 + Step 3 + Step 4
    Etot = Etherm_gas + Epot_gas + Ekin_gas + Epot_part_to_gas
    # Locations of bound and unbound cells
    bound = Etot < 0.0
    unbound = Etot > 0.0
    # Creates a density array that is 0 where the gas is unbound and not background, the actual density where not
    density_bound = data['Zeros']
    density_bound[np.where(bound)] = data['Density'][np.where(bound)]
    # Bound minimum density
    bound_min_density = np.min(density_bound)
    
    return(bound_min_density)

def initialise_arrays():
    current_time = np.ndarray([0],dtype=float) #times

    cycle = np.ndarray([0],dtype=float) #top grid cycles

    Lp_x = np.ndarray([0],dtype=float)  #  particles x angular momentum
    Lp_y = np.ndarray([0],dtype=float)  #  particles y angular momentum
    Lp_z = np.ndarray([0],dtype=float)  #  particles z angular momentum

    Lg_x = np.ndarray([0],dtype=float)  #  gas x angular momentum
    Lg_y = np.ndarray([0],dtype=float)  #  gas y angular momentum
    Lg_z = np.ndarray([0],dtype=float)  #  gas z angular momentum

    Ltot_x = np.ndarray([0],dtype=float)  # tot x angular momentum
    Ltot_y = np.ndarray([0],dtype=float)  # tot y angular momentum
    Ltot_z = np.ndarray([0],dtype=float)  # tot z angular momentum

    return current_time, cycle, Lp_x, Lp_y, Lp_z, Lg_x, Lg_y, Lg_z, Ltot_x, Ltot_y, Ltot_z 

def open_file(file_name, append, particle_number):
    if (append == False): # Overwrite
        output_file = open(file_name, 'w' )

    # Write the first line of information in the file

        header = ("Time(yr), Cycle(#), Lp_x(gcm/s), Lp_y(gcm/s), Lp_z(gcm/s), Lg_x(gcm/s),"
                   "Lg_y(gcm/s), Lg_z(gcm/s), Ltot_x(gcm/s), Ltot_y(gcm/s), Ltot_z(gcm/s)," 
                   "Ltot(gcm/s), Lx_prim(gcm/s), Ly_prim(gcm/s), Lz_prim(gcm/s), ")
                  
        if particle_number >= 2:
            for i in range(particle_number - 1):
                header += ("Lx_comp_" + str(i+1) + "(gcm/s), " 
                           "Ly_comp_" + str(i+1) + "(gcm/s), "
                           "Lz_comp_" + str(i+1) + "(gcm/s), ")
        
        header += "\n"

        output_file.write(header)

        #output_file.write("Time(yr), Cycle(#), Lp_x(gcm/s), Lp_y(gcm/s), Lp_z(gcm/s), Lg_x(gcm/s),"
        #                  "Lg_y(gcm/s), Lg_z(gcm/s), Ltot_x(gcm/s), Ltot_y(gcm/s), Ltot_z(gcm/s)," 
        #                  "Ltot(gcm/s), Lx_prim, Ly_prim, Lz_prim"
        #                   +"\n")

    elif (append == True): # Append
        output_file = open(file_name, 'a' )

    return output_file

def ce_angular_momentum(directory, index, output_file):
    str_index = cef.index2str(index)    

    print("<------->")
    print("READING DIRECTORY " + str_index + ": ", directory)
    print("<------->")

    pf = load(directory)

    # Get length, time and mass units used in the simulation
    lu = pf.parameters['LengthUnits']
    tu = pf.parameters['TimeUnits']
    mu = pf.parameters['MassUnits']

    #  Import all of the data
    ce = pf.h.all_data()

    if (index == initial_path):
        global primary_threshold 
        primary_threshold = ce['BoundMinDensity']
        print("Primary Threshold: ", primary_threshold)

    primary = ce.cut_region(["grid['Density'] >" + str(primary_threshold)])

    if (angular_momentum_wrt_com == True):
        current_time = pf.current_time
        current_cycle = pf.parameters['InitialCycleNumber']
        
        # Create arrays to story the particle positions
        ppx = []
        ppy = [] 
        ppz = [] 
        # Create arrays to story the particle velocities
        pvx = [] 
        pvy = [] 
        pvz = []
        # Create aray to store particle mass
        pm = []

        pdex = ce['particle_index']
        print(ce["particle_position_x"], ce["particle_position_y"], ce["particle_position_z"])        


        # Particles
        for i in range(particle_number):
            # Particle Position (x,y,z)
            ppx.append(ce["particle_position_x"][i] * lu)
            ppy.append(ce["particle_position_y"][i] * lu)
            ppz.append(ce["particle_position_z"][i] * lu)
            # Particle Velocity (x,y,z)
            pvx.append(ce["particle_velocity_x"])
            pvy.append(ce["particle_velocity_y"])
            pvz.append(ce["particle_velocity_z"])
            # Particle Mass
            pm.append(ce["ParticleMass"][i]) 

        print(ppx, ppy, ppz)
        print(pm)

        primary_mass = np.max(pm)
        for i in range(particle_number):
            if pm[i] == primary_mass:
                primary_index = i
                break
            else:
                pass
        
        # Primary Envelope
        gx = primary["x"] * lu
        gy = primary["y"] * lu
        gz = primary["z"] * lu

        gvx = primary["x-velocity"]
        gvy = primary["y-velocity"]
        gvz = primary["z-velocity"]

        mass = primary["CellMass"]
           
        # Get the COM coords
        com_coords = ce.quantities["CenterOfMass"](use_particles=True) * lu
        # Move to the COM of the system (COM has 0 velocity).
        ppx_com = [x - com_coords[0] for x in ppx]
        ppy_com = [y - com_coords[0] for y in ppy]
        ppz_com = [z - com_coords[0] for z in ppz]

        gx_com = gx - com_coords[0]
        gy_com = gy - com_coords[1]
        gz_com = gz - com_coords[2]

        # Angular Momentum wrt COM
        # PARTICLES!
        # Particle Angular Momentum Components
        pLx = [pm[i]*(ppy_com[i]*pvz[i] - ppz_com[i]*pvy[i]) for i in range(len(pm))]
        pLy = [pm[i]*(ppz_com[i]*pvx[i] - ppx_com[i]*pvz[i]) for i in range(len(pm))]
        pLz = [pm[i]*(ppx_com[i]*pvy[i] - ppy_com[i]*pvx[i]) for i in range(len(pm))]

        # Total Angular momentum for each particle
        pL = [(pLx[i]**2 + pLy[i]**2 + pLz[i]**2)**0.5 for i in range(len(pm))]

        # Total Angular momentum for all particles
        pLtot = np.sum(pL)

        # GAS!
        # Gas Angular Momentum Components
        gLx = mass * (gy_com * gvz - gz_com * gvy)
        gLy = mass * (gz_com * gvx - gx_com * gvz)
        gLz = mass * (gx_com * gvy - gy_com * gvx)

        # Total Angular momentum for gas
        gLtot = np.sum((gLx**2 + gLy**2 + gLz**2)**0.5)
       
        current_Lp_x = np.sum(pLx)
        current_Lp_y = np.sum(pLy)
        current_Lp_z = np.sum(pLz)
        current_Lg_x = np.sum(gLx)
        current_Lg_y = np.sum(gLy)
        current_Lg_z = np.sum(gLz)
        current_Ltot_x = current_Lp_x + current_Lg_x
        current_Ltot_y = current_Lp_y + current_Lg_y
        current_Ltot_z = current_Lp_z + current_Lg_z
        current_Ltot = (current_Ltot_x**2 +current_Ltot_y**2 + current_Ltot_z**2)**0.5


    elif (angular_momentum_wrt_com == False):
        current_time = pf.current_time
        current_cycle = pf.parameters['InitialCycleNumber']

        #current output angular momenta values
        current_Lp_x = np.sum(ce['ParticleSpecificAngularMomentumX'] * ce['ParticleMass'])
        current_Lp_y = np.sum(ce['ParticleSpecificAngularMomentumY'] * ce['ParticleMass'])
        current_Lp_z = np.sum(ce['ParticleSpecificAngularMomentumZ'] * ce['ParticleMass'])
        current_Lg_x = np.sum(primary['AngularMomentumX'])
        current_Lg_y = np.sum(primary['AngularMomentumY'])
        current_Lg_z = np.sum(primary['AngularMomentumZ'])
        current_Ltot_x = current_Lp_x + current_Lg_x
        current_Ltot_y = current_Lp_y + current_Lg_y
        current_Ltot_z = current_Lp_z + current_Lg_z
        current_Ltot = (current_Ltot_x**2 +current_Ltot_y**2 + current_Ltot_z**2)**0.5

    # Write the columns of data 

    data = (str(current_time)+" "+str(current_cycle)+" "+str(current_Lp_x)+" "
            +str(current_Lp_y)+" "+str(current_Lp_z)+" "+str(current_Lg_x)+" "
            +str(current_Lg_y)+" "+str(current_Lg_z)+" "+str(current_Ltot_x)+" "
            +str(current_Ltot_y)+" "+str(current_Ltot_z)+" "+str(current_Ltot)+" ")

    data += str(pL[primary_index])[1:-1] + " "

    for i in range(particle_number):
        if i != primary_index:
            data += str(pL[pdex[i]])[1:-1] + " "

    data +="\n"

    output_file.write(data)

#    output_file.write(str(current_time)+" "+str(current_cycle)+" "+str(current_Lp_x)+" "
#                +str(current_Lp_y)+" "+str(current_Lp_z)+" "+str(current_Lg_x)+" "
#                +str(current_Lg_y)+" "+str(current_Lg_z)+" "+str(current_Ltot_x)+" "
#                +str(current_Ltot_y)+" "+str(current_Ltot_z)+" "+str(current_Ltot)+" "
#                +str(pL[0])[1:-1]+" "+str(pL[1])[1:-1]+" "+str(pL[2])[1:-1]+" "+"\n")


if __name__ == "__main__":
    # Read the inlist file and return the values of the variables.
    if len(sys.argv) >= 2:
        inlist_path = sys.argv[1]

        (root_dir, exclude_dir, plot_dir, initial_path, 
         final_path_plus_one, output_file_name, 
         output_file_append, angular_momentum_wrt_com, 
         smoothing_length, particle_number) = read_inlist(inlist_path)
    
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
    output_file = open_file(output_file_name, output_file_append, particle_number)
    
    for index in range(initial_path, final_path_plus_one):
        ce_angular_momentum(root_dir_list[index],  index, output_file)

    #close the written file
    output_file.close()

    print(" ")
    print("<------->")
    print("FINISHED: CE ANGULAR MOMENTUM")
    print("<------->")
