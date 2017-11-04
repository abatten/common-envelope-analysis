from __future__ import absolute_import, division, print_function

import yt.mods as yt

import numpy as np
import sys
import ConfigParser

import cemodules.cefunctions as cef
import cemodules.constants as CONST

def read_inlist(ipath):
    print(" ")
    print("<-------------->")
    print("READING INLIST FILE")
    print("<-------------->")
    inlist_name = ipath.split("/")[-1]
    config = ConfigParser.ConfigParser()
    config.readfp(open(inlist_name, "r"))

    # Read in the config file
    root_dir = config.get("Common", "root_dir")
    exclude_dir = config.get("Common", "exclude_dir")
    plot_dir = config.get("Common", "plot_dir") 
    initial_path = config.getint("Common", "initial_path")
    final_path_plus_one = config.getint("Common", "final_path_plus_one")
    particle_number = config.getint("Common", "particle_number")
    output_file_name = config.get("GasDirection", "output_file_name")

    print("INLIST FILE: " + inlist_name)
    print("ROOT DIRECTORY: " + str(root_dir))
    print("EXCLUDING DIRECTORIES: " + str(exclude_dir)) 
    print("OUTPUT DIRECTORY: " + str(plot_dir)) 
    print("INITIAL PATH: " + str(initial_path))
    print("FINAL PATH PLUS ONE: " + str(final_path_plus_one))
    print("OUTPUT FILE NAME: " + str(output_file_name)) 

    return (root_dir, exclude_dir, plot_dir, initial_path,
            final_path_plus_one, output_file_name, particle_number)


def open_file(file_name, num_particles):
    """
    Creates file and writes the header information
    """
    output_file = open(file_name, "w")

    header = "Time (yr)  Cycle(#)"

    vel_mag_dict ={}
    rad_vel_dict = {}
    perp_vel_dict = {}
    density_dict = {}
    cosi_dict = {}
    particle_vel = {}

    for i in range(num_particles -1):
        vel_mag_dict[i+1] = "Part_%s_VelMag(cm/s)" % (i+1)
        rad_vel_dict[i+1] = "Part_%s_RadVel(cm/s)" % (i+1)
        perp_vel_dict[i+1] = "Part_%s_PerpVel(cm/s)" % (i+1)
        density_dict[i+1] = "Part_%s_Density(g/cm3)" % (i+1)
        cosi_dict[i+1] = "Part_%s_cos(i)" % (i+1)
        particle_vel[i+1] = "Part_%s_Velocity(cm/s)" % (i+1)

        header = "  ".join([header,
                            vel_mag_dict[i+1],
                            rad_vel_dict[i+1],
                            perp_vel_dict[i+1],
                            cosi_dict[i+1],
                            density_dict[i+1],
                            particle_vel[i+1]])

    output_file.write(header + "\n")
    return output_file


def gas_velocity_direction(directory, index, output_file):
    str_index = cef.index2str(index)

    print(" ")
    print("<-------------->")
    print("READ ROOTDIRECTORY " + str_index + ":", directory)
    print("<-------------->")

    pf = yt.load(directory)

    #  Gets the length, time and mass units used in the current simulation
    lu = pf.parameters["LengthUnits"]
    tu = pf.parameters["TimeUnits"]
    mu = pf.parameters["MassUnits"]

    current_cycle = pf.parameters["InitialCycleNumber"]

    year = 365.25 * 24 * 60 * 60
    current_time = pf.current_time * tu / year

    ce = pf.h.all_data()

    pdex = ce["particle_index"]
    prim_index = cef.primary_index(ce)
    prim_coords = np.array(cef.primary_coords(ce, lu))
    num_particles = len(ce["ParticleMassMsun"])

    #print("Prim Index:" , prim_index)
    cell_width = ((pf.h.grid_right_edge - pf.h.grid_left_edge) / pf.h.grid_dimensions)[0][0]
    cell_diag = np.sqrt(3) * cell_width


    #  List of offsets to average over surrounding cells
    offsets = [np.array([0, 0, 0]),
               np.array([-cell_diag, cell_diag, 0]),
               np.array([0, cell_width, 0]),
               np.array([cell_diag, cell_diag, 0]),
               np.array([cell_width, 0, 0]),
               np.array([cell_diag, -cell_diag, 0]),
               np.array([0, -cell_width, 0]),
               np.array([-cell_diag, -cell_diag, 0]),
               np.array([-cell_width, 0, 0]),
               np.array([0, 0, -cell_width]),
               np.array([-cell_diag, cell_diag, -cell_width]),
               np.array([0, cell_width, -cell_width]),
               np.array([cell_diag, cell_diag, -cell_width]),
               np.array([cell_width, 0, -cell_width]),
               np.array([cell_diag, -cell_diag, -cell_width]),
               np.array([0, -cell_width, -cell_width]),
               np.array([-cell_diag, -cell_diag, -cell_width]),
               np.array([-cell_width, 0, -cell_width]),
               np.array([0, 0, cell_width]),
               np.array([-cell_diag, cell_diag, cell_width]),
               np.array([0, cell_width, cell_width]),
               np.array([cell_diag, cell_diag, cell_width]),
               np.array([cell_width, 0, cell_width]),
               np.array([cell_diag, -cell_diag, cell_width]),
               np.array([0, -cell_width, cell_width]),
               np.array([-cell_diag, -cell_diag, cell_width]),
               np.array([-cell_width, 0, cell_width])]


    part_cosi = {}
    part_v_perp = {}
    part_v_orth = {}
    part_density = {}
    part_velmag = {}
    part_velocity = {}

    for i in range(num_particles):
        if i != prim_index:
            part_phys_coords = np.array([ce["particle_position_x"][i] * lu,
                                         ce["particle_position_y"][i] * lu,
                                         ce["particle_position_z"][i] * lu])

            part_code_coords = np.array([ce["particle_position_x"][i],
                               ce["particle_position_y"][i],
                               ce["particle_position_z"][i]])

            particle_velocity = np.array([ce["particle_velocity_x"][i],
                                      ce["particle_velocity_y"][i],
                                      ce["particle_velocity_z"][i]])
            particle_velocity_mag = np.linalg.norm(particle_velocity)

            velocity_list = []
            density_list = []
            for j in range(len(offsets)):
                density = pf.h.find_field_value_at_point(["Density"], 
                                                          part_code_coords + offsets[j])

                velocity = pf.h.find_field_value_at_point(["x-velocity", 
                                                           "y-velocity",
                                                           "z-velocity"],
                                                           part_code_coords + offsets[j])

                velocity_list.append(np.array(velocity))
                density_list.append(density)
            
            velocity_mean = np.mean(velocity_list, axis=0)
            density_mean = np.mean(density_list)
            #print("Density:   ", density_mean)
            
           # velocity_vect = np.array(velocity_mean)
            velocity_mag = np.linalg.norm(velocity_mean)
            velocity_unit = velocity_mean / velocity_mag

            dist_vect = part_phys_coords - prim_coords
            dist_mag = np.linalg.norm(dist_vect)
            dist_unit = dist_vect / dist_mag

            cosinei = np.dot(velocity_unit, dist_unit)
            angle = np.arccos(cosinei)
            sinei = np.sin(angle)

            part_velmag[pdex[i]] = velocity_mag
            part_cosi[pdex[i]] = cosinei
            part_v_orth[pdex[i]] = cosinei * velocity_mag
            part_v_perp[pdex[i]] = sinei * velocity_mag
            part_density[pdex[i]] = density_mean
            part_velocity[pdex[i]] = particle_velocity_mag
            print("ParticleVelocity: ", particle_velocity_mag)
            print("cos(i):           ", cosinei)
            print("sin(i):           ", sinei)
            print("Vel:              ", velocity_mean)
            print("Vel Mag:          ", velocity_mag)
            print("Vel Unit:         ", velocity_unit)
            print("PartCoord:        ", part_phys_coords)
            print("PrimCoord:        ", prim_coords)
            print("Dist:             ", dist_vect)
            print("Dist Mag:         ", dist_mag)
            print("Dist Unit:        ", dist_unit)
            print("Vel Orth:         ", pdex[i], part_v_orth[pdex[i]])
            print("Vel Perp:         ", pdex[i], part_v_perp[pdex[i]])
            print(" ")

    data = str(current_time) + " " + str(current_cycle)

#    print("Pdex: ", pdex)
#    print("Prim: ", prim_index)
    for i in range(len(part_cosi.items())):
#        print(i, "Printing")
#        if i != prim_index:
#        print("Printing Pdex[i]:", pdex[i])
        data = " ".join([data, str(part_velmag.items()[i][1]),
                         str(part_v_orth.items()[i][1]), str(part_v_perp.items()[i][1]), 
                         str(part_cosi.items()[i][1]), str(part_density.items()[i][1]),
                         str(part_velocity.items()[i][1])])
    output_file.write(data)
    output_file.write("\n")



if __name__ == "__main__":
    yt.mylog.disabled = True

    if len(sys.argv) >= 2:
        inlist_path = sys.argv[1]
    else:
        print("INLIST NOT SUPPLIED!!!")
        sys.exit(0)

    
    (root_dir, exclude_dir, plot_dir, initial_path,
    #for grid in pf.h.grids:
    #for grid in pf.h.grids:
     final_path_plus_one, output_file_name, 
     particle_number) = read_inlist(inlist_path)

    root_dir_list = cef.root_sort(root_dir, exclude=exclude_dir)

    output_file_name = plot_dir + output_file_name
    output_file = open_file(output_file_name, particle_number)

    for index in range(initial_path, final_path_plus_one):
        gas_velocity_direction(root_dir_list[index], index, output_file)
    output_file.close()        
