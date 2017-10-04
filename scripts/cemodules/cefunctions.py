import numpy as np
import os
import sys

def distance(point_1, point_2, units=1):
    """
    Calculates the distance between two 3D vectors using pythagoras.

    Args:
        point_1 ([flt, flt, flt]): Coordinates of the first point
        point_2 ([flt, flt, flt]): The coordinates of the second point
        units (flt): Specify the unit conversions. 

    Retuns:
        distance (flt): The distance between point_1 and point_2.
    """

    distance = (((point_2[0]-point_1[0])*units)**2.0
              + ((point_2[1]-point_1[1])*units)**2.0
              + ((point_2[2]-point_1[2])*units)**2.0)**0.5
    
    return distance

def grav_pot(mass1, mass2, radius, smoothed, 
             smoothing_length=3.0, smallest_cell=0):
    """
    Calculates the gravitational potential energy between two objects.
    The potential energy can be smoothed or not. 
    
    For the smoothed potential see M. Ruffert 1993
    
    Args:
        mass1 (flt): Mass of the primary object in grams
        mass2 (flt): Mass of the secondary object in grams
        radius (flt): The distance between the objects in cm
        smoothed (bool): Use a smoothed potential
        smoothing_length (flt): The length scale for the smoothed 
                                potential.
        smallest_cell (flt): The physical size of the smallest cell.

    Returns:
        potential_energy (flt): The gravitational potential energy 
                                  between the two objects
    """

    if smoothed:
        grav_constant = 6.67e-8
        top = grav_constant * mass1 * mass2
        slsc = smoothing_length * smallest_cell
        exponent = (-radius**2.0)/(slsc**2.0)
        bottom = (radius**2.0 + (slsc**2.0 * np.exp(exponent)))**0.5
        potential_energy = - top / bottom

    elif not smoothed:
        grav_constant = 6.67e-8
        top = grav_constant * mass1 * mass2
        bottom = radius
        potential_energy = - top / bottom

    return potential_energy

def index2str(index_input, num_char=4, prepend_char='0'):
    """
    index2str converts an integer to a string to
    represent that integer.

    i.e. 0000, 0001, 0002 etc rather than 1,10,11,12, etc

    Args:
        index_input (int): The index to convert
        num_char (int): The number of characters
        prepend_char (str): The character to prepend at the start.

    Returns:
        index_str (str): The string of the index with the prepended 
                         characters. 

    Examples:
        >>> index2str(4)
            "0004"
        >>> index2str(105)
            "0105"
        >>> index2str(105, num_char=6)
            "000105"
        >>> index2str(40, prepend_char='t')
            "tt40"
    """

    index_str = str(index_input)
    num_to_prepend = num_char - len(index_str)
    new_str_index = []

    for i in range(num_to_prepend):
        new_str_index.append(prepend_char)

    new_str_index.append(index_str)
    index_str = ''.join(new_str_index)

    return(index_str)

def root_sort(root_dir, exclude=[]):
    """
    Finds the paths to all the enzo data, sorts them based on
    directory name and adds them to root_dir_list.

    Args:
        root_dir (str): Path name of the directory to sort
        exclude (list[str]): List of directorys to exclude from 
                             the sort.

    Returns:
        root_dir_list (list[str]): The complete list of common 
                                   envelope data files in root_dir 
                                   sorted.
    """
    print(" ")
    print("<-------------->")
    print("ROOT DIRECTORY " + " : " + root_dir)
    print("<-------------->")
    print(" ")
    print("SORTING ROOT DIRECTORY FILES")
    root_dir_list = []

    for root, dirs, files in os.walk(root_dir):
        if  (root.split("/")[-1] in exclude and 
             root.split("/")[-1] != ''):

            print("EXCLUDING: " + root)
        # Skip the direcories that are listed in exclude_dir
        dirs[:] = [d for d in dirs if d not in exclude]
        files[:] = []  #  Remove all misc files
        current_folder = root
        # We don't want the root directory!!
        if (current_folder != root_dir):
            # Cycles subfolders and files in the current sub-folder
            for sub_root, sub_dirs, sub_files in os.walk(root):
                # Sorts the files in the subfolder to have the file 
                # Pass to yt in position [0]
                sub_files.sort()
                # Appends path of the enzo target file to root_dir_list 
                root_dir_list.append(os.path.join(root, sub_files[0]))
    
    root_dir_list.sort()
    
    return root_dir_list

def primary_coords(ce, lu=1):
    prim_index = primary_index(ce)
    primary_coords = [ce["particle_position_x"][prim_index] * lu,
                      ce["particle_position_y"][prim_index] * lu,
                      ce["particle_position_z"][prim_index] * lu]

    return primary_coords


def primary_index(ce):
    """
    Finds the index of the primary core assuming it has the largest mass.
    """
    particle_masses = ce["ParticleMassMsun"]
    primary_mass = np.max(particle_masses)
    for i in range(len(particle_masses)):
        if particle_masses[i] == primary_mass:
            primary_index = i
            break
        else:
            pass

    return primary_index
