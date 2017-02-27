import numpy as np
import os
import sys

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

def grav_pot(mass1, mass2, radius, smoothed, smoothing_length=3.0, smallest_cell=0):
    """
    Calculates the gravitational potential between two objects.
    The potential can be smoothed or not. For the smoothed potential see M. Ruffert 1993
    
    mass1: Float
    mass2: Float
    radius: Float
    smoothed: Boolean
    smoothing_lengtn: Float
    smallest_cell: Float
    """
    if smoothed:
        grav_constant = 6.67e-8
        top = grav_constant * mass1 * mass2
        slsc = smoothing_length * smallest_cell
        bottom = (radius**2.0 + (slsc**2.0 * np.exp((-radius**2.0)/(slsc**2.0))))**0.5
        potential = - top / bottom

    elif not smoothed:
        grav_constant = 6.67e-8
        top = grav_constant * mass1 * mass2
        bottom = radius
        potential = - top / bottom

    return potential

def index2str(index_input, num_digits=4, append_character='0'):
        '''
        index2str converts an integer to a string to
         represent that integer.

        This is used so that the output files will be sorted in the
        correct order.

        i.e. 0000, 0001, 0002 etc rather than 1,10,11,12, etc

        index_input: Expects an integer.
        num_digits: Expects an integer.
        append_character: Expects a string.


        Examples:
        index2str(4) returns "0004"
        index2str(105) returns "0105"
        index2str(105, num_digits=6) returns "000105"
        index2str(40, append_character='t') returns "tt40"
        '''

        index_str = str(index_input)
        num_to_append = num_digits - len(index_str) # How many characters to add.
        new_str_index = [] # Stores the digits of the new index

        for i in range(num_to_append):
                new_str_index.append(append_character)

        new_str_index.append(index_str)
        index_str = ''.join(new_str_index)

        return(index_str)

def root_sort(root_dir, exclude=[]):
    """
    Finds the paths to all the enzo data, sorts them based on
    directory name and adds them to root_dir_list.


    root_dir: Expects a path name
    exclude: Expects a list of strings.
    """
    print(" ")
    print("<-------------->")
    print("ROOT DIRECTORY " + " : " + root_dir)
    print("<-------------->")
    print(" ")
    print("SORTING ROOT DIRECTORY FILES")
    root_dir_list = []

    for root, dirs, files in os.walk(root_dir):
        if  (root.split("/")[-1] in exclude) and (root.split("/")[-1] != ''):
            print("EXCLUDING: " + root)

        # Skip the direcories that are listed in exclude_dir
        dirs[:] = [d for d in dirs if d not in exclude]
        files[:] = []
        current_folder = root
        # os.walk includes the root directory.
        # We don't want the root directory!!
        if (current_folder != root_dir):
            # cycles subfolders and files in the current sub-folder
            for sub_root, sub_dirs, sub_files in os.walk(root):
                # sorts the files in the subfolder to have the file to pass to yt in position [0]
                sub_files.sort()
                # Appends the path of the enzo target file to root_dir_list 
                root_dir_list.append(os.path.join(root, sub_files[0]))
  # sorts list by directory name
    root_dir_list.sort()
    return root_dir_list

