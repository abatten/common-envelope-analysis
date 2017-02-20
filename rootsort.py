import os
import sys

def root_sort(root_dir, exclude=['NONE']):
    """
    Finds the paths to all the enzo data, sorts them based on
    directory name and adds them to root_dir_list.


    root_dir_list: Expects a list of path names
    exclude: Expects a list of strings.
    """
    print(" ")
    print("<-------------->")
    print("ROOT DIRECTORY " + " : " + root_dir)
    print("<-------------->")
    print(" ")
    print("SORTING ROOT DIRECTORY FILES")
    print(root_dir)
    root_dir_list = []
    # Cycles subfolders and files in the root directory
    for root, dirs, files in os.walk(root_dir):
        # Skip the direcories that are listed in exclude_dir
        for i in range(len(exclude)):
            if  exclude[i] not in root.split("/"):
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
            elif exclude[i] in root.split("/"):
                print("EXCLUDING: " + root)
    # sorts list by directory name
    root_dir_list.sort()
    return root_dir_list
