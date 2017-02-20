from __future__ import absolute_import, division, print_function, unicode_literals

from yt.mods import *

import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys

from index2str import index2str
from rootsort import root_sort

import ConfigParser

################################################################################
################################# User Controls ################################
################################################################################

config = ConfigParser.ConfigParser()
config.readfp(open(r'ce_analysis_inlist'))

root_dir = config.get('Energy Section', 'root_dir')
exclude_dir = config.get('Energy Section', 'exclude_dir')
plot_dir = config.get('Energy Section', 'plot_dir')
initial_path = config.get('Energy Section', 'initial_path')
final_path_plus_one = config.get('Energy Section', 'final_path_plus_one')
output_file_name = config.get('Energy Section', 'output_file_name')
output_file_append = config.get('Energy Section', 'output_file_append')
use_smoothed_potential = config.get('Energy Section', 'use_smoothed_potential')
smoothing_length = config.get('Energy Section', 'smoothing_length')


print(root_dir)
print(exclude_dir)
print(plot_dir)
print(initial_path)
print(final_path_plus_one)
print(output_file_name)
print(output_file_append)
print(use_smoothed_potential)
print(smoothing_length)
################################################################################
################################################################################
################################################################################
