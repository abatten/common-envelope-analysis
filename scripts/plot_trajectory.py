from __future__ import absolute_import, division, print_function, unicode_literals
from yt.mods import *
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import sys
import ConfigParser
from pylab import genfromtxt

import modules.cefunctions as cef

mylog.disabled = True

#txt_file = "/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/particles_position_velocity_new_test.txt"
txt_file = "/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim1/plots/particles_positions_velocitys.txt"
mat0 = genfromtxt(txt_file, skip_header=1)
AU = 1.496e13 #cm
km = 10000 #cm

# Each particle has 6 components x,y,z & vx,vy,vz
# Plus the first two columns of time and cycle.
# Hence there is 6n + 2 columns where n is the number of particles.
num_particles = (len(mat0[0,:]) - 2) / 6

fig, axes = plt.subplots(1, 1, sharex='col', sharey='row')

mat0[:, 8] = (mat0[:, 8] - mat0[0,2]) / AU
mat0[:, 10] = (mat0[:, 10] - mat0[0,2]) / AU

#mat0[:, 14] = (mat0[:, 14] - mat0[0,2]) / AU
#mat0[:, 16] = (mat0[:, 16] - mat0[0,2]) / AU

mat0[:, 2] = (mat0[:, 2] - mat0[0,2]) / AU
mat0[:, 4] = (mat0[:, 4] - mat0[0,4]) / AU

plt.ylabel("$\mathrm{Y-Position\ (AU)}$")
plt.xlabel(("$\mathrm{X-Position\ (AU)}$"))
#row[0].set_yticks(np.arange(mat0[10,6], 1.5*np.max(mat0[:, (index+1)*6])))
plt.plot(mat0[:, 8], mat0[:, 10], 'r-')
#plt.plot(mat0[:, 14], mat0[:, 16], 'g-')
plt.plot(mat0[:, 2], mat0[:, 4], 'b-')

plt.ylim(-3,3)
plt.xlim(-3,3)
fig.subplots_adjust(left=0.1, right=0.98, wspace=0, hspace=0)
plt.savefig("trajectory_jstaff_sim1_2_particles.png")
plt.show()
