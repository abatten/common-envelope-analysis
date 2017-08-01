from yt.mods import *
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from pylab import genfromtxt
import ConfigParser
import argparse


txt_file1 = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/seperations_jstaff_sim2.txt' 
txt_file2 = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/seperations_jstaff_sim2_part2.txt'
#txt_file3 = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/seperations_jstaff_sim2_1particle_restart.txt'
txt_file3 = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/seperations_jstaff_1planet_4yr_restart_149.txt'

data1 = genfromtxt(txt_file1, skip_header=1)
data2 = genfromtxt(txt_file2, skip_header=1)
data3 = genfromtxt(txt_file3, skip_header=1)

plot_loc = txt_file1[: - len(txt_file1.split("/")[-1])]

outfile_name = 'seperations_jstaff_combined_2'
yr = 365.25 * 24 * 60 * 60
Rsun = 6.957e10

plt.plot(data1[:,0], data1[:,2]/Rsun, 'b-', linewidth=2, label=r"$\mathrm{Inner\ Planet}$")
plt.plot(data1[:,0], data1[0:,3]/Rsun,'r-', linewidth=2, label=r"$\mathrm{Outer\ Planet}$")
plt.plot(data2[:,0], data2[:,2]/Rsun, 'r--', linewidth=2, label=r"$\mathrm{Outer\ after\ 5.6yr\ merger}$")
plt.plot(data3[:,0], data3[:,2]/Rsun, 'r:', linewidth=2, label=r"$\mathrm{Outer\ after\ 4.0yr\ merger}$")

plt.xlim(0,15)
plt.xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
plt.ylabel(r"$\mathrm{Seperation}\ (\mathrm{R}_\odot)$", fontsize=20)
plt.legend(loc='best', frameon=False)
plt.tight_layout()
plt.savefig(plot_loc + outfile_name + "_latexed.png")
plt.show()


