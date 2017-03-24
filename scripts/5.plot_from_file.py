from yt.mods import *
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from pylab import genfromtxt
import ConfigParser
import argparse

import cefunctions as cef

def energy_plot(txt_file, smoothed, marked):
    data = genfromtxt(txt_file, skip_header=1);

    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    #smoothed = True
    outfilename = "energy_components_total"

    if smoothed:
        x = np.isnan(data[:,6])
        while True in x:
            # Find NaNs in gas potential and replace with average of previous and following value
            for i in range(len(data[:,6])-1):
                if np.isnan(data[i,6]) and  not np.isnan(data[i+1,6]): # If next number isn't nan, average.
                    data[i,6] = 0.5 * (data[i+1,6] + data[i-1,6])

                elif np.isnan(data[i,6]) & np.isnan(data[i+1,6]): # If next number is nan, use previous.
                    data[i,6] = data[i-1,6]
            # Recalculate total potential energy
            for i in range(len(data[:,3])):
                data[i,3] = data[i,6] + data[i,9]
            # Recalculate Total energy
            for i in range(len(data[:,3])):
                data[i,1] = data[i,2] + data[i,3] + data[i,4]

            x = np.isnan(data[:,6])

        outfilename = outfilename + "_smooth_2"

    ax = plt.subplot(111)

    plt.plot(data[:,0], data[:,1], 'black', linewidth=3.5, label = "Total");
    plt.plot(data[:,0], data[:,2], linewidth=1, label = "Total Kinetic");
    plt.plot(data[:,0], data[:,3], linewidth=1, label = "Total Potential");
    plt.plot(data[:,0], data[:,4], linewidth=1, label = "Total Thermal");
    #plt.plot(data[:,0], data[:,5], label = "Gas Kinetic");
    plt.plot(data[:,0], data[:,6], label = "Gas Potential");
    #plt.plot(data[:,0], data[:,7], label = "Gas Thermal");
    #plt.plot(data[:,0], data[:,8], label = "Particle Kinetic");
    plt.plot(data[:,0], data[:,9], label = "Particle Potential");
    plt.plot(data[:,0], data[:,10], label = "Part-Gas Potential")

    #plt.plot(data[:,0], data[:,12], label = "Particle Thermal");

    plt.ylim(-2e46,1.5e46)
#    plt.xlim(0,14.5)
    plt.xlabel('Time (yr)', fontsize=16)
    plt.ylabel("Energy (ergs)", fontsize=16)

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=1, color='k')

        outfilename = outfilename + "_marked"

    plt.legend(loc='best', frameon=False)
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()

def seperation_plot(txt_file, marked):
    data = genfromtxt(txt_file, skip_header=1);
    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]


    outfilename = "seperations"
    yr = 365.25 * 24 * 60 * 60
    Rsun = 6.957e10 #cm

    for i in range(len(data[0,:])-2): # Skip first two columns.
        plt.plot(data[:,0]/yr, data[:,i+2]/Rsun, linewidth=2);

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=2, color='k')

        outfilename = outfilename + "_marked"

    plt.ylim(0,500)
    plt.xlabel('Time (yr)', fontsize=16)
    plt.ylabel("Seperation " r'($\mathrm{R}_\odot$)', fontsize=16)
    #plt.legend(bbox_to_anchor=(0.7, 0.65),loc=9, frameon=False)
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()

def thermal_plot(txt_file, marked):
    data = genfromtxt(txt_file, skip_header=1);

    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename = "thermal_components"


    ax = plt.subplot(111)

    plt.plot(data[:,0], data[:,4], linewidth=1, label = "Total Thermal");
    plt.plot(data[:,0], data[:,12], linewidth=1, label = "Primary Thermal");
    plt.plot(data[:,0], data[:,13], linewidth=1, label = "Vacuum Thermal");

#    plt.ylim(-2e46,1.5e46)
#    plt.xlim(0,14.5)
    plt.xlabel('Time (yr)', fontsize=16)
    plt.ylabel("Energy (ergs)", fontsize=16)


    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=2, color='k')

        outfilename = outfilename + "_marked"

    plt.legend(loc='best', frameon=False)
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()

def mass_loss_plot(txt_file, marked):
    data = genfromtxt(txt_file, skip_header=1);

    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename = "mass_loss"

    ax1 = plt.subplot(111)
    ax2 = ax1.twinx()

    mass_loss_rate = []
    for i in range(len(data[:,1])-1):
        mass_diff = data[i,1] - data[i+1,1]
        time_diff = data[i+1,0] - data[i,0]

        mass_loss_rate.append(mass_diff/time_diff)

    ax1.plot(data[:,0], data[:,1], 'b', linewidth=1, label = "Mass");
    ax2.plot(data[:-1,0], mass_loss_rate,'r', linewidth=1, label = "Mass loss rate");

#    plt.ylim(-2e46,1.5e46)
    plt.xlim(0,14.5)
    ax1.set_xlabel("Time " + r"($\mathrm{yr}$)", fontsize=16)
    ax1.set_ylabel("Mass "+r"($M_\odot$)", fontsize=16)
    ax2.set_ylabel("Mass loss rate " +r"($M_\odot/\mathrm{yr}$)", fontsize=16)


    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=2, color='k')

        outfilename = outfilename + "_marked"

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='best', frameon=False)
    plt.tight_layout()
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()

def angular_momentum_plot():
    return None


if __name__ == "__main__":

    # Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy", help="Create a plot of the energies from the specified file", action="store_true")
    parser.add_argument("--smoothed", help="Create a smoothed energy plot by replacing NaN's with the average of the values on either side of the,", action="store_true")
    parser.add_argument("--marked", help="Create a marked energy plot by adding vertical lines at x positions. e.g --marked 0 0.1 2 3", type=float, nargs="+")  
    parser.add_argument("--seperation", help="Create a plot of the seperations from the specified file", action="store_true")
    parser.add_argument("--angularmomentum", help="Create a plot of the angular momentum from the specified file", action="store_true")
    parser.add_argument("--thermal", help="Create a plot of the thermal energies from the specified file", action="store_true")
    parser.add_argument("--massloss", help="Create a plot of the mass loss and the mass loss rate from the specified file", action="store_true")
    parser.add_argument("txts", help="The text files to read.", nargs='*')
    
    args = parser.parse_args()

    if args.energy:
        for i in range(len(args.txts)):
            energy_plot(args.txts[i], args.smoothed, args.marked)

    if args.seperation:
        for i in range(len(args.txts)):
            seperation_plot(args.txts[i], args.marked)

    if args.thermal:
        for i in range(len(args.txts)):
            thermal_plot(args.txts[i], args.marked)

    if args.massloss:
        for i in range(len(args.txts)):
            mass_loss_plot(args.txts[i], args.marked)

