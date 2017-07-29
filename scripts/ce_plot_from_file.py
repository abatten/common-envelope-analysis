from yt.mods import *
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from pylab import genfromtxt
import ConfigParser
import argparse

import cemodules.cefunctions as cef

def energy_plot(txt_file, smoothed, marked):
    data = genfromtxt(txt_file, skip_header=1);

    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    #smoothed = True
    outfilename = "energy_components_plot"

    if smoothed:
        nan_list = np.isnan(data[:,6])
        while True in nan_list:
            # Find NaNs in gas potential and replace with 
            # average of previous and following value
            for i in range(len(data[:,6])-1):
                # If next number isn't nan then average.
                if np.isnan(data[i,6]) and  not np.isnan(data[i+1,6]):
                    data[i,6] = 0.5 * (data[i+1,6] + data[i-1,6])
                # If next number is nan, use previous
                elif np.isnan(data[i,6]) & np.isnan(data[i+1,6]): 
                    data[i,6] = data[i-1,6]
            # Recalculate total potential energy
            for i in range(len(data[:,3])):
                data[i,3] = data[i,6] + data[i,9] + data[i,10]
            # Recalculate Total energy
            for i in range(len(data[:,3])):
                data[i,1] = data[i,2] + data[i,3] + data[i,4]

            nan_list = np.isnan(data[:,6])

        outfilename = outfilename + "_smooth"

    ax = plt.subplot(111)

    ergs_unit = 10.0**46.0  # Normalising unit

    plt.plot(data[:,0], data[:,1]/ergs_unit, "k", linewidth=3.5,
             label = r"$\mathrm{E_{Total}}$")
    plt.plot(data[:,0], data[:,2]/ergs_unit, "b-", 
             label = r"$\mathrm{K_{Total}}$")
    plt.plot(data[:,0], data[:,3]/ergs_unit, "g-", 
             label = r"$\phi_{\mathrm{Total}}$");
    plt.plot(data[:,0], data[:,4]/ergs_unit, "r-", 
             label = r"$\mathrm{U_{Total}}$")
    plt.plot(data[:,0], data[:,5]/ergs_unit, "cyan", 
             label = r"$\mathrm{K_{Gas}}$")
    plt.plot(data[:,0], data[:,8]/ergs_unit, "magenta", 
             label = r"$\mathrm{K_{Particles}}$")
    plt.plot(data[:,0], data[:,6]/ergs_unit, "gold", 
             label = r"$\phi_{\mathrm{GG}}$")
    plt.plot(data[:,0], data[:,9]/ergs_unit, "hotpink", 
             label = r"$\phi_{\mathrm{PP}}$")
    plt.plot(data[:,0], data[:,10]/ergs_unit, "purple", 
             label = r"$\phi_{\mathrm{PG}}$")

    plt.ylim(-5,1.5)  # Change to make look good
    plt.xlim(0,14.5)  # Time axis
    plt.xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
    plt.ylabel(r"$\mathrm{Energy}\ (10^{46}\ \mathrm{ergs})$", fontsize=20)

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=1, color='k')

        outfilename = outfilename + "_marked"
    plt.tight_layout()
    plt.legend(ncol=3, loc='lower left', frameon=False)
    plt.savefig(plot_loc + outfilename + ".png")
    plt.savefig(plot_loc + outfilename + ".svg") 
    plt.show()

def seperation_plot(txt_file, marked):
    data = genfromtxt(txt_file, skip_header=1);

    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]

    outfilename = "seperations_jstaff_sim2_part2"
    yr = 365.25 * 24 * 60 * 60
    Rsun = 6.957e10 #cm

    for i in range(len(data[0,:])-2): # Skip first two columns.
        plt.plot(data[:,0], data[:,i+2]/Rsun, linewidth=2);

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=2, color='k')

        outfilename = outfilename + "_marked"

#    plt.ylim(0,500)
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
    plt.plot(data[:,0], data[:,4], linewidth=1, label = "Total Thermal")
    plt.plot(data[:,0], data[:,12], linewidth=1, label = "Primary Thermal")
    plt.plot(data[:,0], data[:,13], linewidth=1, label = "Vacuum Thermal")

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
        # Calculate the change in mass
        mass_diff = data[i,1] - data[i+1,1]
        # Calculate the change in time
        time_diff = data[i+1,0] - data[i,0]
        # Mass-loss rate = dM/dt
        mass_loss_rate.append(mass_diff/time_diff)

    ax1.plot(data[:,0], data[:,1], 'b', 
             label = r"$\mathrm{M_{env}}$");
    ax2.plot(data[:-1,0], mass_loss_rate,'r', 
             label = r"$\mathrm{\dot{M}_{env}}$");

#    plt.ylim(-2e46,1.5e46)
    plt.xlim(0,14.5)
    ax1.set_xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
    ax1.set_ylabel(r"$\mathrm{Mass}\ (\mathrm{M}_\odot)$", fontsize=20)
    ax2.set_ylabel(r"$\mathrm{Mass\ loss\ rate}\ (\mathrm{M}_\odot/\mathrm{yr})$", fontsize=20)
    ax1.set_ylim(0.2,0.35)

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=2, color='k')

        outfilename = outfilename + "_marked"

    h1, l1 = ax1.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1+h2, l1+l2, loc='best', prop={'size':18}, frameon=False)
    plt.tight_layout()
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()

def angular_momentum_plot(txt_file):
    mat0 = genfromtxt(txt_file, skip_header=1)
    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename = "angular_momentum_test_run"

    plt.plot(mat0[:,0], mat0[:,11],'k', linewidth=3, label = "Total $L$")
    plt.plot(mat0[:,0], mat0[:,2], 'r--', linewidth=2, label = "Particle $L_x$")
    plt.plot(mat0[:,0], mat0[:,3], 'r-.', linewidth=1,label = "Particle $L_y$")
    plt.plot(mat0[:,0], mat0[:,4], 'r:', linewidth=2, label = "Particle $L_z$")
    plt.plot(mat0[:,0], mat0[:,5], 'b--', linewidth=2, label = "Gas $L_x$")
    plt.plot(mat0[:,0], mat0[:,6], 'b-.', linewidth=1, label = "Gas $L_y$")
    plt.plot(mat0[:,0], mat0[:,7], 'b:', linewidth=2, label = "Gas $L_z$")
    #plt.plot(mat0[:,0], mat0[:,8], label = "Total $L_x$")
    #plt.plot(mat0[:,0], mat0[:,9], label = "Total $L_y$")
    #plt.plot(mat0[:,0], mat0[:,10], label = "Total $L_z$")

    #plt.ylim(-0.3e47,1e47)
    plt.xlabel('$\mathrm{Time (yr)}$', fontsize=16)
    plt.ylabel("$L\ (g\ cm\ s^{-1})$", fontsize=16)
#    plt.legend(loc='best', frameon=False); 
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show();

def position_and_velocities_plot(txt_file):
    mat0 = genfromtxt(txt_file, skip_header=1)
    AU = 1.496e13 #cm
    km = 10000 #cm
    # Each particle has 6 components x,y,z & vx,vy,vz
    # Plus the first two columns of time and cycle.
    # Hence there is 6n + 2 columns where n is the number of particles.
    num_particles = (len(mat0[0,:]) - 2) / 6

    fig, axes = plt.subplots(2, num_particles, sharex='col', sharey='row')

    irow=0
    for row in axes:
        for index in range(num_particles): 
            if irow == 0:
                mat0[:, (index+1)*6] = ((mat0[:, (index + 1) * 6] -
                                         mat0[0, (index + 1 )* 6]) / AU)
                row[0].set_ylabel("$\mathrm{Z-displacement\ (AU)}$")
                #row[0].set_yticks(np.arange(mat0[10,6], 1.5*np.max(mat0[:, (index+1)*6])))

            elif irow == 1:
                mat0[:, (index+1)*6+irow] = mat0[:, (index+1)*6+irow] / km
                row[0].set_ylabel("$\mathrm{Z-Velocity\ (km/s)}$")
            
            row[index].plot(mat0[:,0], mat0[:,(index+1)*6 + irow])
            
            row[index].set_xlabel("$\mathrm{Time\ (yr)}$")
            row[index].set_xticks(np.arange(mat0[0,0], mat0[-1,0], 3))
#            row[0].set_ylabel("$\mathrm{Z-displacement\ (cm)}$")
#            row[0].set_ylabel("$\mathrm{Z-Velocity\ (km/s)}$")
        irow+=1

    fig.subplots_adjust(left=0.1, right=0.98, wspace=0, hspace=0)
    plt.savefig("testing2.png")
    plt.show();



if __name__ == "__main__":
# Command line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("--energy", 
                       help="Create a plot of the energies from the \
                       specified file", action="store_true")
    parser.add_argument("--smoothed",
                       help="Create a smoothed energy plot by replacing NaN's \
                       with the average of the values on either side,", 
                       action="store_true")
    parser.add_argument("--marked", 
                        help="Create a marked energy plot by adding vertical \
                        lines at x positions. e.g --marked 0 0.1 2 3", 
                        type=float, nargs="+")  
    parser.add_argument("--seperation", 
                        help="Create a plot of the seperations from the \
                        specified file", action="store_true")
    parser.add_argument("--angularmomentum", 
                        help="Create a plot of the angular momentum \
                        from the specified file", action="store_true")
    parser.add_argument("--thermal", 
                        help="Create a plot of the thermal energies from the \
                        specified file", action="store_true")
    parser.add_argument("--massloss", 
                        help="Create a plot of the mass loss & the mass loss \
                        rate from the specified file", action="store_true")
    parser.add_argument("--posvel", 
                        help="Create plots of the positions and velocities \
                        of all the particles", action="store_true")
    parser.add_argument("txts", help="The text files to read.")


    
    args = parser.parse_args()

    if args.energy:
        energy_plot(args.txts, args.smoothed, args.marked)

    if args.seperation:
        seperation_plot(args.txts, args.marked)

    if args.thermal:
        thermal_plot(args.txts, args.marked)

    if args.massloss:
        mass_loss_plot(args.txts, args.marked)

    if args.angularmomentum:
         angular_momentum_plot(args.txts)

    if args.posvel:
         position_and_velocities_plot(args.txts)

