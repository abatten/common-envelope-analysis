from yt.mods import *
from pylab import genfromtxt
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset


import numpy as np
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
    plt.plot(data[:,0], data[:,2]/ergs_unit, "#377eb8", 
             label = r"$\mathrm{K_{Total}}$")
    plt.plot(data[:,0], data[:,3]/ergs_unit, "#4daf4a", 
             label = r"$\phi_{\mathrm{Total}}$");
    plt.plot(data[:,0], data[:,4]/ergs_unit, "#e41a1c", 
             label = r"$\mathrm{U_{Total}}$")
    plt.plot(data[:,0], data[:,5]/ergs_unit, "#f781bf", 
             label = r"$\mathrm{K_{Gas}}$")
    plt.plot(data[:,0], data[:,8]/ergs_unit, "#984ea3", 
             label = r"$\mathrm{K_{Particles}}$")
    plt.plot(data[:,0], data[:,6]/ergs_unit, "#CCCC00", 
             label = r"$\phi_{\mathrm{GG}}$")
    plt.plot(data[:,0], data[:,9]/ergs_unit, "#ff7f00", 
             label = r"$\phi_{\mathrm{PP}}$")
    plt.plot(data[:,0], data[:,10]/ergs_unit, "#a65628", 
             label = r"$\phi_{\mathrm{PG}}$")

    plt.ylim(-5,1.5)  # Change to make look good
    plt.xlim(5.6,23.5)  # Time axis
    plt.xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
    plt.ylabel(r"$\mathrm{Energy}\ (10^{46}\ \mathrm{ergs})$", fontsize=20)

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=1, color='k')

        outfilename = outfilename + "_marked"
    plt.tight_layout()
    plt.legend(ncol=3, loc='lower left', frameon=False)
    plt.savefig(plot_loc + outfilename + ".png")
    plt.savefig(plot_loc + outfilename + ".pdf")
    plt.savefig(plot_loc + outfilename + ".eps")
    plt.show()

def seperation_plot(txt_file, marked):
    data = genfromtxt(txt_file, skip_header=1);

    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]

    outfilename = "seperations_plot"
    yr = 365.25 * 24 * 60 * 60
    Rsun = 6.957e10 #cm

    colours = ["b", "r", "g"]

    fig, ax = plt.subplots(figsize=[8,6])

    for i in range(len(data[0,:])-2): # Skip first two columns.
        ax.plot(data[:,0], data[:,i+2]/Rsun, colours[i], linewidth=2);

    if marked:
        for i in range(len(marked)):
            plt.axvline(x=marked[i], ymin=0, ymax=1, linewidth=2, color='k')

        outfilename = outfilename + "_marked"

    #plt.ylim(0,2600)
    plt.xlabel(r"$\mathrm{Time\ (yr)}$", fontsize=20)
    plt.ylabel(r"$\mathrm{Seperation}\ (\mathrm{R}_\odot$)", fontsize=20)
    #plt.legend(bbox_to_anchor=(0.7, 0.65),loc=9, frameon=False)
    
#    axins = zoomed_inset_axes(ax, 6, loc=4)
#    axins.semilogy(data[:,0], data[:,2]/Rsun, colours[0], linewidth=1)

#    x1,x2,y1,y2 = 0.7, 1.5, 0, 100

#    axins.set_xlim(x1, x2)
#    axins.set_ylim(y1, y2)

#    plt.xticks(visible=False)
#    plt.yticks(visible=False)

    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
#    mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")

#    plt.draw()
    plt.savefig(plot_loc + outfilename + ".png")
    plt.savefig(plot_loc + outfilename + ".eps")
    plt.savefig(plot_loc + outfilename + ".pdf")
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
#    plt.xlim(0,14.5)
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
    plt.savefig(plot_loc + outfilename + ".pdf")
    plt.savefig(plot_loc + outfilename + ".eps")
    plt.show()

def angular_momentum_plot(txt_file):
    data0 = genfromtxt(txt_file, skip_header=1)
    # Find the location of the input file to save plot to.
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename = "angular_momentum_test_run"

    plt.plot(data0[:,0], data0[:,11],'k', linewidth=3, label = "Total $L$")
    plt.plot(data0[:,0], data0[:,2], 'r--', linewidth=2, label = "Particle $L_x$")
    plt.plot(data0[:,0], data0[:,3], 'r-.', linewidth=1,label = "Particle $L_y$")
    plt.plot(data0[:,0], data0[:,4], 'r:', linewidth=2, label = "Particle $L_z$")
    plt.plot(data0[:,0], data0[:,5], 'b--', linewidth=2, label = "Gas $L_x$")
    plt.plot(data0[:,0], data0[:,6], 'b-.', linewidth=1, label = "Gas $L_y$")
    plt.plot(data0[:,0], data0[:,7], 'b:', linewidth=2, label = "Gas $L_z$")
    #plt.plot(data0[:,0], data0[:,8], label = "Total $L_x$")
    #plt.plot(data0[:,0], data0[:,9], label = "Total $L_y$")
    #plt.plot(data0[:,0], data0[:,10], label = "Total $L_z$")

    #plt.ylim(-0.3e47,1e47)
    plt.xlabel('$\mathrm{Time (yr)}$', fontsize=16)
    plt.ylabel("$L\ (g\ cm\ s^{-1})$", fontsize=16)
#    plt.legend(loc='best', frameon=False); 
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()

def position_and_velocities_plot(txt_file):
    data0 = genfromtxt(txt_file, skip_header=1)

    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename = "positions_and_velocities"

    AU = 1.496e13 #cm
    km = 10000 #cm
    # Each particle has 6 components x,y,z & vx,vy,vz
    # Plus the first two columns of time and cycle.
    # Hence there is 6n + 2 columns where n is the number of particles.
    num_particles = (len(data0[0,:]) - 2) / 6

    fig, axes = plt.subplots(2, num_particles, sharex='col', sharey='row')

    irow=0
    for row in axes:
        for index in range(num_particles): 
            if irow == 0:
                data0[:, (index+1)*6] = ((data0[:, (index + 1) * 6] -
                                         data0[0, (index + 1 )* 6]) / AU)
                row[0].set_ylabel("$\mathrm{Z-displacement\ (AU)}$")
                #row[0].set_yticks(np.arange(data0[10,6], 1.5*np.max(data0[:, (index+1)*6])))

            elif irow == 1:
                data0[:, (index+1)*6+irow] = data0[:, (index+1)*6+irow] / km
                row[0].set_ylabel("$\mathrm{Z-Velocity\ (km/s)}$")
            
            row[index].plot(data0[:,0], data0[:,(index+1)*6 + irow])
            
            row[index].set_xlabel("$\mathrm{Time\ (yr)}$")
            row[index].set_xticks(np.arange(data0[0,0], data0[-1,0], 3))
#            row[0].set_ylabel("$\mathrm{Z-displacement\ (cm)}$")
#            row[0].set_ylabel("$\mathrm{Z-Velocity\ (km/s)}$")
        irow+=1

    fig.subplots_adjust(left=0.1, right=0.98, wspace=0, hspace=0)
    plt.savefig(plot_dir + outfilename + ".png")
    plt.show()

def gravodrag_plot(txt_file):
    data0 = genfromtxt(txt_file, skip_header=1)
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename = "gravodrag_plot"

    plt.plot(data0[:,0], data0[:,2],'b', linewidth=1, label = "Outer Particle")
    #plt.plot(data0[:,0], data0[:,15],'g', linewidth=1, label = "Inner Particle")

    plt.ylim(1e31,2.5e33)
    plt.xlim(0,14.5)
    plt.xlabel('$\mathrm{Time (yr)}$', fontsize=20)
    plt.ylabel(r"$\mathrm{Gravodrag\ (dynes)}$", fontsize=20)
    plt.legend(loc='best', frameon=False)
    plt.tight_layout()
    plt.savefig(plot_loc + outfilename + ".png")
    plt.show()


def centre_of_mass_plot(txt_file):
    data0 = genfromtxt(txt_file, skip_header=1)
    data1 = genfromtxt("/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/centre_of_mass.txt", skip_header=1)
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfilename1 = "centre_of_mass_xy_plot"
    outfilename2 = "centre_of_mass_time_plot"

    #plt.plot(data0[:, 0], data0[:, 2] - 0.5, "b:", label="COMx")
    #plt.plot(data0[:, 0], data0[:, 3] - 0.5, "r:", label="COMy")
    #plt.plot(data0[:, 0], data0[:, 4] - 0.5, "g:", label="COMz")
    plt.title("Centre of Mass")
    plt.plot(data0[:,2]-0.5, data0[:,3]-0.5, "b", label="After Reduced Mass")
    plt.plot(data0[0,2]-0.5, data0[0,3]-0.5, "bp")
    plt.plot(data0[-1,2]-0.5, data0[-1,3]-0.5, "b+")
    plt.plot(data1[:57,2]-0.5, data1[:57,3]-0.5, "g", label="Before Reduced Mass")
    plt.plot(data1[0,2]-0.5, data1[0,3]-0.5, "gp")
    plt.plot(data1[56,2]-0.5, data1[56,3]-0.5, "g+")
    plt.legend(loc="best", frameon=False)
    plt.xlim(-0.3,-0.05)
    plt.ylim(-0.16,-0.05)
    plt.xlabel("X-cood")
    plt.ylabel("Y-Cood")
    plt.savefig(plot_loc + outfilename1 + ".png")
    plt.show()

    plt.plot(data0[:, 0], data0[:, 4] - 0.5, "g:", linewidth=2, label="COMz-After")
    plt.plot(data0[:, 0], data0[:, 3] - 0.5, "r:", linewidth=2, label="COMy-After")
    plt.plot(data0[:, 0], data0[:, 2] - 0.5, "b:", linewidth=2, label="COMx-After")
    plt.plot(data1[:57, 0], data1[:57, 4] - 0.5, "g", linewidth=2, label="COMz-Before")
    plt.plot(data1[:57, 0], data1[:57, 3] - 0.5, "r", linewidth=2, label="COMy-Before")
    plt.plot(data1[:57, 0], data1[:57, 2] - 0.5, "b", linewidth=2, label="COMx-Before")
    
    plt.ylabel("Coordinate of Centre of Mass")
    plt.xlabel(r"$\mathrm{Time\ (yr)}$")
    plt.legend(loc="best", frameon=False)
    plt.savefig(plot_loc +outfilename2 + ".png")
    plt.show()


def gas_direction_plot(txt_file):
    data = genfromtxt(txt_file, skip_header=1)
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfile_name = "gas_direction_plot_with_Particles" 

    km = 10**5

    plt.plot(data[:, 0], data[:, 2] / km, "k:", linewidth=2, label=r"$\mathrm{Inner}\ |v_g|$")
    plt.plot(data[:, 0], data[:, 3] / km, "b:", label=r"$\mathrm{Inner}\ v_{\parallel}$")
    plt.plot(data[:, 0], data[:, 4] / km, "r:", label=r"$\mathrm{Inner}\ v_{\bot}$")
    #plt.plot(data[:, 0], data[:, 7] / km , "g:", label=r"$\mathrm{Inner}\ |v_p|}$")
    plt.plot(data[:, 0], data[:, 8] / km, "k", linewidth=2, label=r"$\mathrm{Outer}\ |v_g|$") 
    plt.plot(data[:, 0], data[:, 9] / km, "b", label=r"$\mathrm{Outer}\ v_{\parallel}$") 
    plt.plot(data[:, 0], data[:, 10] / km, "r", label=r"$\mathrm{Outer}\ v_{\bot}$") 
    plt.plot(data[:, 0], data[:, 13] / km, "g", linewidth=2, label=r"$\mathrm{Outer}\ |v_p|}$")

    plt.ylabel(r"$\mathrm{Velocity\ \left(km\ s^{-1}\right)}$", fontsize=20)
    plt.xlabel(r"$\mathrm{Time\ \left( yr \right)}$", fontsize=20)
    plt.xlim(0,14.3)
    plt.legend(loc="best", frameon=False)
    plt.savefig(plot_loc + outfile_name + ".png")
    plt.savefig(plot_loc + outfile_name + ".pdf")
    plt.savefig(plot_loc + outfile_name + ".eps")
    plt.show()


def orbit_params_plot(txt_file):
    data = genfromtxt(txt_file, skip_header=1)
    plot_loc = txt_file[: - len(txt_file.split("/")[-1])]
    outfile_name = "orbital_parameters_plot"

    f, axarr = plt.subplots(2, sharex=True)
    
    
    
    #axarr[0].plot(x, y)
    #axarr[0].set_title('Sharing X axis')
    #axarr[1].scatter(x, y)

    #ax1 = plt.subplot(111)
    #ax2 = plt.subplot(111)

    semi_major_scale = 10**14
   
    axarr[0].plot(data[:, 0], data[:, 2] / semi_major_scale, "b", label=r"$\mathrm{Inner\ Semi-Major\ Axis}$")
    axarr[0].plot(data[:, 0], data[:, 11] / semi_major_scale, "r", label=r"$\mathrm{Outer\ Semi-Major\ Axis}$")

    axarr[1].plot(data[:,0], data[:,4], "b", label=r"$\mathrm{Inner\ Eccentricity}$")
    axarr[1].plot(data[:,0], data[:,13], "r", label=r"$\mathrm{Outer\ Eccentricity}$")

    axarr[1].set_xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
    axarr[0].set_ylabel(r"$\mathrm{a}\ (\mathrm{10^{14}\ cm})$", fontsize=20)
    axarr[1].set_ylabel(r"$\mathrm{e}$", fontsize=20)
    axarr[0].set_ylim(1/1000,1)
    axarr[0].set_xlim(0,14.3)
    axarr[1].set_xlim(0,14.3)

    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)

    h1, l1 = axarr[0].get_legend_handles_labels()
    h2, l2 = axarr[1].get_legend_handles_labels()

    axarr[0].legend(h1, l1, loc='best', prop={'size':14}, frameon=False)
    axarr[1].legend(h2, l2, loc='best', prop={'size':14}, frameon=False)

    plt.show()

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
    parser.add_argument("--gravodrag",
                        help="Create plot of gravitational drag", action="store_true")
    parser.add_argument("--gasdir", 
                        help="Creates plot of gas direction", action="store_true")
    parser.add_argument("--orbitparams",
                        help="Creates plot of orbital parameters", action="store_true")
    parser.add_argument("--centreofmass", help="Create plot of center of mass", action="store_true")
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

    if args.gravodrag:
        gravodrag_plot(args.txts)

    if args.centreofmass:
        centre_of_mass_plot(args.txts)

    if args.gasdir:
        gas_direction_plot(args.txts)

    if args.orbitparams:
        orbit_params_plot(args.txts)
