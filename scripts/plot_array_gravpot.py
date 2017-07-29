from yt.mods import *
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid


plot_dir = '/disks/ceres/makemake/acomp/abatten/jstaff/sim2/plots/'
root_dir = '/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj/'
directory_list = ['DD0000/CE0008','DD0005/CE0005', 'DD0015/CE0015', 'DD0020/CE0020', 'DD0025/CE0025','DD0030/CE0030', 'DD0035/CE0035', 'DD0060/CE0060', 'DD0140/CE0140']

fig = plt.figure()

grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (3, 3),
                axes_pad = 0.05,
                label_mode = "L",
                share_all = True,
                cbar_location="right",
                cbar_mode="single",
                cbar_size="3%",
                cbar_pad="0%")


for i, directory_list in enumerate(directory_list):
    # Load the data and create a single plot
    ds = load(root_dir + directory_list) # load data
    p = SlicePlot(ds, 'z', 'Grav_Potential')
    p.annotate_particles(1.0, p_size=25.0, marker='o', col='black')
    # Ensure the colorbar limits match for all plots
    p.set_zlim('Grav_Potential', -9e-1, 0)

    # This forces the ProjectionPlot to redraw itself on the AxesGrid axes.
    plot = p.plots['Grav_Potential']
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    p._setup_plots()

plt.savefig(plot_dir+'multiplot_3x3_Grav_Potential_array.png')



#pf = load('/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj/DD0005/CE0005')
#common_envelope = pf.h.all_data()
# Create binned profiles to evaluate velocity, density, pressure, temperature as a function of radius

#sp = SlicePlot(pf, 'z', "Density", width = 1.0) #centered in the center of the axis
#sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particles
#sp.annotate_text((0.7,1.05), "time = " + str(pf.current_time) + "yr") #annotates the current simulation time on the plot
#   sp.annotate_image_line((0.5, 0.0), (0.5, 1.0), plot_args={'linewidth':1,'color':'black'})
#   sp.annotate_image_line((0.0, 0.5), (1.0, 0.5), plot_args={'linewidth':1,'color':'black'})
#sp.save(plot_dir + "test_array" + ".png")

