from yt.mods import *
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

#plot_dir = '/disks/ceres/makemake/acomp/abatten/tauris/plots/'
plot_dir = '/disks/ceres/makemake/acomp/abatten/jstaff/sim2/plots/'
#root_dir = '/disks/ceres/makemake/acomp/riaconi/rxi552/Tauris_vandenHeuvel/128_4lev_25AUbox_273Rsun/run/'

root_dir = '/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj/'

directory_list = ['DD0000/CE0008', 'DD0005/CE0005','DD0015/CE0015', 'DD0020/CE0020', 'DD0025/CE0025', 'DD0030/CE0030', 'DD0035/CE0035', 'DD0060/CE0060','DD0140/CE0140' ]

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
    sp = SlicePlot(ds, 'z', "ThermalEnergy", width = 1.0)
    #sp.annotate_velocity(factor=10, normalize=False) #overplots velocity field vectors
    sp.annotate_particles(1.0, p_size=50.0, marker='o', col='black') #overplots the projection on the axis of the particLes
    sp.annotate_text((0.1,0.88), "t = " + str("{0:.2f}".format(ds.current_time)) + "yr")

    # Ensure the colorbar limits match for all plots
    sp.set_zlim('ThermalEnergy', 7e11, 3e15)

    # This forces the slice to redraw itself on the AxesGrid axes.
    plot = sp.plots['ThermalEnergy']
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    sp._setup_plots()

plt.savefig(plot_dir+'multiplot_3x3_z_thermal_array_jstaff.png')
plt.show()

