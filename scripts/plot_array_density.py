from yt.mods import *
import numpy as np
import math
import os
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

#plot_dir = '/disks/ceres/makemake/acomp/abatten/tauris/plots/'
plot_dir = '/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/'
#root_dir = '/disks/ceres/makemake/acomp/riaconi/rxi552/Tauris_vandenHeuvel/128_4lev_25AUbox_273Rsun/run/'
#root_dir = '/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj/1part/'
root_dir = '/disks/ceres/makemake/acomp/jstaff/raijin/enzoamr-2.3/1Moz0.02late-2/128/4levels/2part/run210Mj/'

#directory_list = ['DD0170/CE0170','DD0190/CE0190', 'DD0200/CE0200', 'DD0210/CE0210', 'DD0220/CE0220', 'DD0230/CE0230', 'DD0240/CE0240', 'DD0250/CE0250','DD0260/CE0260' ]
directory_list = ['DD0000/CE0008','DD0010/CE0010', 'DD0020/CE0020', 
                  'DD0025/CE0025', 'DD0030/CE0030','DD0040/CE0040',
                  'DD0050/CE0050', 'DD0060/CE0060']

fig = plt.figure()

grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                nrows_ncols = (4, 2),
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
    pf = SlicePlot(ds, 'z', 'Density')
#    pf.annotate_velocity(factor=16, normalize=False)
    #pf.annotate_streamlines('particle_velocity_x', 'particle_velocity_y')
    pf.annotate_particles(1.0, p_size=25.0, marker='o', col='black')
    pf.annotate_text((0.1,0.88), "t = " + str("{0:.2f}".format(ds.current_time)) + "yr", text_args={'color':'white'})

    # Ensure the colorbar limits match for all plots
    pf.set_zlim('Density', 1e-12, 1e-4)

    # This forces the slice to redraw itself on the AxesGrid axes.
    plot = pf.plots['Density']
    plot.figure = fig
    plot.axes = grid[i].axes
    plot.cax = grid.cbar_axes[i]

    # Finally, this actually redraws the plot.
    pf._setup_plots()

plt.savefig('multiplot_4x2_z_density_array_jstaff.png', bbox_inches='tight',pad_inches=0)
plt.show()
