from pylab import genfromtxt
from matplotlib import pyplot as plt

sep_file = "/disks/ceres/makemake/acomp/tr/phantomruns3/adam_trinary/separation_vs_time_trinary.ev"
outfilename = "separations_phantom"
plot_loc = "/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/"

data = genfromtxt(sep_file, skip_header=1)

yscale = 1

lunits = 1#6.960e10 / yscale
tunits = 1.594e3 / (365.25 * 24 * 60 * 60)

ax = plt.plot(111)

plt.plot(data[:, 0] * tunits, data[:, 4] * lunits, "b", linewidth=2, label = r"$\mathrm{Inner\ Planet}}$")
plt.plot(data[:, 0] * tunits, data[:, 8] * lunits, "r", linewidth=2, label = r"$\mathrm{Outer\ Planet}$")


plt.ylim(0,4.5)
plt.xlim(0,11.5)

plt.xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
plt.ylabel(r"$\mathrm{Separation}\ (\mathrm{R_\odot})$", fontsize=20)
plt.tight_layout()
plt.legend(ncol=1, loc='lower left', frameon=False)

plt.savefig(plot_loc + outfilename + ".png")
plt.savefig(plot_loc + outfilename + ".pdf")
plt.savefig(plot_loc + outfilename + ".eps")

plt.show()


