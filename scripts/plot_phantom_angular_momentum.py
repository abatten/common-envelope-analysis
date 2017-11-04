from pylab import genfromtxt
from matplotlib import pyplot as plt

ang_mom_file = "/disks/ceres/makemake/acomp/tr/phantomruns3/adam_trinary/energy.ev"
outfilename = "angular_momentum_plot_phantom"
plot_loc = "/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/"

data = genfromtxt(ang_mom_file, skip_header=1)

yscale = 1e29

AMunits = 1.2478e30 / yscale
tunits = 1.594e3 / (365.25 * 24 * 60 * 60)

ax = plt.plot(111)

plt.plot(data[:, 0] * tunits, data[:, 17] * AMunits, "k", linewidth=2, label = r"$\mathrm{J_{Total}}$")
plt.plot(data[:, 0] * tunits, data[:, 20] * AMunits, "#e7298a", linewidth=1.5, label = r"$\mathrm{J_{Particle}}$")
plt.plot(data[:, 0] * tunits, (data[:, 18]+data[:, 19]) * AMunits, "#117713", linewidth=1.5, label = r"$\mathrm{J_{Gas}}$")

plt.ylim(0,3.7)
plt.xlim(0,11.5)

plt.xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
plt.ylabel(r"$\mathrm{Angular\ Momentum}\ (10^{27}\ \mathrm{g\ cm^2s^{-1}})$", fontsize=20)
plt.tight_layout()
plt.legend(ncol=1, loc='lower left', frameon=False)

plt.savefig(plot_loc + outfilename + ".png")
plt.savefig(plot_loc + outfilename + ".pdf")
plt.savefig(plot_loc + outfilename + ".eps")

plt.show()

