from pylab import genfromtxt
from matplotlib import pyplot as plt

energy_file = "/disks/ceres/makemake/acomp/tr/phantomruns3/adam_trinary/energy.ev"
outfilename = "energy_components_plot_phantom"
plot_loc = "/disks/ceres/makemake/acomp/abatten/masters/jstaff/sim2/plots/"

data = genfromtxt(energy_file, skip_header=1)

yscale = 1e46

Eunits = 3.793e48 / yscale
tunits = 1.594e3 / (365.25 * 24 * 60 * 60)

ax = plt.plot(111)

plt.plot(data[:, 0] * tunits, data[:, 1] * Eunits, "k", linewidth=3.5, label = r"$\mathrm{E_{Total}}$")
plt.plot(data[:, 0] * tunits, data[:, 3] * Eunits, "#377eb8", label = r"$\mathrm{K_{Total}}$")
plt.plot(data[:, 0] * tunits, data[:, 2] * Eunits, "#4daf4a", label = r"$\phi_\mathrm{{Total}}$")
plt.plot(data[:, 0] * tunits, data[:, 4] * Eunits, "#e41a1c", label = r"$\mathrm{U_{Total}}$")
plt.plot(data[:, 0] * tunits, (data[:, 11] + data[:, 12]) * Eunits, "#f781bf", label = r"$\mathrm{K_{Gas}}$")
plt.plot(data[:, 0] * tunits, data[:, 6] * Eunits, "#984ea3", label = r"$\mathrm{K_{Particles}}$")
plt.plot(data[:, 0] * tunits, data[:, 15] * Eunits, "#CCCC00", label = r"$\phi_\mathrm{{GG}}$")
plt.plot(data[:, 0] * tunits, data[:, 5] * Eunits, "#ff7f00", label = r"$\phi_\mathrm{{PP}}$")
plt.plot(data[:, 0] * tunits, data[:, 16] * Eunits, "#a65628", label = r"$\phi_{\mathrm{PG}}$")


plt.ylim(-3.5,1.5)
plt.xlim(0,11.5)

plt.xlabel(r"$\mathrm{Time}\ (\mathrm{yr})$", fontsize=20)
plt.ylabel(r"$\mathrm{Energy}\ (10^{46}\ \mathrm{ergs})$", fontsize=20)
plt.tight_layout()
plt.legend(ncol=3, loc='lower left', frameon=False)

plt.savefig(plot_loc + outfilename + ".png")
plt.savefig(plot_loc + outfilename + ".pdf")
plt.savefig(plot_loc + outfilename + ".eps")

plt.show()




