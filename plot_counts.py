# setup
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

binedges = np.array([15,25,50,100,350])
hdul = fits.open('~/data/mirror/swift/00096727251/bat/rate/sw00096727251brtms.lc')
plot_data = hdul[1].data
print(plot_data.columns)
print(plot_data)


# creating the figure
fig, axes = plt.subplots(nrows = 2, ncols = 1, sharex = True)

# plot 1: split by energy bands
for i in range(4):
    axes[0].plot(plot_data['TIME'], plot_data['COUNTS'][:, i], ".", label = f"{binedges[i]} - {binedges[i+1]} keV")
axes[0].legend()
axes[0].set(ylabel = "COUNTS", xlabel = "TIME (secs)")
y_points = np.array(plot_data['COUNTS'])

# # plot 2: counts vs. time for 15-50 keV energy bands
counts = np.sum(plot_data['COUNTS'][:,0:2], axis=-1)
axes[1].plot(plot_data['TIME'], counts, ".")
axes[1].set(ylabel = f"COUNTS ({binedges[0]} - {binedges[2]} keV)")

# plot 3: rate vs. time for 15-50 keV energy bands
#dt = np.median(np.diff(plot_data['TIME']))    # size of the timebin
#rate = np.sum(plot_data['COUNTS'][:,0:2], axis=-1) / dt
#axes[1].plot(plot_data['TIME'], rate, ".")
#axes[1].set(ylabel = f"RATE (counts {binedges[0]} - {binedges[2]} keV) / s")


plt.show()
