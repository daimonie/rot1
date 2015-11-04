#This file will test a intensity-coloured 2d plot with a contour to show.
#I want that sort of plot for my presentation :)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

r, phi = np.meshgrid( np.linspace(0, 1, 1000), np.linspace(-2*np.pi, 2*np.pi, 1000))

x = r * np.cos(phi)
y = r * np.sin(phi)

z = np.cos( r * np.cos(phi)) - np.cos( r * np.sin(phi))
#from http://matplotlib.org/examples/images_contours_and_fields/pcolormesh_levels.html

# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array. 
levels = MaxNLocator(nbins=20).tick_values(z.min(), z.max())


# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('summer')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig, (ax0, ax1) = plt.subplots(nrows=2)

im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
fig.colorbar(im, ax=ax0)
ax0.set_title('pcolormesh with levels')


# contours are *point* based plots, so convert our bound into point
# centers
cf = ax1.contourf(x, y, z, levels=levels,
                  cmap=cmap)
fig.colorbar(cf, ax=ax1)
ax1.set_title('contourf with levels')

# adjust spacing between subplots so `ax1` title and `ax0` tick labels
# don't overlap
fig.tight_layout()

plt.show()