#This file will test a intensity-coloured 2d plot with a contour to show.
#I want that sort of plot for my presentation :)
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator

r, phi = np.meshgrid( np.linspace(0, 1, 1000), np.linspace(-2*np.pi, 2*np.pi, 1000))

x = r * np.cos(phi)
y = r * np.sin(phi)

z = r * np.sin(phi)



# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('summer')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
im = ax0.pcolormesh(x, y, z, cmap=cmap, norm=norm)
fig.colorbar(im, ax=ax0)