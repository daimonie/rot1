import numpy as np;
from mpl_toolkits.mplot3d import axes3d;
import matplotlib.pyplot as plt;
from matplotlib import cm

#equation (1), harlinger
#>>> nx, ny = (3, 2)
#>>> x = np.linspace(0, 1, nx)
#>>> y = np.linspace(0, 1, ny)
#>>> xv, yv = meshgrid(x, y)

numberOfCells = 1;
resolution = numberOfCells*100;

k_x, k_y = np.meshgrid(
   np.linspace(-numberOfCells, numberOfCells, resolution),
   np.linspace(-numberOfCells, numberOfCells, resolution)
);
relativeDelta = (np.cos(k_x*2*np.pi) - np.cos(k_y*2*np.pi))**4+1;

fig = plt.figure();
ax = fig.add_subplot(111,projection='3d');

offset = numberOfCells+2;
limit  = 20;

ax.plot_surface(k_x, k_y, relativeDelta, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(k_x, k_y, relativeDelta, zdir='z', offset=-offset, cmap=cm.coolwarm)
cset = ax.contour(k_x, k_y, relativeDelta, zdir='x', offset=-offset, cmap=cm.coolwarm)
cset = ax.contour(k_x, k_y, relativeDelta, zdir='y', offset=offset, cmap=cm.coolwarm)
 
ax.set_xlabel('k_x')
ax.set_xlim(-(numberOfCells+1) , (numberOfCells+1))
ax.set_ylabel('k_y')
ax.set_ylim(-(numberOfCells+1) , (numberOfCells+1) )
ax.set_zlabel('Delta/Delta_0')
ax.set_zlim(-offset , limit)

plt.show()
