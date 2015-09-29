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

epsilon = 0.3;

oscillation = np.cos(2*np.pi*k_x) - np.cos(2*np.pi*k_y);

relativeDeltaReal = (1-epsilon)*oscillation;
relativeDeltaImag = epsilon*2*np.sin(2*np.pi*k_x)*np.sin(2*np.pi*k_y);

fig = plt.figure();
ax = fig.add_subplot(121,projection='3d');

offset = numberOfCells+2;
limit  = 20;

ax.plot_surface(k_x, k_y, relativeDeltaReal, rstride=8, cstride=8, alpha=0.3)
cset = ax.contour(k_x, k_y, relativeDeltaReal, zdir='z', offset=-offset, cmap=cm.coolwarm)
cset = ax.contour(k_x, k_y, relativeDeltaReal, zdir='x', offset=-offset, cmap=cm.coolwarm)
cset = ax.contour(k_x, k_y, relativeDeltaReal, zdir='y', offset=offset, cmap=cm.coolwarm)
 
ax.set_xlabel('k_x')
ax.set_xlim(-(numberOfCells+1) , (numberOfCells+1))
ax.set_ylabel('k_y')
ax.set_ylim(-(numberOfCells+1) , (numberOfCells+1) )
ax.set_zlabel('Delta/Delta_0')
ax.set_zlim(-offset , limit)

plt.title("Real");


ax2 = fig.add_subplot(122,projection='3d');
 
ax2.plot_surface(k_x, k_y, relativeDeltaImag, rstride=8, cstride=8, alpha=0.3)
cset = ax2.contour(k_x, k_y, relativeDeltaImag, zdir='z', offset=-offset, cmap=cm.coolwarm)
cset = ax2.contour(k_x, k_y, relativeDeltaImag, zdir='x', offset=-offset, cmap=cm.coolwarm)
cset = ax2.contour(k_x, k_y, relativeDeltaImag, zdir='y', offset=offset, cmap=cm.coolwarm)
 
ax2.set_xlabel('k_x')
ax2.set_xlim(-(numberOfCells+1) , (numberOfCells+1))
ax2.set_ylabel('k_y')
ax2.set_ylim(-(numberOfCells+1) , (numberOfCells+1) )
ax2.set_zlabel('Delta/Delta_0')
ax2.set_zlim(-offset , limit)

plt.title("Imaginary");
plt.show()
