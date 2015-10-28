import numpy as np;
#Plotting.
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

N = 10;
x = np.arange(0,100, 100/N)
y = np.arange(0,100, 100/N)

z = np.sin(y);

f = lambda xx, yy: ( (xx-z[x=xx])**2 + yy**2)**0.5

[xgrid, ygrid] = np.meshgrid(x,y);
values = f(xgrid,ygrid)
 
ax = fig.add_subplot(111, projection='3d')  
ax.view_init(50, 80) 
plt.xlabel("x")
plt.ylabel("y") 
plt.title("f(x,y) using z")
ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, antialiased=False) 
plt.show();