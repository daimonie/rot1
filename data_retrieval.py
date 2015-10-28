import numpy as np;
#Plotting.
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
 
generate_data = lambda pp: np.abs(np.sin(pp)) + 1.0;

N = 100;
####################################################################################################################### 
dr = 2.0/N;
dphi = 2*np.pi/N;
radiusArray = np.arange(0, 2.0+dr, dr)
phiArray = np.arange(0,2*np.pi+dphi, dphi)

dataArray = generate_data(phiArray); #z is a placeholder for an array of data. The number of data points is appropriate for y.shape

radius, phi = np.meshgrid(radiusArray, phiArray);

# pp is for phi 
# rr is for radius

retrieve = lambda pp: map(lambda pp: dataArray[pp==phiArray][0], phi)
f = lambda rr, pp: (1.0 + rr * (np.cos(pp)**2 + np.sin(pp)**2))*retrieve(phi)


x = radius * np.cos(phi);
y = radius * np.sin(phi);

surface = f(x, y)
#######################################################################################################################
#plotting
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(111, projection='3d')  
ax.view_init(50, 80) 
plt.xlabel("x")
plt.ylabel("y") 
plt.title("f(x,y) using z")
ax.plot_surface(x,y, surface, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, antialiased=False) 
plt.show();