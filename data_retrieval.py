import numpy as np;
#Plotting.
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt

Heaviside =  lambda xx: 1.0 * (xx>0)
generate_data = lambda pp: 0.24*(np.abs(np.sin(4*pp))+pp)

N = 40;
####################################################################################################################### 
dr = 2.0/N;
dphi = 2*np.pi/N;
radiusArray = np.arange(0, 2.0+dr, dr)
phiArray = np.arange(0,2*np.pi+dphi, dphi)

dataArray = generate_data(phiArray); #z is a placeholder for an array of data. The number of data points is appropriate for y.shape

radius, phi = np.meshgrid(radiusArray, phiArray);

#umptieth attempt at fixing this
def populate(array):
	population = phi*0.
	for i in range(N):
		population[:, i] = array;
	return population
# pp is for phi 
# rr is for radius
retrieve = lambda pp: populate(dataArray)
f = lambda rr, pp: 1.+3.*(1+rr)**(-2)
f3 = lambda rr, pp: Heaviside( retrieve(pp)- rr )
f2 = lambda rr, pp: f(rr,pp) * f3(rr,pp)


x = radius * np.cos(phi)
y = radius * np.sin(phi)

surface1 = f(radius, phi)
surface2 = f2(radius, phi)
surface3 = f3(radius, phi)
 
#######################################################################################################################
#plotting
fig = plt.figure(figsize=(15,15))
ax = fig.add_subplot(311, projection='3d')  
ax.view_init(50, 90) 
plt.xlabel("x")
plt.ylabel("y") 
plt.title("f(x,y)")
ax.plot_surface(x,y, surface1, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, antialiased=False) 

ax = fig.add_subplot(312, projection='3d')  
ax.view_init(60, 160) 
plt.xlabel("x")
plt.ylabel("y") 
plt.title("f(x,y) using z")
ax.plot_surface(x,y, surface2, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, antialiased=False)  


ax = fig.add_subplot(313, projection='3d')  
ax.view_init(60, 160) 
plt.xlabel("r")
plt.ylabel("phi") 
plt.title("Heaviside(retrieve - rr)")
ax.plot_surface(radius, phi, surface3, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, antialiased=False) 
plt.show();