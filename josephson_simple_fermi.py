import numpy as np;
import scipy as sp;
from scipy.constants import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt 

hbar =  physical_constants["Planck constant over 2 pi"][0]
mass =  physical_constants["electron mass"][0]
ev =  physical_constants["electron volt"][0]
lattice = 3.787 * physical_constants["Angstrom star"][0] #(From Wikipedia)


eta = hbar**2 / 2 / mass / ev / lattice**2;
deltas = 3e-3;

k_F = 3*np.pi/4;
#linspace is really slow. I'm serious. Let's not use it :-)
#arange is about seven times faster. See test.py
N = 100
deltam = np.arange(0, 2*deltas, 2*deltas/N)
flux = np.arange(-5, 5, 10.0/N)
angle = np.arange(0,2*np.pi,2*np.pi/N);

flux, deltam, angle = np.meshgrid(flux, deltam, angle)

current = np.multiply(np.multiply(np.cos(k_F*np.cos(angle)), deltam), np.sin(flux)) 

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(flux,deltam,current, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
 
plt.show()