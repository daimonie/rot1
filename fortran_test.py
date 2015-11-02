#Regular calculation modules.
import numpy as np 
import scipy as sp 
#Allows a debug-output stream.
import sys as sys 
#Physical constants list.
from scipy.constants import *
#Time differences.
import time as time  
#Command line arguments.
import argparse as argparse 
#Plotting.
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
from laosto import *
from conductivity import * 

N		= 15
band		= 0
#Some physical constants.
hbar	=  physical_constants["Planck constant over 2 pi"][0]
mass	=  physical_constants["electron mass"][0]
ev	=  physical_constants["electron volt"][0]
lattice	= 3.787 * physical_constants["Angstrom star"][0] #(From Wikipedia)

eta	= hbar**2 / 2 / mass / ev / lattice**2
deltaS	= 3e-3
fluxnum = 2.0 
#small changes
dFlux = 4.0*fluxnum*np.pi/N
dDelta = 10*deltaS/N
dPhi = 2*np.pi/N

#Make our arrays of parameters.
fluxArray	= np.arange(-fluxnum*2*np.pi,	fluxnum*2*np.pi+dFlux,	dFlux)
deltaArray	= np.arange(2*deltaS/N, 	10*deltaS+dDelta, 	dDelta)
phiArray	= np.arange(0,			2*np.pi+dPhi, 		dPhi) 


#Calculate fermi surface

kF0 = [0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi, -1, -1]
system = LAOSTO(mu=0, H=30, theta=np.pi/4, g=5, gL=1, tl1=340, tl2=340, th=12.5, td=12.5, dE=60, dZ=15, dSO=5)


kFermis, indices = kF_angle(system, phiArray, kF0)
start = indices[band, 1]
end = indices[band, 2]

fermiLevel = kFermis[start:end,:]
fermiSurface = ((fermiLevel**2).sum(axis=1))**0.5

#Parameters, small changes, arrays; these are here because they depend on the fermi surface calculation, which itself depends on phiArray
kFermi		= np.max(fermiSurface)
dK		= kFermi/N
kArray		= np.arange(1e-15, kFermi+2*dK, dK) 
dkdphi		= dK*dPhi;



flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);



startTime = time.time();
 
#time to introduce the fortran function and get a move on



#end of test
print "Time elapsed is %2.3f [s]." % (time.time() - startTime);