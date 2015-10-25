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

#Commandline arguments instruction.
parser	= argparse.ArgumentParser(prog="Josephson.Py",
  description = "This file calculates the josephson current, given a microscopic mechanis, scattering mechanism and a resolution. Note that this file also has inspection modes.")  
parser.add_argument('-n', '--number', help='Resolution onf the calculation.', default = 10, action='store', type = int)  
parser.add_argument('-d', '--gapfunction', help='Resolution onf the calculation.', default = 10, action='store', type = int)  
parser.add_argument('-m', '--mechanism', help='Resolution onf the calculation.', default = 10, action='store', type = int)   
parser.add_argument('-p', '--plot', help='Resolution onf the calculation.', default = 10, action='store', type = int)  
parser.add_argument('-f', '--filename', help='Sets the filename. Only saves when given.', default = "default.png", action='store', type = str) 
args	= parser.parse_args() 


N		= args.number 		#Resolution
gapfunction	= args.gapfunction 	#Gap function
mechanism	= args.mechanism 	#Scattering
plotMode	= args.plot 		#Plot Mode
filename	= args.filename 	#Filename if we want to save

#Some physical constants.
hbar	=  physical_constants["Planck constant over 2 pi"][0]
mass	=  physical_constants["electron mass"][0]
ev	=  physical_constants["electron volt"][0]
lattice	= 3.787 * physical_constants["Angstrom star"][0] #(From Wikipedia)


eta	= hbar**2 / 2 / mass / ev / lattice**2
deltaS	= 3e-3
kFermi	= 3*np.pi/4
fluxnum = 2.0 
#small changes

dFlux = 4.0*fluxnum*np.pi/N
dDelta = 10*deltaS/N
dK = kFermi/N
dPhi = 2*np.pi/N

dkdphi = dK*dPhi;
#Make our arrays of parameters.
fluxArray	= np.arange(-fluxnum*2*np.pi,fluxnum*2*np.pi, dFlux)
deltaArray	= np.arange(2*deltaS/N, 10*deltaS, dDelta)
kArray		= np.arange(kFermi/N, kFermi, dK)
phiArray	= np.arange(0,2*np.pi, dPhi) 
#We define the figure here so that the different modes can assign labels/titles
fig = plt.figure(figsize=(15,15))
#Define Lambda functions.
#Gap Function.
if gapfunction == 0: #d_{x^2-y^2}
	Delta	= lambda dd, kk, pp: dd * ( np.cos(kk*np.cos(pp)) - np.cos(kk*np.sin(pp)))
elif gapfunction == 1: #d_{xy}
	Delta = lambda dd, kk, pp: dd*np.sin(kk*np.cos(pp))*np.sin(kk*np.sin(pp))
else:
	raise Exception("Unknown gap function.") 
#The tunnel 'function' is required for the other functions.
Heaviside = lambda xx: 1.0 * (xx>0) 
if mechanism == 0: #constant rate
	tunnel = lambda kk, pp: 1.0 
elif mechanism == 1: #A slice of k-space.
	tunnel = lambda kk, pp: Heaviside(pp)*Heaviside(pp - np.pi/2)*Heaviside(kk-kFermi/2)
else:
	raise Exception("Unknown tunnel function.") 
#Energies.
Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
#Current
dCurrent	= lambda ff, dd, kk, pp:tunnel(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(np.angle(Delta(dd,kk,pp)) + ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
#Pre-plotting
ax = fig.add_subplot(111, projection='3d')  
ax.view_init(50, 80) 
#Mode switching
if plotMode == 0: #Plot tunnel function in k-space
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	z = tunnel(x, y)
	
	plt.xlabel("k_x")
	plt.ylabel("k_y")
	plt.title("Inspecting tunneling matrix, N=%d, d=%d, m=%d, p=%d" % (N,gapfunction, mechanism, plotMode))
elif plotMode == 1: 
	flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
	z = dCurrent(flux,delta, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, deltaArray)
	
	zero = np.max(z[np.abs(x)<0.5*dFlux])/np.max(z)
	
	plt.xlabel("k_x")
	plt.ylabel("k_y")
	plt.title("Inspecting Josephson current,zero-flux/max ratio %2.3e, N=%d, d=%d, m=%d, p=%d" % (zero, N,gapfunction, mechanism, plotMode))
elif plotMode == 2: 
	k, phi = np.meshgrid(kArray,phiArray);
	z = dCurrent(np.pi/2,0.3*deltaS, k, phi)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	plt.xlabel("k_x")
	plt.ylabel("k_y")
	plt.title("Inspecting Josephson current element, N=%d, d=%d, m=%d, p=%d" % (N,gapfunction, mechanism, plotMode))
else:
	raise Exception("Unknown plot mode.");  
ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, antialiased=False)

if filename != "default.png":	
	fig.savefig(filename)
else:
	plt.show()