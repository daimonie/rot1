import numpy as np;
import scipy as sp;
import sys as sys;
from scipy.constants import *
import time as time;
#params are nice...
import argparse as argparse
#plotting
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
#script
#arguments
parser	= argparse.ArgumentParser(prog="python main.py",
  description = "Takes Delta(k) and plots phase/amp graphs. Use -help to see more information."); 
parser.add_argument('-n', '--number', help='Resolution onf the calculation.', default = 10, action='store', type = int);
parser.add_argument('-s', '--switch', help='Setting this to false shows dI(k,phi).', default = 0, action='store', type = int);
parser.add_argument('-a', '--parama', help='Sets a parameter for some function (debug).', default = 0.5, action='store', type = float);
parser.add_argument('-b', '--paramb', help='Sets a parameter for some function (debug).', default = 0.3, action='store', type = float);
args	= parser.parse_args();


N	= args.number;
switch	= args.switch;
parama	= args.parama;
paramb	= args.paramb;
#some constants
hbar	=  physical_constants["Planck constant over 2 pi"][0]
mass	=  physical_constants["electron mass"][0]
ev	=  physical_constants["electron volt"][0]
lattice	= 3.787 * physical_constants["Angstrom star"][0] #(From Wikipedia)


eta	= hbar**2 / 2 / mass / ev / lattice**2
deltaS	= 3e-3
kFermi	= 3*np.pi/4
fluxnum = 2.0;
#linspace is really slow. I'm serious. Let's not use it :-)
#arange is about seven times faster. See test.py

flux	= np.arange(-fluxnum*np.pi,fluxnum*np.pi, 2.0*fluxnum*np.pi/N)
delta	= np.arange(2*deltaS/N, deltaS, 2*deltaS/N)
k	= np.arange(kFermi/N,kFermi, kFermi/N)
phi	= np.arange(0,2*np.pi, 2*np.pi/N)

dkdphi	= kFermi * 2.0*fluxnum*np.pi/ N**2

flux, delta, k, phi = np.meshgrid(flux,delta,k,phi);
if switch == 0:
	# lambda functions
	Delta		= lambda dd, kk, pp: dd * ( np.cos(kk*np.cos(pp)) - np.cos(kk*np.sin(pp)))
	Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
	EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
	currentKPhi	= lambda ff, dd, kk, pp: dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
	
	print "You have chosen the current as a function of applied flux and Delta-parameter of x2-y2."
	current = currentKPhi(flux,delta, k, phi).sum(axis=-1).sum(axis=-1)

	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d') 
	ax.plot_surface(flux[..., 0, 0], delta[..., 0, 0], current, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False)
	plt.xlabel("Flux")
	plt.ylabel("Delta_0")
	plt.show()
	
elif switch == 1: #Let's do d_{xy}!
	# lambda functions
	Delta		= lambda dd, kk, pp: dd * ( np.sin(kk*np.cos(pp))*np.sin(kk*np.sin(pp)))
	Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
	EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
	currentKPhi	= lambda ff, dd, kk, pp: dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
	
	print "You have chosen the current as a function of applied flux and Delta-parameter of xy."
	current = currentKPhi(flux,delta, k, phi).sum(axis=-1).sum(axis=-1)

	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d') 
	ax.plot_surface(flux[..., 0, 0], delta[..., 0, 0], current, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False)
	plt.xlabel("Flux")
	plt.ylabel("Delta_0")
	plt.show()
else:
	print "You have chosen the current in k-space, where Flux is set to pi*%2.3f and Delta-parameter is set to DeltaS*%2.3f." % (parama,paramb)
	#currentKPhi plot. Is apparently needed :-)
	k	= np.arange(kFermi/N,kFermi, kFermi/N)
	phi	= np.arange(0,2*np.pi, 2*np.pi/N)
	k, phi = np.meshgrid(k,phi)
	current = currentKPhi(np.pi*parama,deltaS*paramb,k,phi);
 
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d') 
	ax.plot_surface(k, phi, current)
	plt.show()