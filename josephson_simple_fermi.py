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
parser.add_argument('-t', '--tunnel', help='Sets the tunneling model.', default = 0, action='store', type = int);
parser.add_argument('-a', '--parama', help='Sets a parameter for some function (debug).', default = 1.0, action='store', type = float);
parser.add_argument('-b', '--paramb', help='Sets a parameter for some function (debug).', default = 1.0, action='store', type = float);
parser.add_argument('-f', '--filename', help='Sets the filename. Only saves when given.', default = "default.png", action='store', type = str);
args	= parser.parse_args();


N	= args.number;
switch	= args.switch;
tunnel	= args.tunnel;
filename= args.filename;
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
if N > 100:
	raise Exception("You know this is going to crash...");
flux	= np.arange(-fluxnum*2*np.pi,fluxnum*2*np.pi, 2.0*fluxnum*np.pi/N)
delta	= np.arange(2*deltaS/N, 10*deltaS, 2*deltaS/N)
k	= np.arange(kFermi/N,kFermi, kFermi/N)
phi	= np.arange(0,2*np.pi, 2*np.pi/N)

dkdphi	= kFermi * 2.0*fluxnum*np.pi/ N**2

Heaviside = lambda xx: 1 * (xx>0)


fig = plt.figure(figsize=(15,15))
	
	
if tunnel == 0:#Constant rate; no tunneling effects
	tunnelMatrix = lambda kk, pp: 1;
elif tunnel==1:#
	tunnelMatrix = lambda kk, pp:  Heaviside(np.pi*0.85-pp)*Heaviside(pp - 0.35*np.pi)*Heaviside(kk-kFermi/2.0)
if switch == 0:

	flux, delta, k, phi = np.meshgrid(flux,delta,k,phi);
	# lambda functions
	Delta		= lambda dd, kk, pp: dd * ( np.cos(kk*np.cos(pp)) - np.cos(kk*np.sin(pp)))
	Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
	EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
	currentKPhi	= lambda ff, dd, kk, pp:tunnelMatrix(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
	
	print "You have chosen the current as a function of applied flux and Delta-parameter of x2-y2, with tunneling model %d." % tunnel
	current = currentKPhi(flux,delta, k, phi).sum(axis=-1).sum(axis=-1)

	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
 
	ax = fig.add_subplot(111, projection='3d') 
	ax.plot_surface(flux[..., 0, 0], delta[..., 0, 0], current, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False)
	
	
	
	ax.view_init(10, 80)
	plt.xlabel("Flux")
	plt.ylabel("Delta_0")
	plt.title("CPR for JJ s-wave into d-wave x^2-y^2, tunnel mode %d" % tunnel) 
elif switch==1:
	k, phi = np.meshgrid(k, phi)
	
	Delta		= lambda dd, kk, pp: dd * ( np.cos(kk*np.cos(pp)) - np.cos(kk*np.sin(pp)))
	Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
	EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
	currentKPhi	= lambda ff, dd, kk, pp: tunnelMatrix(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
	
	print "[x2-y2] You have chosen the current in k-space, where Flux is set to pi*%2.3f and Delta-parameter is set to DeltaS*%2.3f, with tunneling model %d." % (parama,paramb, tunnel)
	#currentKPhi plot. Is apparently needed :-)  
	current = currentKPhi(np.pi*parama,deltaS*paramb,k,phi);
 
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
 
	ax = fig.add_subplot(111, projection='3d') 
	
	
	x = k*np.cos(phi);
	y = k*np.sin(phi);
	
	ax.plot_surface(x, y, current, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False) 
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title("CPR term under sumfor JJ s-wave into d-wave x^2-y^2, tunnel mode %d" % tunnel) 
elif switch == 2: #Let's do d_{xy}!

	flux, delta, k, phi = np.meshgrid(flux,delta,k,phi)
	# lambda functions
	Delta		= lambda dd, kk, pp: dd * ( np.sin(kk*np.cos(pp))*np.sin(kk*np.sin(pp)))
	Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
	EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
	currentKPhi	= lambda ff, dd, kk, pp: tunnelMatrix(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
	
	print "You have chosen the current as a function of applied flux and Delta-parameter of xy, with tunneling model %d." % tunnel
	current = currentKPhi(flux,delta, k, phi).sum(axis=-1).sum(axis=-1)

	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
 
	ax = fig.add_subplot(111, projection='3d') 
	ax.plot_surface(flux[..., 0, 0], delta[..., 0, 0], current, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False)
	
	
	ax.view_init(10, 80)
	plt.xlabel("Flux")
	plt.ylabel("Delta_0")
	plt.title("CPR for JJ s-wave into d-wave xy, tunnel mode %d" % tunnel) 
elif switch==3:

	k, phi = np.meshgrid(k,phi)
	Delta		= lambda dd, kk, pp: dd * ( np.sin(kk*np.cos(pp))*np.sin(kk*np.sin(pp)))
	Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
	EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
	currentKPhi	= lambda ff, dd, kk, pp: tunnelMatrix(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(ff) /( (Energy(kk,dd)+EnergyR(kk,dd,pp)) * (Energy(kk,dd)*EnergyR(kk,dd,pp)))
	
	print "[xy] You have chosen the current in k-space, where Flux is set to pi*%2.3f and Delta-parameter is set to DeltaS*%2.3f, with tunneling mode %d." % (parama,paramb, tunnel)
	#currentKPhi plot. Is apparently needed :-)  
	current = currentKPhi(np.pi*parama,deltaS*paramb,k,phi);
 
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt
 
	ax = fig.add_subplot(111, projection='3d') 
	
	x = k*np.cos(phi);
	y = k*np.sin(phi);
	
	ax.plot_surface(x, y, current, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False) 
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title("CPR term under sumfor JJ s-wave into d-wave xy, tunnel mode %d" % tunnel)
elif switch == 4: #Tunnel matrix
	k, phi = np.meshgrid(k,phi)
	tunnelfdkp = lambda kk, pp: tunnelMatrix(kk,pp) 
	
	tunnelFunc = tunnelfdkp(k,phi)
	
	print "Inspecting tunnelingMatrix %d." % tunnel 

	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.pyplot as plt

	ax = fig.add_subplot(111, projection='3d')  
	ax.view_init(50, 80)
	
	x = k*np.cos(phi);
	y = k*np.sin(phi);
	
	ax.plot_surface(x, y, tunnelFunc, rstride=int(parama), cstride=int(paramb), cmap=cm.summer,linewidth=0, antialiased=False)
	plt.xlabel("x")
	plt.ylabel("y")
	plt.title("Inspecting tunnelingMatrix %d" % tunnel) 
	
	
if filename != "default.png":	
	fig.savefig(filename)
else:
	plt.show()