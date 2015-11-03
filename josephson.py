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
#Commandline arguments instruction.
parser	= argparse.ArgumentParser(prog="Josephson.Py",
  description = "This file calculates the josephson current, given a microscopic mechanis, scattering mechanism and a resolution. Note that this file also has inspection modes.")  
parser.add_argument('-n', '--number', help='Resolution of the calculation.', default = 10, action='store', type = int)  
parser.add_argument('-d', '--gapfunction', help='Gap function.', default = 10, action='store', type = int)  
parser.add_argument('-m', '--mechanism', help='Scattering mechanism.', default = 10, action='store', type = int)   
parser.add_argument('-p', '--plot', help='Plotting mode. 0 tunnel matrix, 1 josephson-surface, 2 josephson-current, 3 gap function, 4 fermi surface, 5 tunnel matrix times gap times fermi.', default = 10, action='store', type = int)  
parser.add_argument('-k', '--fermi', help='Fermi Surface.', default = 10, action='store', type = int)  
parser.add_argument('-b', '--band', help='Fermi Surface band.', default = 10, action='store', type = int)  
parser.add_argument('-f', '--filename', help='Sets the filename. Only saves when given.', default = "default.png", action='store', type = str) 
parser.add_argument('-q', '--fortran', help='Use fortran for calculation of fermi surface.', default = True, action='store', type = bool) 
parser.add_argument('-s', '--silent', help='Do not plot, do not save', default = False, action='store', type = bool) 
args	= parser.parse_args() 


N		= args.number 		#Resolution
gapfunction	= args.gapfunction 	#Gap function
mechanism	= args.mechanism 	#Scattering
plotMode	= args.plot 		#Plot Mode
filename	= args.filename 	#Filename if we want to save
fermi		= args.fermi 		#Type of Fermi surface
band		= args.band	 	#Band number. There are typically 4 bands as far as I can see.
useFortran	= args.fortran		#Silent. Don't plot, don't save.
silent		= args.silent	 	#Silent. Don't plot, don't save.
 
startTime = time.time();

if filename != "default.png":
	print "Saving Figure (%s) for -d %d -m %d -p %d -n %d -k %d -b %d." % (filename, gapfunction, mechanism, plotMode, N, fermi, band)
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
kFermi	= np.max(fermiSurface)
dK = kFermi/N
kArray		= np.arange(1e-15, kFermi+2*dK, dK) 
dkdphi = dK*dPhi;


#We define the figure here so that the different modes can assign labels/titles
fig = plt.figure(figsize=(15,15))
#Define Lambda functions.
Heaviside 	= lambda xx: 1.0 * (xx>0) 
Dirac 		= lambda xx: 1.0 * (xx==0) 
#Gap Function.
if gapfunction == 0: #d_{x^2-y^2}
	Delta	= lambda dd, kk, pp: dd * ( np.cos(kk*np.cos(pp)) - np.cos(kk*np.sin(pp)))
elif gapfunction == 1: #d_{xy}
	Delta = lambda dd, kk, pp: dd*np.sin(kk*np.cos(pp))*np.sin(kk*np.sin(pp))
elif gapfunction == 2: #s-wave
	Delta = lambda dd, kk, pp: dd
else:
	raise Exception("Unknown gap function.") 
#By far, the easiest way to incorporate the Fermi Surface is by just letting it tag along
#	with the tunnel functionality.
#	the only requirement is that the Fermi surface is at least as high as the maximum calculated length.
#	This requires us to calculate the fermi surface and, sadly, remake the kFermi, dK and kArray variables. 
#	It is, however, the smallest step.
if band >= 4:
	raise Exception("There are only four bands");

if fermi == 0:
	Fermi = lambda kk, pp: np.max(fermiSurface);
elif fermi == 1:   
	if useFortran: 
		import populate as fpopulate  
		def Fermi (kk, pp): 
			if len(kk.shape) == 4:
				return fpopulate.populate.fermi_integrand(first=kk.shape[0], second=kk.shape[1], third=kk.shape[2], fourth=kk.shape[3],
					fermi_surface=fermiSurface, angle_array=phiArray, angle=pp, radius=kk);
			elif len(kk.shape) == 2: 
				return fpopulate.populate.fermi_contour(first=kk.shape[0], second=kk.shape[1],
					fermi_surface=fermiSurface, angle_array=phiArray, angle=pp, radius=kk);
	else:
		def Fermi( kk, pp):	
			raise Exception("This function is being made into a fortran function");
			# I am sorry to say that this is the best implementation I could find. My python needs some work.
			# Anyway, if have a new mode with a new meshgrid, you'll need to add a clause.
			fermi = pp; 
			
			popTime = time.time(); 
			if len(fermi.shape) == 4: 
				for i in range(0,fermi.shape[0]):
					for j in range(0,fermi.shape[1]):
						for ii in range(0,fermi.shape[2]):
							for jj in range(0,fermi.shape[3]):
								#fermi[i, j, ii, jj] = 1;
								if (fermiSurface[pp[i,j, ii, jj]==phiArray][0]  > kk[i,j, ii, jj]):
									fermi[i,j, ii, jj] = 1.  
								else:
									fermi[i,j, ii, jj] = 0. 
			elif len(fermi.shape) == 2:
				#print "Using Fermi (N+1,N+2)";
				#Principle has been tested
				for i in range(0, fermi.shape[0]): 
					for j in range(0,fermi.shape[1]): 
						#print "(%d,%d): %s %s" % (i, j, fermiSurface[pp[i,j]==phiArray][0], kk[i,j])
						if (fermiSurface[pp[i,j]==phiArray][0]  > kk[i,j]):
							fermi[i,j] = 1.  
						else:
							fermi[i,j] = 0.
			else: 
				raise Exception("Tuple too large.");
			
			print "Time taken for population is [%2.3f] s." % (time.time() - popTime);
			#print fermi
			return fermi;
else:
	raise Exception("Unknown fermi surface requested.")

#The tunnel 'function' is required for the other functions.
if mechanism == 0: #constant rate 
	tunnel = lambda kk, pp: 1.0 * Fermi(kk, pp)
elif mechanism == 1: #A slice of k-space.
	tunnel = lambda kk, pp:  Heaviside(np.pi/4.-pp) * Fermi(kk, pp)
else:
	raise Exception("Unknown tunnel function.") 
#Energies.
Energy		= lambda kk, dd: (eta*kk**2 + dd**2)**0.5
EnergyR		= lambda kk, dd, pp: (eta*kk**2 + Delta(kk,dd,pp)**2)**0.5
#Current
dCurrent	= lambda ff, dd, kk, pp:tunnel(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(np.angle(Delta(dd,kk,pp)) + ff) /( (Energy(kk,deltaS)+EnergyR(kk,dd,pp)) * (Energy(kk,deltaS)*EnergyR(kk,dd,pp)))
#Pre-plotting
ax = fig.add_subplot(111, projection='3d')  
ax.view_init(50, 80) 
#Mode switching 
title = "Unknown Mode." 
if plotMode == 0: #Plot tunnel function in k-space
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	 
	z = tunnel(k, phi)
	
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Tunneling Matrix";
elif plotMode == 1:  
	
	flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
	 
	z = dCurrent(flux,delta, k, phi).sum(axis=-1).sum(axis=-1) 
	 
	
	x,y = np.meshgrid(fluxArray, deltaArray)
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\Delta^0_m$")   
	title = "Current";
elif plotMode == 2:  
	ax.view_init(0, 90) 
	#It's just mode 1 with a different view angle
	
	flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
	z = dCurrent(flux,delta, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, deltaArray) 
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\Delta^0_m$")   
	title = "Current_k_phi";
elif plotMode == 3:
	
	ax.view_init(30, 30) 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	z = Delta(1, k, phi);
	
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Gap function";
elif plotMode == 4:
	ax.view_init(30, 30) 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	z =  Fermi(k, phi);
	
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Fermi disc";
elif plotMode == 5:
	
	ax.view_init(30, 30) 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	deltaTunnel = lambda kk, pp: tunnel(kk,pp) * Delta(deltaS, kk, pp)
	
	z = deltaTunnel(k, phi);
	
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Gap function times tunnel"; 
else:
	raise Exception("Unknown plot mode.");   

stride = int(N/100);
plt.title("%s, N=%d, d=%d, m=%d, p=%d, k=%d, b=%d" % (title, N,gapfunction, mechanism, plotMode, fermi, band))
ax.set_aspect('equal')
ax.plot_surface(x, y, z, rstride=stride, cstride=stride, cmap=cm.summer,linewidth=0.1)  


print "Elapsed time %2.3f" % (time.time() - startTime)
if silent:
	print "Silent Mode; no plotting, no saving.";
	print "Title [%s]." % title;
	plt.close();
if filename != "default.png":	
	fig.savefig(filename)
else:
	plt.show()