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
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
#LAO/STO model from Leiden
from laosto import *
from conductivity import *
#Commandline arguments instruction.
parser	= argparse.ArgumentParser(prog="Josephson.Py",
  description = "This file calculates the josephson current, given a microscopic mechanis, scattering mechanism and a resolution. Note that this file also has inspection modes.")  
parser.add_argument('-n', '--number', help='Resolution of the calculation.', default = 10, action='store', type = int)  
parser.add_argument('-d', '--gapfunction', help='Gap function.', default = 10, action='store', type = int)  
parser.add_argument('-m', '--mechanism', help='Scattering mechanism.', default = 10, action='store', type = int)   
parser.add_argument('-p', '--plot', help='Plotting mode.', default = 10, action='store', type = int)  
parser.add_argument('-k', '--fermi', help='Fermi Surface.', default = 10, action='store', type = int)  
parser.add_argument('-b', '--band', help='Fermi Surface band.', default = 10, action='store', type = int)  
parser.add_argument('-f', '--filename', help='Sets the filename. Only saves when given.', default = "default.png", action='store', type = str) 
parser.add_argument('-s', '--silent', help='Do not plot, do not save.', default = False, action='store', type = bool) 
parser.add_argument('-v', '--verbose', help='Output data for use with, say, gnuplot.', default = False, action='store', type = bool) 
parser.add_argument('--scatter', help='Use a scatterplot instead of a surface plot.', default = False, action='store', type = bool) 
parser.add_argument('--alpha', help='First view angle.', default = 1000, action='store', type = int)  
parser.add_argument('--beta', help='Second view angle.', default = 1000, action='store', type = int)  
parser.add_argument('--gamma', help='Plot transparency.', default = 0.9, action='store', type = float)  
parser.add_argument('--magnetic', help='Alter magnetic field.', default = 0.0, action='store', type = float)  
args	= parser.parse_args() 


N		= args.number
gapfunction	= args.gapfunction
mechanism	= args.mechanism
plotMode	= args.plot
filename	= args.filename
fermi		= args.fermi
band		= args.band
scatter		= args.scatter
alpha		= args.alpha
beta		= args.beta
gamma		= args.gamma
silent		= args.silent
magnetic	= args.magnetic
 
startTime = time.time();

print >> sys.stderr, "The plot will be viewed for angle (%2.3f, %2.3f) and with alpha set to %2.3f" % (alpha, beta, gamma)

if filename != "default.png":
	print >> sys.stderr, "Saving Figure (%s) for -d %d -m %d -p %d -n %d -k %d -b %d." % (filename, gapfunction, mechanism, plotMode, N, fermi, band)
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
phiArray	= np.linspace(0,		2*np.pi, N) 


#Calculate fermi surface

laoTime = time.time();
kF0 = [0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi, -1, -1]
system = LAOSTO(mu=0, H=magnetic, theta=np.pi/4, g=5, gL=1, tl1=340, tl2=340, th=12.5, td=12.5, dE=60, dZ=15, dSO=5)


kFermis, indices = kF_angle(system, phiArray, kF0)
start = indices[band, 1]
end = indices[band, 2]

fermiLevel = kFermis[start:end,:]
fermiSurface = ((fermiLevel**2).sum(axis=1))**0.5
print >> sys.stderr, "Magnetic Field H=%2.3f for the LAOSTO model." % magnetic
print >> sys.stderr, "Time taken for LAO/STO model is [%2.3f] s."% (time.time() - laoTime);
#Parameters, small changes, arrays; these are here because they depend on the fermi surface calculation, which itself depends on phiArray
kFermi	= np.max(fermiSurface)
dK 	= kFermi/N
kArray	= np.linspace(1e-15, kFermi, N) 
dkdphi  = dK*dPhi;

#Define Lambda functions.
Heaviside 	= lambda xx: 1.0 * (xx>0) 
Dirac 		= lambda xx: 1.0 * (xx==0) 
#Gap Function.
if gapfunction == 0: 
	print >> sys.stderr, "Selected d_{x^2-y^2} wave";
	Delta	= lambda dd, kk, pp: dd * ( np.cos(kk*np.cos(pp)) - np.cos(kk*np.sin(pp)))
elif gapfunction == 1:  
	print >> sys.stderr, "Selected d_{xy} wave";
	Delta = lambda dd, kk, pp: dd*np.sin(kk*np.cos(pp))*np.sin(kk*np.sin(pp))
elif gapfunction == 2:  
	print >> sys.stderr, "Selected s^0-wave (equiv. leading order)";
	Delta = lambda dd, kk, pp: dd + pp*0.
elif gapfunction == 3:  
	print >> sys.stderr, "Selected s^1-wave (equiv. next leading order)";
	Delta = lambda dd, kk, pp: deltaS + dd * ( np.cos(kk*np.cos(pp)) + np.cos(kk*np.sin(pp)))
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
	import populate as fpopulate  
	def Fermi (kk, pp):  
		popTime = time.time(); 
		if len(kk.shape) == 4:
			result = fpopulate.populate.fermi_integrand(first=kk.shape[0], second=kk.shape[1], third=kk.shape[2], fourth=kk.shape[3],
				fermi_surface=fermiSurface, angle_array=phiArray, angle=pp, radius=kk);
		elif len(kk.shape) == 2: 
			result =  fpopulate.populate.fermi_contour(first=kk.shape[0], second=kk.shape[1],
				fermi_surface=fermiSurface, angle_array=phiArray, angle=pp, radius=kk);
		print "Time taken for f.population is [%2.3f] s." % (time.time() - popTime);
		return result
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
EnergyPart 	= lambda kk,dd,pp: 1./( (Energy(kk,deltaS)+EnergyR(kk,dd,pp)) * (Energy(kk,deltaS)*EnergyR(kk,dd,pp)))
#Current
dCurrent	= lambda ff, dd, kk, pp:tunnel(kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(-np.angle(Delta(dd,kk,pp)) + ff) * EnergyPart(kk,dd,pp)
#Pre-plotting
flatWorld = True

if plotMode == 6:
	flatWorld = False


if flatWorld:
		
	#We define the figure here so that the different modes can assign labels/titles
	fig = plt.figure(figsize=(15,15))
	ax = fig.add_subplot(111, projection='3d')  
	ax.view_init(50, 80) 
else:
	fig, ax = plt.subplots()
#Mode switching 
title = "Unknown Mode." 
if plotMode == 0: #Plot tunnel function in k-space
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	 
	z = tunnel(k, phi)
	
	ax.set_aspect('equal'); #can be used for k-space, but this wants x-axis = y-axis
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Tunneling Matrix";
elif plotMode == 1:  
	
	flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
	 
	z = dCurrent(flux,delta, k, phi).sum(axis=-1).sum(axis=-1) 
	 
	
	x,y = np.meshgrid(fluxArray, deltaArray)
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\Delta^0$") 
	title = "Current";
elif plotMode == 2:  
	ax.view_init(0, -90) 
	#It's just mode 1 with a different view angle
	
	flux, delta, k, phi = np.meshgrid(fluxArray,deltaArray,kArray,phiArray);
	z = dCurrent(flux,delta, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, deltaArray) 
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\Delta^0$") 
	title = "Current_k_phi";
elif plotMode == 3:
	
	ax.view_init(30, 30) 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	z = Delta(1, k, phi);
	
	ax.set_aspect('equal'); #can be used for k-space, but this wants x-axis = y-axis
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Gap function";
elif plotMode == 4:
	ax.view_init(30, 30) 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	z =  Fermi(k, phi);
	
	ax.set_aspect('equal'); #can be used for k-space, but this wants x-axis = y-axis
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Fermi disc";
elif plotMode == 5:
	ax.view_init(40, -30) 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	
	deltaTunnel = lambda kk, pp: tunnel(kk,pp) * Delta(deltaS, kk, pp)
	
	z = deltaTunnel(k, phi);
	
	ax.set_aspect('equal'); #can be used for k-space, but this wants x-axis = y-axis
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Gap function times tunnel";  
elif plotMode == 6: 
	k, phi = np.meshgrid(kArray,phiArray)
	
	x = k * np.cos(phi)
	y = k * np.sin(phi)
	 
	
	z =  Delta(deltaS, k, phi);
	 
	ax.set_aspect('equal'); #can be used for k-space, but this wants x-axis = y-axis
	plt.xlabel("$k_x$")
	plt.ylabel("$k_y$")
	title = "Gap function times tunnel";  
#Modes 7-8-9-10 are different, because of the weird parameters being used.
#They are meant to investigate how sensitive the measurement is to the selection.
#Mode 7 is the template, so it will use the block tunnel that could've been done
#the 'original' way.
elif plotMode == 7:  
	#ax.view_init(89, -91) top
	ax.view_init(60, -70)   
	
	#FD block
	## from scatter_plot
	##  y1 =  heaviside(mu1 + mu/2 - x)  - heaviside(mu1 - mu/2 - x)
	
	mu = np.pi/4.*3.
	width = np.pi/4.
	sigmaArray = np.linspace(mu/100.,mu,N)
	dd = deltaS
	#w for weird
	wunnel = lambda ss,mm, kk,pp: Heaviside( mm + width - pp) - Heaviside(mm - width-pp)
	wCurrent = lambda ff, ss, kk, pp: Fermi(kk, pp)*wunnel(ss,mu,kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(-np.angle(Delta(dd,kk,pp)) + ff) *EnergyPart(kk,dd,pp)

	flux, sigma, k, phi = np.meshgrid(fluxArray,sigmaArray,kArray,phiArray);
	z = wCurrent(flux,sigma, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, sigmaArray) 
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\sigma$") 
	title = "Block relaxation";
elif plotMode == 8:  
	#ax.view_init(89, -91) top
	ax.view_init(60, -70)    
	
	#FD block
	## from scatter_plot
	nfd = lambda xx, sigma, mu: (np.exp( (xx-mu)/sigma) + 1)**(-1)
	##y2 =  nfd(-x, sigma, -(mu2+mu/2)) * nfd(x, sigma, mu2-mu/2)
	
	mu = np.pi/4.*3.
	width = np.pi/4.
	sigmaArray = np.linspace(mu/100.,2*mu,N)
	dd = deltaS
	#w for weird
	wunnel = lambda ss,mm, kk,pp: nfd(-pp, ss, -(mu+width)) * nfd(pp, ss, mu-width)
	wCurrent = lambda ff, ss, kk, pp: Fermi(kk, pp)*wunnel(ss,mu,kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(-np.angle(Delta(dd,kk,pp)) + ff) *EnergyPart(kk,dd,pp)

	flux, sigma, k, phi = np.meshgrid(fluxArray,sigmaArray,kArray,phiArray);
	z = wCurrent(flux,sigma, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, sigmaArray) 
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\sigma$") 
	title = "Fermi-Dirac relaxation";
elif plotMode == 9:  
	#ax.view_init(89, -91) top
	ax.view_init(60, -70)  
	
	#FD block
	## from scatter_plot 
	gaussian = lambda xx, sigma, mu: np.exp( - (xx-mu)**2 / (2 * sigma**2)) 
	
	mu = np.pi/4.*3.
	width = np.pi/4.
	sigmaArray = np.linspace(mu/100.,2*mu,N)
	dd = deltaS
	#w for weird
	wunnel = lambda ss,mm, kk,pp: gaussian(pp,ss,mu);
	wCurrent = lambda ff, ss, kk, pp: Fermi(kk, pp)*wunnel(ss,mu,kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(-np.angle(Delta(dd,kk,pp)) + ff) *EnergyPart(kk,dd,pp)

	flux, sigma, k, phi = np.meshgrid(fluxArray,sigmaArray,kArray,phiArray);
	z = wCurrent(flux,sigma, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, sigmaArray) 
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\sigma$") 
	title = "Gaussian relaxation";
elif plotMode == 10:  
	#ax.view_init(89, -91) top
	ax.view_init(60, -70)  
	
	#FD block
	## from scatter_plot 
	lorentzian = lambda xx, sigma, mu: sigma**2 / ( (xx-mu)**2 + sigma**2)
	
	mu = np.pi/4.*3.
	width = np.pi/4.
	sigmaArray = np.linspace(mu/100.,2*mu,N)
	dd = deltaS
	#w for weird
	wunnel = lambda ss,mm, kk,pp: lorentzian(pp,ss,mu);
	wCurrent = lambda ff, ss, kk, pp: Fermi(kk, pp)*wunnel(ss,mu,kk,pp)*dkdphi*np.abs(Delta(dd, kk, pp))*deltaS*np.sin(-np.angle(Delta(dd,kk,pp)) + ff) *EnergyPart(kk,dd,pp)

	flux, sigma, k, phi = np.meshgrid(fluxArray,sigmaArray,kArray,phiArray);
	z = wCurrent(flux,sigma, k, phi).sum(axis=-1).sum(axis=-1) 
	
	
	x,y = np.meshgrid(fluxArray, sigmaArray) 
	
	plt.xlabel("$\Phi$")
	plt.ylabel("$\sigma$") 
	title = "Lorentzian relaxation"; 
else:
	raise Exception("Unknown plot mode.");   


if alpha < 360 and beta < 360:
	ax.view_init(alpha, beta) 	
	
plt.title("%s, N=%d, d=%d, m=%d, p=%d, k=%d, b=%d" % (title, N,gapfunction, mechanism, plotMode, fermi, band)) 

if scatter:
	ax.scatter(x, y, z,c=50.*(z+np.min(z))/np.max(z), cmap=cm.winter);
elif plotMode == 6: #sorry, this one is atypical.
 	levels = MaxNLocator(nbins=40).tick_values(z.min(), z.max())

	cmap = plt.get_cmap('winter')
	norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

	cf = ax.contourf(x, y, z, levels=levels, cmap=cmap)
	fig.colorbar(cf, ax=ax)   
	
	xf = fermiSurface * np.cos(phiArray)
	yf = fermiSurface * np.sin(phiArray)
	
	ax.plot(xf,yf, 'r--');
elif gamma > 0.95:
	ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.summer,linewidth=0) 
else:
	ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.summer,linewidth=0, alpha=gamma) 

print >> sys.stderr, "Elapsed time %2.3f" % (time.time() - startTime)
if silent:
	print >> sys.stderr, "Silent Mode; no plotting, no saving.";
	print >> sys.stderr, "Title [%s]." % title;
	plt.close();
if filename != "default.png":	
	fig.savefig(filename)
else:
	plt.show()