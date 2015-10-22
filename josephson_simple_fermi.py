import numpy as np;
import scipy as sp;
from scipy.constants import *
#params are nice...
import argparse as argparse
#plotting
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
#script
#arguments
parser = argparse.ArgumentParser(prog="python main.py",
  description = "Takes Delta(k) and plots phase/amp graphs. Use -help to see more information."); 
parser.add_argument('-n', '--number', help='Resolution onf the calculation.', default = 10, action='store', type = int);
args = parser.parse_args();


N = args.number;
#some constants
hbar =  physical_constants["Planck constant over 2 pi"][0]
mass =  physical_constants["electron mass"][0]
ev =  physical_constants["electron volt"][0]
lattice = 3.787 * physical_constants["Angstrom star"][0] #(From Wikipedia)


eta = hbar**2 / 2 / mass / ev / lattice**2;
deltas = 3e-3;

kFermi = 3*np.pi/4;
#linspace is really slow. I'm serious. Let's not use it :-)
#arange is about seven times faster. See test.py

#Due to time constraints, we're going for dirty for loops. Sorry. 
for flux in np.arange(-5, 5, 10.0/N):
	for delta in np.arange(2*deltas/N, 2*deltas, 2*deltas/N):
		current = 0;
		for phi in np.arange(0,N, 2 * np.pi / N): 
			for k in np.arange(0, kFermi, kFermi/N):
				kx = kFermi * np.cos(phi);
				ky = kFermi * np.sin(phi); 
				
				Delta = delta * (np.cos(kx) - np.cos(ky));
				
				norm = np.abs(Delta)
				arg = np.angle(Delta);
			
				energyL = ( eta*k**2 + deltas**2);
				energyR = ( eta*k**2 + norm**2);
				
				current += deltas*norm/(energyL+energyR)/energyL/energyR*np.sin(arg+flux*2*np.pi)
		print "%2.3e\t%2.3e\t%2.3e" % (flux,delta,current)
	print "" 	
			
			
			
			
			
			
			
			
			
