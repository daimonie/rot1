#File received from Beenakker group
from __future__ import division
from math import pi, sin, cos, exp
import numpy as np
import tinyarray as ta
from scipy.optimize import newton
from scipy.linalg import eigvalsh
import sys
import time

from laosto import *
from conductivity import *


def fermi_surface(ss, pp, kk):
    kFs, indices = kF_angle(system, phis, kF0)  
    return kFs, indices
 

for magnetic in np.arange(0,30,0.1):
	print "Magnetic field H=%2.3f" % magnetic 
	res = 50
	stime = time.time();
	
	kF0 = [0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi, -1, -1]
	phis = np.linspace(0., 2*pi, res, endpoint=False)

	system = LAOSTO(mu=0, H=magnetic, theta=np.pi/4, g=5, gL=1, tl1=340, tl2=340, th=12.5,
		td=12.5, dE=60, dZ=15, dSO=5)
	kFermis, indices = fermi_surface(system, phis, kF0) 
	dt = time.time() - stime; 
	xtotal = 0.
	ytotal = 0.
	for b in range(0,4):
		start = indices[b, 1]
		end = indices[b, 2]

		fermiLevel = kFermis[start:end,:]
				
		fermiSurface = ((fermiLevel**2).sum(axis=1))**0.5
		
		x = fermiSurface * np.cos(phis);
		y = fermiSurface * np.cos(phis);
		
		x = x.sum(axis=0)
		y = y.sum(axis=0) 
		
		xtotal += x
		ytotal += y
	xtotal /= 4.
	ytotal /= 4

	print "\t%d gives %2.3f" % (res, ( xtotal**2 + ytotal**2 )**0.5) 		
