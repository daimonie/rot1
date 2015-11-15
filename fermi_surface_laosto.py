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


#res = 5;
sres = 10
eres = 1000
for res in range(sres,eres): 
	stime = time.time();
	
	kF0 = [0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi, -1, -1]
	phis = np.linspace(0., 2*pi, res, endpoint=False)

	system = LAOSTO(mu=0, H=0, theta=np.pi/4, g=5, gL=1, tl1=340, tl2=340, th=12.5,
		td=12.5, dE=60, dZ=15, dSO=5)
	kFermis, indices = fermi_surface(system, phis, kF0) 
	dt = time.time() - stime;
	print "Calculated [%d] in [%2.3f] seconds" % (res, dt)
	for b in range(0,4):
		start = indices[b, 1]
		end = indices[b, 2]

		fermiLevel = kFermis[start:end,:]
				
		fermiSurface = ((fermiLevel**2).sum(axis=1))**0.5
		
		x = fermiSurface * np.cos(phis);
		y = fermiSurface * np.cos(phis);
		
		x = x.sum(axis=0)
		y = y.sum(axis=0)
		print "\t%d\t%2.3f\t%2.3f\t%2.3f" % (b, x, y, ( x**2 + y**2 )**0.5)
		
