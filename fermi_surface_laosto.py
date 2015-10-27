from __future__ import division
from math import pi, sin, cos, exp
import numpy as np
import tinyarray as ta
from scipy.optimize import newton
from scipy.linalg import eigvalsh
import sys
 

from laosto import *
from conductivity import *

res = 5;

kF0 = [0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi, -1, -1]
phis = np.linspace(0., 2*pi, res, endpoint=False)

system = LAOSTO(mu=0, H=30, theta=np.pi/4, g=5, gL=1, tl1=340, tl2=340, th=12.5,
               td=12.5, dE=60, dZ=15, dSO=5)

def fermi_surface(system):
    kFs, indices = kF_angle(system, phis, kF0)  
    return kFs, indices

print "\tphi\tk_x\tk_y\tk"
Fermi, Indices = fermi_surface(system);
for n in range(0, Indices.shape[0]):
	start = Indices[n,1]
	end = Indices[n,2]
	print "%d-th band from %d to %d.\n" % (n, start, end)
	for i in range(start,end):
		kx = Fermi[i, 0]
		ky = Fermi[i, 1]
		print "\t%2.3e\t%2.3e\t%2.3e\t%2.3e" % (phis[i-n*res], kx, ky, (kx**2+ky**2)**0.5)