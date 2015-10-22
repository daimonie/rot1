from __future__ import division
from math import pi, sin, cos, exp
import numpy as np
import tinyarray as ta
from scipy.optimize import newton
from scipy.linalg import eigvalsh
import sys

sys.path.append('/home/bovenzi/boltzmann')

from laosto import *
from conductivity import *

kF0 = [0.5*pi, 0.5*pi, 0.5*pi, 0.5*pi, -1, -1]
phis = np.linspace(0., 2*pi, 200, endpoint=False)

system = LAOSTO(mu=0, H=30, theta=np.pi/4, g=5, gL=1, tl1=340, tl2=340, th=12.5,
               td=12.5, dE=60, dZ=15, dSO=5)

def fermi_surface(system):
    kFs, indices = kF_angle(system, phis, kF0)  
    return kFs, indices




    
