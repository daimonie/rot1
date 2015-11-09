#For the Presentation, a few plots for teh scatterers
import numpy as np
import matplotlib.pyplot as plt


maximus = 2*np.pi;
N = 1000;
x = np.arange(0, maximus, maximus/N)

transparency = 0.25


heaviside = lambda xx: 1. * (xx>=0)

setzero = lambda xx: 1. - heaviside(xx -  maximus*(1.-1./N))

nfd = lambda xx, sigma, mu: (np.exp( (xx-mu)/sigma) + 1)**(-1)
gaussian = lambda xx, sigma, mu: np.exp( - (xx-mu)**2 / (2 * sigma**2))
lorentzian = lambda xx, sigma, mu: sigma**2 / ( (xx-mu)**2 + sigma**2)

mu = np.pi/6.
sigma = mu/15.

mu1 = mu * 1.
mu2 = mu * 2.
mu3 = mu * 3.
mu4 = mu * 4.

#Note: Fill has to end at zero.
y1 =  heaviside(mu1 + mu/2 - x)  - heaviside(mu1 - mu/2 - x)
y2 =  nfd(-x, sigma, -(mu2+mu/2)) * nfd(x, sigma, mu2-mu/2)
y3 =  lorentzian(x,sigma,mu3)
y4 =  gaussian(x,sigma,mu4)

#make them even
y1 = y1 + y1[::-1]
y2 = y2 + y2[::-1]
y3 = y3 + y3[::-1]
y4 = y4 + y4[::-1]

#normalisation
y1 = y1 / np.max(y1)
y2 = y2 / np.max(y2)
y3 = y3 / np.max(y3)
y4 = y4 / np.max(y4)


zeroArray = setzero(x)
zeroArray = np.multiply(zeroArray[::-1], zeroArray)

plt.fill(x, y1 * zeroArray, 'r', alpha=transparency, label='Block');
plt.fill(x, y2 * zeroArray, 'g', alpha=transparency, label='Fermi-Dirac block');
plt.fill(x, y3 * zeroArray, 'b', alpha=transparency, label='Lorentzian');
plt.fill(x, y4 * zeroArray, 'm', alpha=transparency, label='Gaussian'); 
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.show()