from __future__ import division
from math import pi, sin, cos, exp
import numpy as np
import tinyarray as ta
from scipy.optimize import newton
from scipy.linalg import eigvalsh


class SimpleNamespace(object):
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


def n_F(E, T):
    """Fermi occupation number.

    Parameters:
    -----------
    E : float
        energy in meV
    T : float
        temperature in K
    """
    # k_B * Kelvin = 0.086173 meV
    if T == 0:
        return E < 0
    try:
        return 1/(exp(E/(T*0.086173)) + 1)
    except OverflowError:
        return 0


class LAOSTO(object):
    def __init__(self, **kwargs):
        self.set_params(**kwargs)
        self.Nbands = 6


    def H(self, kx, ky):
        """3-band LAO/STO Hamiltonian.

        Parameters:
        -----------
        kx : float
            x-lattice momentum in units of 1/lattice constant
        ky : float
           y-lattice momentum in units of 1/lattice constant

        """
        p = self.p
        epsxy = 2*p.tl1*(2-cos(kx)-cos(ky)) - p.dE - p.mu
        epsxz = 2*p.tl2*(1-cos(kx)) + 2*p.th*(1-cos(ky)) - p.mu
        epsyz = 2*p.th*(1-cos(kx)) + 2*p.tl2*(1-cos(ky)) - p.mu
        delta = 2*p.td*sin(kx)*sin(ky)
        dZ_x = p.dZ*sin(kx)
        dZ_y = p.dZ*sin(ky)
        # mu_B * Tesla = 0.058 meV
        B = 0.058*p.H
        Bx = B*cos(p.theta)*cos(p.alpha)
        By = B*sin(p.theta)*cos(p.alpha)
        Bz = B*sin(p.alpha)

        BxS = p.g*Bx/2
        ByS = p.g*By/2
        BzS = p.g*Bz/2

        BxL = p.gL*Bx
        ByL = p.gL*By
        BzL = p.gL*Bz

        return ta.array(
            [[       epsxy - BzS,     -BxS + 1j*ByS,            1j*dZ_y - 1j*BxL,                  1j/2*p.dSO,            1j*ByL + 1j*dZ_x,                  -1/2*p.dSO],
             [     -BxS - 1j*ByS,       epsxy + BzS,                  1j/2*p.dSO,            1j*dZ_y - 1j*BxL,                   1/2*p.dSO,            1j*ByL + 1j*dZ_x],
             [ -1j*dZ_y + 1j*BxL,       -1j/2*p.dSO,                 epsxz - BzS,               -BxS + 1j*ByS, delta + 1j/2*p.dSO - 1j*BzL,                           0],
             [       -1j/2*p.dSO, -1j*dZ_y + 1j*BxL,               -BxS - 1j*ByS,                 epsxz + BzS,                           0, delta - 1j/2*p.dSO - 1j*BzL],
             [ -1j*ByL - 1j*dZ_x,         1/2*p.dSO, delta - 1j/2*p.dSO + 1j*BzL,                           0,                 epsyz - BzS,               -BxS + 1j*ByS],
             [        -1/2*p.dSO, -1j*ByL - 1j*dZ_x,                           0, delta + 1j/2*p.dSO + 1j*BzL,               -BxS - 1j*ByS,                 epsyz + BzS]], complex)


    def dHdkx(self, kx, ky):
        """Derivative of the Hamiltonian with respect to kx.

        Parameters:
        -----------
        kx : float
            x-lattice momentum in units of 1/lattice constant
        ky : float
            y-lattice momentum in units of 1/lattice constant
        """

        p = self.p
        depsxy_dkx = 2*p.tl1*sin(kx)
        depsxz_dkx = 2*p.tl2*sin(kx)
        depsyz_dkx = 2*p.th*sin(kx)
        ddelta_dkx = 2*p.td*cos(kx)*sin(ky)
        ddZ_x = p.dZ*cos(kx)

        return ta.array(
            [[        depsxy_dkx,             0,            0,             0,         1j*ddZ_x,             0],
             [                 0,    depsxy_dkx,            0,             0,                0,      1j*ddZ_x],
             [                 0,             0,   depsxz_dkx,             0,       ddelta_dkx,             0],
             [                 0,             0,            0,    depsxz_dkx,                0,    ddelta_dkx],
             [        - 1j*ddZ_x,             0,   ddelta_dkx,             0,       depsyz_dkx,             0],
             [                 0,    - 1j*ddZ_x,            0,    ddelta_dkx,                0,    depsyz_dkx]], complex)


    def dHdky(self, kx, ky):
        """Derivative of the Hamiltonian with respect to ky.

        Parameters:
        -----------
        kx : float
            x-lattice momentum in units of 1/lattice constant
        ky : float
            y-lattice momentum in units of 1/lattice constant
        p : SimpleNamespace of floats
        """
        p = self.p
        depsxy_dky = 2*p.tl1*sin(ky)
        depsxz_dky = 2*p.th*sin(ky)
        depsyz_dky = 2*p.tl2*sin(ky)
        ddelta_dky = 2*p.td*sin(kx)*cos(ky)
        ddZ_y = p.dZ*cos(ky)

        return ta.array(
            [[        depsxy_dky,             0,     1j*ddZ_y,             0,                0,             0],
             [                 0,    depsxy_dky,            0,      1j*ddZ_y,                0,             0],
             [        - 1j*ddZ_y,             0,   depsxz_dky,             0,       ddelta_dky,             0],
             [                 0,    - 1j*ddZ_y,            0,    depsxz_dky,                0,    ddelta_dky],
             [                 0,             0,   ddelta_dky,             0,       depsyz_dky,             0],
             [                 0,             0,            0,    ddelta_dky,                0,    depsyz_dky]], complex)


    def set_params(self, mu=0, H=0, theta=0, alpha=0, g=2, gL=1, tl1=340,
                   tl2=340, th=12.5, td=12.5, dE=60, dZ=15, dSO=5):
        """Parameters for t2g Hamiltonian for LAO/STO

        Parameters
        ----------
        mu : float
            chemical potential of the t2g electrons
        H : float
            in-plane magnetic-field in Tesla
        theta : float
            angle between the magnetic field and the x-axis
        alpha : float
            out-of-plane angle of the magnetic field
        g : float
            Lande g-factor
        gL : float
            "fake Lande g-factor" for B dot L coupling,
            where L is the orbital mometum of d_xy, xz, yz
        tl1 : float
            inverse light mass of the xy-band
        tl2 : float
            inverse light mass of the upper bands
        th : float
            inverse heavy mass
        td : float
            hybridization of the d_xz and d_yz bands
        dE : float
            split-of energy of the d_xy band
        dSO : float
            atomic spin-orbit splitting
        dZ : float
            interface inversion breaking
        """
        self.p = SimpleNamespace(mu=mu, H=H, theta=theta, alpha=alpha, g=g,
                                 gL=gL, tl1=tl1, tl2=tl2, th=th, td=td, dE=dE,
                                 dZ=dZ, dSO=dSO)


    def calc_n(self, T, n_cut=1e-4, Nk=100, kmax=pi/4):
        """integrate carrier density

        Parameters:
        -----------
        params : SimpleNamespace
            Hamiltonian parameters
        T : float
            temperature in Kelvin
        n_cut : float
            cutoff density. if > 0 integration limit kmax is chosen
            such that n(kmax) < n_cut. Nk is chosen accordingly.
        Nk : int
            Number of descrete k-values in each direction.
            ignored if n_cut > 0
        kmax : float
            integration limit in k-space.
            ignored if n_cut > 0
        """
        if n_cut > 0:
            def n_kmax(k):
                Emin = min(eigvalsh(self.H(-k, 0))[0],
                           eigvalsh(self.H(k, 0))[0],
                           eigvalsh(self.H(0, -k))[0],
                           eigvalsh(self.H(0, k))[0])
                return n_F(Emin, T)

            if n_kmax(0.15*pi) < n_cut:
                Nk, kmax = 150, 0.15*pi
            elif n_kmax(0.2*pi) < n_cut:
                Nk, kmax = 200, 0.2*pi
            elif n_kmax(0.25*pi) < n_cut:
                Nk, kmax = 250, 0.25*pi
            elif n_kmax(0.3*pi) < n_cut:
                Nk, kmax = 300, 0.3*pi
            elif n_kmax(0.35*pi) < n_cut:
                Nk, kmax = 350, 0.35*pi
            elif n_kmax(0.40*pi) < n_cut:
                Nk, kmax = 400, 0.4*pi
            elif n_kmax(0.50*pi) < n_cut:
                Nk, kmax = 500, 0.5*pi
            else:
                Nk, kmax = 1000, pi

        n = 0

        kx_list = np.linspace(-kmax, kmax, Nk)
        ky_list = np.linspace(-kmax, kmax, Nk)
        deltak = 2*kmax/(Nk-1)
        for kx in kx_list:
            for ky in ky_list:
                E_list = eigvalsh(self.H(kx, ky))
                n += sum([n_F(E_list[i], T) for i in range(6)])

        n *= deltak**2 / (4*pi**2)

        return n


    def calc_mu(self, n, T, tol=1e-4, maxiter=50):
        def offset(mu, n):
            self.p.mu = mu
            return (self.calc_n(T, n_cut=tol)-n)/n

        self.p.mu = newton(offset, self.p.mu, args=(n,),
                           tol=1e-4, maxiter=50)
        return self.p.mu
