#File received from Beenakker Group
from __future__ import division
import numpy as np
import tinyarray as ta
from scipy.linalg import eigvalsh, eigh
from scipy.optimize import newton
from math import pi, sin, cos, exp

# Marching squares algorithm
from scipy.optimize import brentq

marching_step_dict = {
    0b0000: ta.array([1, 0]),
    0b0001: ta.array([-1, 0]),
    0b0010: ta.array([0, -1]),
    0b0100: ta.array([1, 0]),
    0b1000: ta.array([0, 1]),
    0b0011: ta.array([-1, 0]),
    0b0110: ta.array([0, -1]),
    0b1100: ta.array([1, 0]),
    0b1001: ta.array([0, 1]),
    0b0111: ta.array([-1, 0]),
    0b1110: ta.array([0, -1]),
    0b1101: ta.array([1, 0]),
    0b1011: ta.array([0, 1])
}


def marching_step(cell, func, middle, d_ij):
    """Return the direction to go to reach the next cell in the marching
    squares algorithm. `func`, `middle` and `d_ij` are only accessed
    when they are necessary to resolve saddle-point ambiguities
    (i.e. ``cell == 0b0101`` or ``cell == 0b1010``).

    """
    try:
        return marching_step_dict[cell]
    except KeyError:
        if cell == 0b0101:
            if func(*middle) < 0:
                if d_ij == (0, 1):
                    return ta.array([1, 0])
                elif d_ij == (0, -1):
                    return ta.array([-1, 0])
                else:
                    raise RuntimeError("I shouldn't have come here from this"
                                       "direction...")
            else:
                if d_ij == (0, 1):
                    return ta.array([-1, 0])
                elif d_ij == (0, -1):
                    return ta.array([1, 0])
                else:
                    raise RuntimeError("I shouldn't have come here from this "
                                       "direction...")
        elif cell == 0b1010:
            if func(*middle) < 0:
                if d_ij == (1, 0):
                    return ta.array([0, -1])
                elif d_ij == (-1, 0):
                    return ta.array([0, 1])
                else:
                    raise RuntimeError("I shouldn't have come here from this "
                                       "direction...")
            else:
                if d_ij == (1, 0):
                    return ta.array([0, 1])
                elif d_ij == (-1, 0):
                    return ta.array([0, -1])
                else:
                    raise RuntimeError("I shouldn't have come here from this "
                                       "direction...")
        else:
            raise RuntimeError("cell " + bin(cell) + " shouldn't happen ...")


class MarchingCellInterpolate(object):
    def __init__(self, ij, kxs, kys, func):
        self.ij = ta.array(ij)
        self.kxs = kxs
        self.kys = kys
        self.func = func
        self.d_ij = None

        i, j = ij
        self.values = (self.func(kxs[i], kys[j]),
                       self.func(kxs[i+1], kys[j]),
                       self.func(kxs[i+1], kys[j+1]),
                       self.func(kxs[i], kys[j+1]))
        self.cell = 0
        for i in xrange(4):
            self.cell += (self.values[i] < 0) << i


    def next(self):
        kxs, kys = self.kxs, self.kys
        d_ij = marching_step(self.cell, self.func,
                             (0.5 * (kxs[self.ij[0]] + kxs[self.ij[0]+1]),
                              0.5 * (kys[self.ij[1]] + kys[self.ij[1]+1])),
                             self.d_ij)
        ij = self.ij + d_ij

        i, j = ij
        cell_slow = (((self.func(kxs[i], kys[j]) < 0) << 0) +
                     ((self.func(kxs[i+1], kys[j]) < 0) << 1) +
                     ((self.func(kxs[i+1], kys[j+1]) < 0) << 2) +
                     ((self.func(kxs[i], kys[j+1]) < 0) << 3))

        if d_ij[0] == 1:
            self.values = (self.values[1], self.func(kxs[i+1], kys[j]),
                           self.func(kxs[i+1], kys[j+1]), self.values[2])
            kx = kxs[i]
            ky = kys[j] - self.values[0] * ((kys[j+1] - kys[j]) /
                                            (self.values[3] - self.values[0]))
        elif d_ij[0] == -1:
            self.values = (self.func(kxs[i], kys[j]), self.values[0],
                           self.values[3], self.func(kxs[i], kys[j+1]))
            kx = kxs[i+1]
            ky = kys[j] - self.values[1] * ((kys[j+1] - kys[j]) /
                                            (self.values[2] - self.values[1]))
        elif d_ij[1] == 1:
            self.values = (self.values[3], self.values[2],
                           self.func(kxs[i+1], kys[j+1]),
                           self.func(kxs[i], kys[j+1]))
            kx = kxs[i] - self.values[0] * ((kxs[i+1] - kxs[i]) /
                                            (self.values[1] - self.values[0]))
            ky = kys[j]
        else:
            self.values = (self.func(kxs[i], kys[j]),
                           self.func(kxs[i+1], kys[j]),
                           self.values[1], self.values[0])
            kx = kxs[i] - self.values[3] * ((kxs[i+1] - kxs[i]) /
                                            (self.values[2] - self.values[3]))
            ky = kys[j+1]

        cell = 0
        for i in xrange(4):
            cell += (self.values[i] < 0) << i
        assert cell == cell_slow

        self.cell = cell
        self.ij = ij
        self.d_ij = d_ij

        return kx, ky


class MarchingCellBrent(object):
    def __init__(self, ij, kxs, kys, func):
        self.ij = ta.array(ij)
        self.kxs = kxs
        self.kys = kys
        self.func = func
        self.d_ij = None

        i, j = ij
        self.cell = (((self.func(kxs[i], kys[j]) < 0) << 0) +
                     ((self.func(kxs[i+1], kys[j]) < 0) << 1) +
                     ((self.func(kxs[i+1], kys[j+1]) < 0) << 2) +
                     ((self.func(kxs[i], kys[j+1]) < 0) << 3))


    def next(self):
        kxs, kys = self.kxs, self.kys
        d_ij = marching_step(self.cell, self.func,
                             (0.5 * (kxs[self.ij[0]] + kxs[self.ij[0]+1]),
                              0.5 * (kys[self.ij[1]] + kys[self.ij[1]+1])),
                             self.d_ij)
        ij = self.ij + d_ij

        i, j = ij
        cell_slow = (((self.func(kxs[i], kys[j]) < 0) << 0) +
                     ((self.func(kxs[i+1], kys[j]) < 0) << 1) +
                     ((self.func(kxs[i+1], kys[j+1]) < 0) << 2) +
                     ((self.func(kxs[i], kys[j+1]) < 0) << 3))

        kx, ky = None, None

        if d_ij[0] == 1:
            cell = (((self.cell & 0b0010) >> 1) + ((self.cell & 0b0100) << 1) +
                    ((self.func(kxs[i+1], kys[j]) < 0) << 1) +
                    ((self.func(kxs[i+1], kys[j+1]) < 0) << 2))

            if self.cell != 0b0000:
                kx = kxs[i]
                ky = brentq(lambda y: self.func(kx, y), kys[j], kys[j+1])
        elif d_ij[0] == -1:
            cell = (((self.cell & 0b0001) << 1) + ((self.cell & 0b1000) >> 1) +
                    ((self.func(kxs[i], kys[j]) < 0) << 0) +
                    ((self.func(kxs[i], kys[j+1]) < 0) << 3))

            if self.cell != 0b0000:
                kx = kxs[i+1]
                ky = brentq(lambda y: self.func(kx, y), kys[j], kys[j+1])
        elif d_ij[1] == 1:
            cell = (((self.cell & 0b0100) >> 1) + ((self.cell & 0b1000) >> 3) +
                    ((self.func(kxs[i+1], kys[j+1]) < 0) << 2) +
                    ((self.func(kxs[i], kys[j+1]) < 0) << 3))

            if self.cell != 0b0000:
                ky = kys[j]
                kx = brentq(lambda x: self.func(x, ky), kxs[i], kxs[i+1])
        else:
            cell = (((self.cell & 0b0001) << 3) + ((self.cell & 0b0010) << 1) +
                    ((self.func(kxs[i], kys[j]) < 0) << 0) +
                    ((self.func(kxs[i+1], kys[j]) < 0) << 1))

            if self.cell != 0b0000:
                ky = kys[j+1]
                kx = brentq(lambda x: self.func(x, ky), kxs[i], kxs[i+1])

        assert cell == cell_slow

        self.cell = cell
        self.ij = ij
        self.d_ij = d_ij

        return kx, ky


def find_contours(kx1, kx2, ky1, ky2, Nx, Ny, E, refinement=10,
                  method='brent'):
    """Find the all contours where the function E(kx, ky) is 0.

    This function finds all contours within a rectangular region determined
    by kxrange and kyrange (see below for details). It is guaranteed to find
    any contour that is larger than the grid spacing specified there. This sets
    the limit for the fineness of the initial mesh, on which the function E is
    evaluated.

    After finding existing contours, the contours themselves can be computed on
    a finer mesh (using the parameter `refinement`) evaluating E only for mesh
    points near the contour, using the marching squares algorithm.

    Parameters
    ----------
    kx1 : float
    kx2 : float
    ky1 : float
    ky2 : float
        Search for contours within [`kx1`, `kx2`] x [`ky1`, ky2`]. The contours
        must lie completely within this area.
    Nx : int
    Ny : int
        Number of points in x and y-direction, respectively, in the initial
        grid.
    E : function of two parameters
        `E` = E(kx, ky) is a function taking two floats, and returning one
        float. The contours where this function is 0 are returned.
    refinement : int
        The factor with which Nx and Ny are multiplied for the final, fine
        grid. Defaults to 10.
    method : ("brent", "interpolate")
        The method with which the intersection of the contour with the
        individual cells is found. "brent" finds the points on the contour up
        to numerical precision using Brent's method (some advanced Newton-type
        root finding), "interpolate" uses linear interpolation. Defaults to
        "brent".

    Returns
    -------
    contours : list of numpy arrays
        List of all individual contours. Each contour is returned as a numpy
        array of shape (Ncontour, 2), where Ncontour is the number of points
        in the individual contour.
    """

    if method == 'brent':
        MarchingCell = MarchingCellBrent
    elif method == 'interpolate':
        MarchingCell = MarchingCellInterpolate

    # First, consider the full range of k-space with the coarse grid
    # in order to find existing contours

    kx_range = np.linspace(kx1, kx2, Nx)
    ky_range = np.linspace(ky1, ky2, Ny)

    gridvals = np.vectorize(E)(*np.meshgrid(kx_range, ky_range,
                                            indexing='ij')) < 0
    cells = ((gridvals[:-1, :-1] << 0) + (gridvals[1:, :-1] << 1) +
             (gridvals[1:, 1:] << 2) + (gridvals[:-1, 1:] << 3))

    def find_one_contour(cells):
        # start with a particular cell (needed for later, so that the
        # refinement works easily)
        boundary_cells = np.logical_or(cells == 0b0010,
                                       cells == 0b0110,
                                       cells == 0b1110).nonzero()

        try:
            ij = start = tuple(np.array(boundary_cells)[:, 0])
        except IndexError:
            return None

        d_ij = None
        while True:
            d_ij = marching_step(cells[ij], E,
                                 (0.5 * (kx_range[ij[0]] + kx_range[ij[0]+1]),
                                  0.5 * (ky_range[ij[1]] + ky_range[ij[1]+1])),
                                 d_ij)
            if cells[ij] not in (0b0101, 0b1010):
                cells[ij] = 0b0000
            ij = tuple(ta.array(ij) + d_ij)

            if ij == start:
                break

        return start

    contour_starts = []
    start = find_one_contour(cells)
    while start:
        contour_starts.append(start)
        start = find_one_contour(cells)

    # Having found one point of each contour, now find the individual
    # contours on a refined grid, computing function values only along
    # the contour (as opposed to all of space as in the first step)

    kx_range = np.linspace(kx1, kx2, (Nx-1) * refinement + 1)
    ky_range = np.linspace(ky1, ky2, (Ny-1) * refinement + 1)
    contours = []
    for start in contour_starts:
        kFs = []
        i, j = start
        cell = MarchingCell((i * refinement, j * refinement),
                            kx_range, ky_range, E)

        while cell.cell == 0b0000:
            cell.next()
        start = cell.ij

        while True:
            kFs.append(cell.next())
            if cell.ij == start:
                break
        contours.append(np.array(kFs))

    return contours


def kF_marchingsquares(sys, kx1, kx2, ky1, ky2, Nx, Ny, refinement=10, E=0):
    """Determine the Fermi surface using the marching squares algorithm

    Parameters:
    -----------
    params : SimpleNamespace
        Hamiltonian parameters
    kx1 : float
    kx2 : float
    ky1 : float
    ky2 : float
        Find contours in area [`kx1`, `kx2`] x [`ky1`, `ky2`].
    Nx : int
    Ny : int
        Number of points used to initially discretize k-space in kx and
        ky-direction. Note that only contours will be found for which at least
        1 point lies inside the contour. Hence, finding a contour is in general
        only guaranteed if (kx2-k1)/Nx and/or (ky2-ky1)/Ny are smaller than the
        "diameter" of the contour.
    refinement : int
        The factor with which Nx and Ny are multiplied for the final, fine
        grid. Defaults to 10.
    E : float
        energy H(kF, params) = E
    tol : float
        tolarance for the convergence of the Newton iteration
    maxiter : int
        maximum number of iterations per kF

    Returns:
    --------
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    """
    def zero(kx, ky, n_band):
        return np.sort(eigvalsh(sys.H(kx, ky)))[n_band] - E

    kFs = []
    bands = []
    Nbands = sys.Nbands

    for n_band in xrange(Nbands):
        contours = find_contours(kx1, kx2, ky1, ky2, Nx, Ny,
                                 lambda kx, ky: zero(kx, ky, n_band),
                                 refinement)

        kFs.extend(contours)
        bands.extend([n_band] * len(contours))

    index = np.zeros(shape=(len(kFs), 3), dtype = int)
    index[:, 0] = bands
    index[1:, 1] = [len(kFs[i]) for i in xrange(len(bands)-1)]
    index[:, 2] = [len(kFs[i]) for i in xrange(len(bands))]
    index[:, 1:3] = np.cumsum(index[:, 1:3], axis=0)

    if len(kFs):
        return np.concatenate(kFs), index
    else:
        return np.zeros(shape=(0, 2), dtype=float), index


def kF_angle(sys, phis, kF0, E=0, tol=1e-5, maxiter=30):
    """determine Fermi momenta through Newton iteration. Only works for
    bands that have a single minimum at k=0.

    Parameters:
    -----------
    sys : BoltzmannSystem
        Hamiltonian
    phis : np.array or list of np.arrays
        if np.array: array of discritization angles for which kF is determined,
        if list: list of np.arrays with angles for each band.
    kF0 : list of floats
        initial values for Newton iteration
    E : float
        energy H(kF, params) = E
    tol : float
        tolarance for the convergence of the Newton iteration
    maxiter : int
        maximum number of iterations per kF

    Returns:
    --------
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    """
    def zero(k, phi, n_band):
        return np.sort(eigvalsh(sys.H(k*cos(phi), k*sin(phi))))[n_band] - E

    def newton_kF(kF0, phi, n_band):
        return newton(zero, kF0, args=(phi, n_band), tol=tol, maxiter=maxiter)

    bands = np.arange(sys.Nbands)[eigvalsh(sys.H(0, 0)) < E]
    kFs = []

    if type(phis) is np.ndarray:
        phis = [phis for i in xrange(sys.Nbands)]

    for n_b in bands:
        # kFs at phis[0]
        kF0n = kF0[n_b] if kF0[n_b] > 0 else 0.01*pi
        kF = np.zeros(shape=(len(phis[n_b]),), dtype=float)

        oldkF = kF[0] = kF0n  # just as a start for the first Newton iteration
        for i, phi in enumerate(phis[n_b]):
            try:
                kF[i] = newton_kF(oldkF, phi, n_b)
            except RuntimeError:
                try:
                    kF[i] = newton_kF(kF[0]/10, phi, n_b)
                except RuntimeError:
                    try:
                        kF[i] = newton_kF(kF[0]*10, phi, n_b)
                    except RuntimeError:
                        raise RuntimeError(
                            "failed to find k_F for band {0}".format(n_b))
            oldkF = kF[i]

        kFs.append(np.vstack([kF * np.cos(phis[n_b]),
                              kF * np.sin(phis[n_b])]).T)

    indices = np.zeros(shape=(len(kFs), 3), dtype = int)
    indices[:, 0] = bands
    indices[:, 1] = [sum([len(phis[bands[j]]) for j in xrange(i)])
                     for i in xrange(len(bands))]
    indices[:, 2] = [sum([len(phis[bands[j]]) for j in xrange(i+1)])
                     for i in xrange(len(bands))]

    if len(kFs):
        return np.concatenate(kFs), indices
    else:
        return np.zeros(shape=(0, 2), dtype=float), indices


def norm(k):
    """norm of a vector"""
    return np.sqrt(np.sum(k**2))


def dl(kFs, indices):
    """(approximate) line elements

    Parameters:
    -----------
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.

    Returns:
    --------
    dls : np.array((len(kFs),)) [i, n]
        list of line elements at angle phis[i], band bands[n]
    """

    dls = np.zeros(shape=(len(kFs),), dtype=float)

    for n_b, start, stop in indices:
        dls[start:stop-1] = np.sqrt((kFs[start:stop-1, 0] -
                                     kFs[start+1:stop, 0])**2 +
                                    (kFs[start:stop-1, 1] -
                                     kFs[start+1:stop, 1])**2)
        dls[stop-1] = np.sqrt((kFs[stop-1, 0] - kFs[start, 0])**2 +
                              (kFs[stop-1, 1] - kFs[start, 1])**2)

    return dls


def Es_psis(sys, kFs, indices):
    """Energies and eigenstates.

    Parameters:
    -----------
    sys : BoltzmannSystem
        system definition. Must implement the Hamiltonian ``H(kx, ky)``.
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.

    """

    Es = np.zeros(shape=(len(kFs),), dtype=float)
    psis = np.zeros(shape=(len(kFs), sys.Nbands), dtype='complex128')

    for n_b, start, stop in indices:
        for i in xrange(start, stop):
            kx = kFs[i, 0]
            ky = kFs[i, 1]
            Ens, psins = eigh(sys.H(kx, ky))
            Es[i] = Ens[n_b]
            psis[i, :] = psins[:, n_b]

    return Es, psis


def v(sys, kFs, indices, method="feynman-hellman", dk=1e-5):
    """Quasi-particle velocities

    Parameters:
    -----------
    sys : BoltzmannSystem
        system definition. Must implement the Hamiltonian ``H(kx, ky)``.
        Additionally, if ``method == "feynman-hellman``, it must also
        implement the derivatives of the Hamiltonian with respect to
        kx and ky: ``dHdkx(kx, ky)`` and ``dHdky(kx, ky)``, respectively.
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    method : ("feynman-hellman", "numerical")
        method to compute the velocities. "feynman-hellman" uses the
        derivative of the Hamiltonian, and is the recommended method.
        "numerical" uses numerical differentation. Defaults to
        "feynman-hellman"
    dk : float
        finite difference over which velocity is estimated. Only used
        if ``method == "feynman-hellman"``

    Returns:
    --------
    v : np.array((len(kFs), 2)) [i, n, (vx, vy)]
        list of velocities (vx, vy) at momenta kFs[i]
    """

    vs = np.zeros(shape=(len(kFs), 2), dtype=float)

    if method == "feynman-hellman":
        Es, psis = Es_psis(sys, kFs, indices)

        for n_b, start, stop in indices:
            for i in xrange(start, stop):
                kx = kFs[i, 0]
                ky = kFs[i, 1]
                psic = psis[i].conj()
                psip = psis[i]
                A = sys.dHdkx(kx, ky)
                B = sys.dHdky(kx, ky)
                vs[i, 0] = psic.dot(np.dot(A, psip)).real
                vs[i, 1] = psic.dot(np.dot(B, psip)).real

    elif method == "numerical":
        for n_b, start, stop in indices:
            for i in xrange(start, stop):
                kx = kFs[i, 0]
                ky = kFs[i, 1]

                E = eigvalsh(sys.H(kx, ky))[n_b]
                Edkx = eigvalsh(sys.H(kx + dk, ky))[n_b]
                Edky = eigvalsh(sys.H(kx, ky + dk))[n_b]

                vs[i, 0] = (Edkx-E)/dk
                vs[i, 1] = (Edky-E)/dk

    else:
        raise ValueError("Unknown method " + str(method))

    return vs


def deltar(sys, kFs, indices, dk=1e-5):
    """Quasi-particle side-jumps delta r_k'k

    Parameters:
    -----------
    sys : BoltzmannSystem
        system definition. Must implement the Hamiltonian ``H(kx, ky)``.
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    dk : float
        finite difference over which k-space derivatives are estimated.

    Returns:
    --------
    deltars : np.array((len(kFs), len(kFs), 2)) [i, j, (dr_x, dr_y)]
        side jumps (dr_x, dr_y) from momenta kFs[j] to kFs[i]
    """
    # deltar is obtained from gauge invariant Pancharatnam phase Phi_qkk'
    # see Sinitsyn et al. PRB 73, 075318 (2006)

    # states and derivatives
    psi = Es_psis(sys, kFs, indices)[1]
    dkx = np.zeros((len(kFs), 2))
    dkx[:, 0] = dk
    dky = np.zeros((len(kFs), 2))
    dky[:, 1] = dk
    psi_dkx = Es_psis(sys, kFs+dkx, indices)[1]
    psi_dky = Es_psis(sys, kFs+dky, indices)[1]

    # overlaps
    psic_dot_psi = np.einsum("ik,jk->ij", psi.conj(), psi)
    psic_dot_psi_dkx = np.einsum("ik,jk->ij", psi.conj(), psi_dkx)
    psic_dot_psi_dky = np.einsum("ik,jk->ij", psi.conj(), psi_dky)

    psic_dot_psi_dkx_row = np.diag(psic_dot_psi_dkx)
    psic_dot_psi_dkx_col = np.reshape(psic_dot_psi_dkx_row, (len(kFs), 1))

    psic_dot_psi_dky_row = np.diag(psic_dot_psi_dky)
    psic_dot_psi_dky_col = np.reshape(psic_dot_psi_dky_row, (len(kFs), 1))

    # discret derivatives of Phi_qkk' yield Phi_q+dk,kk'/dk in the
    # limit to k/k' we have:
    arg_q_to_k_x = np.angle(psic_dot_psi_dkx_col.conj() * psic_dot_psi
                            * psic_dot_psi_dkx.transpose())
    arg_q_to_kp_x = np.angle(psic_dot_psi_dkx.conj() * psic_dot_psi
                             * psic_dot_psi_dkx_row)

    arg_q_to_k_y = np.angle(psic_dot_psi_dky_col.conj() * psic_dot_psi
                            * psic_dot_psi_dky.transpose())
    arg_q_to_kp_y = np.angle(psic_dot_psi_dky.conj() * psic_dot_psi
                             * psic_dot_psi_dky_row)

    # delta r_k'k
    deltar = np.zeros((len(kFs), len(kFs), 2))
    deltar[:, :, 0] = 1/dk*(arg_q_to_k_x + arg_q_to_kp_x)
    deltar[:, :, 1] = 1/dk*(arg_q_to_k_y + arg_q_to_kp_y)

    return deltar


def amplitudes(sys, kFs, indices, dls, vFs, xi, nI=1/25**2, U=10, skew=0):
    """Scattering amplitudes

    Parameters:
    -----------
    sys : BoltzmannSystem
        system definition. Must implement the Hamiltonian ``H(kx, ky)``.
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    vFs : np.array((len(phis),)) [i]
        list of quasi-particle Fermi-velocities at momentum kFs[i]
    dls : np.array((len(phis),)) [i]
        list of Fermi-surface line elements at momentum kFs[i]
    xi : float
        correlation length of the scattering potential
        in units of the lattice constant
    nI : float
        impurity density nI*a^2
    U : float
        amplitude of the disorder-distribution U_i in U*[-1, 1+skew]
    skew : float
        skew-component of the disorder-distribution U_i in U*[-1, 1+skew]

    Returns:
    --------
    qs : np.array((len(kFs), len(kFs))) [j, i]
        scattering amplitudes from kFs[i], to kFs[j]
    """

    psis = Es_psis(sys, kFs, indices)[1]
    k_min_kp_sq = (np.subtract.outer(kFs[:, 0], kFs[:, 0])**2 +
                   np.subtract.outer(kFs[:, 1], kFs[:, 1])**2)
    V_kpk = np.exp(-xi**2/4 * k_min_kp_sq)
    psic_dot_psi = np.einsum("ik,jk->ij", psis.conj(), psis)

    # k-indpendent scattering amplitude to q0
    A0 = 2*pi**3/3 * U**2 * nI * xi**4 if xi > 0 else U**2 #* 10
    # 0th order symmetric amplitude
    q0 = A0 * V_kpk**2 * abs(psic_dot_psi)**2

    # k-independent scattering amplitude to q1
    A1 = -4*pi**5 * nI * U**3 * skew * xi**6
    # sum of intermediate momenta
    sum_A_q = np.einsum("iq,qj->ji",
                        dls/vFs/(2*pi)**2 * psic_dot_psi * V_kpk,
                        psic_dot_psi * V_kpk)
    # 1th order skew-symmetric amplitude
    q1 = A1 * V_kpk * (psic_dot_psi * sum_A_q).imag

    qs = q0 + q1
    np.fill_diagonal(qs, 0.0)
    return qs


def lambda_k(sys, kFs, indices, qs, vs, vFs, dls):#, deltars):
    """vector mean-free paths

    Parameters:
    -----------
    sys : BoltzmannSystem
        system definition. Must implement the Hamiltonian ``H(kx, ky)``.
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    qs : np.array((len(kFs), len(kFs))) [j, i]
        scattering amplitudes from kFs[i], to kFs[j]
    vs : np.array((len(kFs), 2)) [i, (vx, vy)]
        list of quasi-particle velocities (vx, vy) at momentum kFs[i]
    vFs : np.array((len(phis),)) [i]
        list of quasi-particle Fermi-velocities at momentum kFs[i]
    dls : np.array((len(phis),)) [i]
        list of Fermi-surface line elements at momentum kFs[i]
    deltars : np.array((len(kFs), len(kFs), 2)) [i, j, (dr_x, dr_y)]
        side jumps (dr_x, dr_y) from momenta kFs[j] to kFs[i]

    Returns:
    --------
    lambda_x : np.array((len(kFs),)) [i]
        vmfp in x-direction at momentum kFs[i]
    lambda_y : np.array((len(kFs),)) [i]
        vmfp in y-direction at momentum kFs[i]
    """

    # rate equation
    rates = np.zeros(shape=(len(kFs), len(kFs)), dtype=float)
    # off-diagonal elements: scattering-in rates
    for i in xrange(len(kFs)):
        rates[i, :] = -dls[i]/vFs[i]/(2*pi)**2 * qs[i]
    # diagonal elements: scattering-out rates
    for i in xrange(len(kFs)):
        rates[i, i] = np.sum(dls/vFs/(2*pi)**2 * qs[i])

    # side jump velocity
    #v_sj = np.einsum("ij,ijx->ix", dls/vFs/(2*pi)**2 * qs, deltars)

    # normalization condition is <lambda>_SF = 0
    # we chose the basis of the rate equation where < >_SF = sum(lambda)
    vx_tld = dls/vFs * (vs[:, 0])# + v_sj[:, 0])
    vy_tld = dls/vFs * (vs[:, 1])# + v_sj[:, 1])

    # orbital shift of quasiparticle distribution
    # e/hbar = mu_B/(t_e*a^2)
    # where t_e = hbar^2/(2*m_e*a^2) = 238 meV
    # mu_B*Tesla = 0.058 meV

    bz =  sys.p.B * sin(sys.p.alpha)# * 0.058/238
    vodl = vFs/dls

    orbital = np.zeros(rates.shape)
    #for band, start, stop in indices:
     #   orbital[[range(start, stop-1), range(start+1, stop)]] = \
      #      bz/2 * vodl[start+1: stop]
       # orbital[[range(start+1, stop), range(start, stop-1)]] = \
        #    -bz/2 * vodl[start: stop-1]
        #orbital[stop-1, start] = bz/2 * vodl[start]
        #orbital[start, stop-1] = -bz/2 * vodl[stop-1]

    # find a least square solution
    lambda_x, residuals_x, rank_x, s_x = np.linalg.lstsq(rates - orbital,
                                                         vx_tld, rcond=1e-8)
    lambda_y, residuals_y, rank_y, s_y = np.linalg.lstsq(rates - orbital,
                                                         vy_tld, rcond=1e-8)

    # normalize lambda (tensordot(rates, dls/vFs) = 0)
    lambda_x = lambda_x - np.sum(lambda_x)/np.sum(dls/vFs) * dls/vFs
    lambda_y = lambda_y - np.sum(lambda_y)/np.sum(dls/vFs) * dls/vFs

    return (lambda_x, lambda_y, rates, orbital, residuals_x, residuals_y,
            rank_x, rank_y, s_x, s_y)


def tau(sys, kFs, indices, vs, dls, xi, nI, U, skew):
    """Ziman relaxation time

    Parameters:
    -----------
    kFs : np.array(:, 2)
        array of Fermi momenta (kx, ky) for all contours.
    indices : np.array(Ncontour, 3)
        band index and indices in kFs for the Ncontour individual contours.
        `indices[i, 0]` is the band index of the i-th contour, and the
        momenta of the i-th contour are `kFs[index[i, 1]:index[i, 2], :]`.
    vs : np.array((len(kFs), 2)) [i, (vx, vy)]
        list of quasi-particle velocities (vx, vy) at momentum kFs[i]
    dls : np.array((len(phis),)) [i]
        list of Fermi-surface line elements at momentum kFs[i]
    xi : float
        correlation length of the scattering potential
        in units of the lattice constant
    nI : float
        impurity density nI*a^2
    U : float
        amplitude of the disorder-distribution U_i in U*[-1, 1+skew]
    skew : float
        skew-component of the disorder-distribution U_i in U*[-1, 1+skew]

    Retruns:
    --------
    taus : np.array((len(kFs),)) [i]
        Ziman relaxation time at momentu kFs[i]
    """

    vxs = vs[:, 0]
    vys = vs[:, 1]
    vFs = np.sqrt(vxs**2 + vys**2)

    amps = amplitudes(sys, kFs, indices, dls, vFs, xi, nI, U, skew)

    cos_v = np.zeros(shape=(len(kFs), len(kFs)), dtype=float)
    v_dot_vp = np.outer(vxs, vxs) + np.outer(vys, vys)
    vF_vFp = np.outer(vFs, vFs)
    cos_v[:, :] = 1 - v_dot_vp/vF_vFp

    taus = np.zeros(shape=(len(kFs),), dtype=float)
    for i in xrange(len(kFs)):
        taus[i] = 1/np.sum(dls/vFs * amps[i] * cos_v[i])

    return taus


def calc_sigma(sys, xi=0, nI=1/25**2, U=10, skew=0, E=0,
               kF_method="angle", kF_kwargs=None, v_kwargs=None,
               ziman=False):
    """Calculate conductivities sigma_xx, yy, xy

    Parameters:
    -----------
    sys : BoltzmannSystem
        system definition. Must implement the Hamiltonian ``H(kx, ky)``.
        Additionally, if v_kwargs is None, or v_kwargs specifies the
        method "feynman-hellman, also the derivatives of Hamiltonian with
        respect to kx and ky must be implemented: ``dHdkx(kx, ky)`` and
        ``dHdky(kx, ky)``, respectively.
    xi : float
        correlation length of the scattering potential
        in units of the lattice constant a
    nI : float
        impurity density nI*a^2
    U : float
        amplitude of the disorder-distribution U_i in U*[-1, 1+skew]
    skew : float
        skew-component of the disorder-distribution U_i in U*[-1, 1+skew]
    E : float
        energy
    kF_method : ("marchingsquares", "angle")
        method for computing the Fermi surface.
        "marchingsquares" - scan the Fermis surface using the marching squares
        algorithm (i.e. using a square grid, but with a smart way of
        sampling points only close to the contour). Recommended
        for a general Fermi surface.
        "angle" - scan the Fermi surface with equidistant steps in angles.
        Only works for Fermi surfaces that have a minimum at k=0.
    kF_kwargs : dictionary or None
        keyword arguments of kF_marchingsquares() or kF_angle() (depending
        on kF_method)
    v_kwargs : dictionary
        keyword arguments of v().
    ziman : bool
        whether to compute the conductivity in Ziman approximation

    Returns:
    --------
    (sigmaxx, sigmayy, sigmaxy, sigmayx) : tuple of floats
        conductivities
    """

    if kF_method == "marchingsquares":
        kFs, indices = kF_marchingsquares(sys, E=E, **kF_kwargs)
    elif kF_method == "angle":
        kFs, indices = kF_angle(sys, E=E, **kF_kwargs)
    else:
        raise ValueError("Unknown kF_method " + str(kF_method))

    if len(kFs) == 0:
        return 0., 0., 0., 0.

    vs = v(sys, kFs, indices, **v_kwargs)
    vxs = vs[:, 0]
    vys = vs[:, 1]
    vFs = np.sqrt(vxs**2 + vys**2)
    dls = dl(kFs, indices)

    if ziman:
        # use Ziman approxmation (does not include orbital part!)
        taus = tau(sys, kFs, indices, vs, dls, xi, nI, U, skew)
        sigma_xx = np.sum(dls/vFs * taus * vxs**2)
        sigma_yy = np.sum(dls/vFs * taus * vys**2)
        sigma_xy = np.sum(dls/vFs * taus * vxs*vys)

        return sigma_xx, sigma_yy, sigma_xy, sigma_xy
    # solve Bolzman equation numerically
    qs = amplitudes(sys, kFs, indices, dls, vFs, xi, nI, U, skew)
    # deltars seems small but calc is unstable, so it is safer
    # to ignore them for now!
    #deltars = np.zeros((len(kFs), len(kFs), 2))
    #deltars = deltar(sys, kFs, indices)
    lam_x, lam_y = lambda_k(sys, kFs, indices, qs, vs, vFs, dls)[:2]#, deltars)[:2]

    sigma_xx = np.sum(vxs * lam_x)/(2*pi)**2
    sigma_yy = np.sum(vys * lam_y)/(2*pi)**2
    sigma_xy = np.sum(vxs * lam_y)/(2*pi)**2
    sigma_yx = np.sum(vys * lam_x)/(2*pi)**2

    return sigma_xx, sigma_yy, sigma_xy, sigma_yx


def finite_T_sigma(T, DeltaE, N_E, sys, xi=0, nI=1/25**2, U=10, skew=0,
                   kF_method="marchingsquares", kF_kwargs=None, v_kwargs=None,
                   ziman=False):
    """calculate temperature averaged conductivties

    Parameters:
    -----------
    T : float
        temperature in Kelvin
    DeltaE : float
        limit of the energy window linspace(-DeltaE, DeltaE, N_E)
    N_E : int
        number of energies in the window
    others see ``calc_sigma``

    Returns:
    --------
    (sigmaxx, sigmayy, sigmaxy) : tuple of floats
        conductivities
    """

    def dn_FdE(E):
        return 1/(T*0.086173) * 1/(exp(E/(T*0.086173))
                                   + 2 + exp(-E/(T*0.086173)))

    sigmaxx, sigmayy, sigmaxy, sigmayx = (0, 0, 0, 0)

    for E in np.linspace(-DeltaE/2, +DeltaE/2, N_E):
        sxx, syy, sxy, syx = calc_sigma(sys, xi, nI, U, skew, E, kF_method,
                                        kF_kwargs, v_kwargs, ziman)
        sigmaxx += sxx * dn_FdE(E) * DeltaE/N_E
        sigmayy += syy * dn_FdE(E) * DeltaE/N_E
        sigmaxy += sxy * dn_FdE(E) * DeltaE/N_E
        sigmayx += syx * dn_FdE(E) * DeltaE/N_E

    return sigmaxx, sigmayy, sigmaxy, sigmayx


def rho(sxx_syy_sxy_syx):
    """resistivites (inverts conductivities)

    Paramters:
    ----------
    sxx_syy_sxy_syx : tuple of floats
        conductivities

    Returns:
    --------

    (rhoxx, rhoyy, rhoxy) : tuple of floats
        resitivities
    """

    sxx, syy, sxy, syx = sxx_syy_sxy_syx
    rhoxx = syy/(sxx*syy - sxy*syx)
    rhoyy = sxx/(sxx*syy - sxy*syx)
    rhoxy = -sxy/(sxx*syy - sxy*syx)
    rhoyx = -syx/(sxx*syy - sxy*syx)

    return rhoxx, rhoyy, rhoxy, rhoyx
