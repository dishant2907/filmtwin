"""Zernike decomposition of wafer maps.

Zernike polynomials are an orthogonal set of shapes on a disk. Any wafer map
can be written as a weighted sum of them, and each weight tells you how much
of one basic shape is in the map:

    piston      overall mean (removed, since we fit the normalised map)
    tilt x/y    one side thick, the other thin       -> leveling, pumping asymmetry
    defocus     bowl or dome                         -> centre-to-edge balance
    astigmatism saddle                               -> two-fold asymmetry
    coma        lopsided bowl                        -> off-centre hardware
    spherical   ring / W shape                       -> sharp edge or centre features

Coefficients use OSA/ANSI indexing and RMS normalisation, so each coefficient
equals that shape's RMS contribution in the units of the map (here, percent
of mean thickness). That makes them directly comparable.
"""
from __future__ import annotations

from math import factorial

import numpy as np

# (n, m, name) up to radial order 4
TERMS = [
    (0, 0, "piston"),
    (1, -1, "tilt_y"),
    (1, 1, "tilt_x"),
    (2, -2, "astig_45"),
    (2, 0, "defocus"),
    (2, 2, "astig_0"),
    (3, -3, "trefoil_y"),
    (3, -1, "coma_y"),
    (3, 1, "coma_x"),
    (3, 3, "trefoil_x"),
    (4, -4, "quad_y"),
    (4, -2, "astig2_45"),
    (4, 0, "spherical"),
    (4, 2, "astig2_0"),
    (4, 4, "quad_x"),
]
NAMES = [t[2] for t in TERMS]


def radial(n: int, m: int, rho):
    m = abs(m)
    out = np.zeros_like(rho, dtype=float)
    for k in range((n - m) // 2 + 1):
        c = ((-1) ** k * factorial(n - k)
             / (factorial(k) * factorial((n + m) // 2 - k) * factorial((n - m) // 2 - k)))
        out += c * rho ** (n - 2 * k)
    return out


def zernike(n: int, m: int, rho, theta):
    """RMS-normalised Zernike polynomial."""
    norm = np.sqrt(2 * (n + 1) / (1 + (m == 0)))
    R = radial(n, m, rho)
    if m >= 0:
        return norm * R * np.cos(m * theta)
    return norm * R * np.sin(-m * theta)


def basis(xy_mm, r_max_mm: float = 147.0, terms=TERMS):
    """Design matrix (n_sites, n_terms) for the given site positions."""
    rho = np.hypot(xy_mm[:, 0], xy_mm[:, 1]) / r_max_mm
    theta = np.arctan2(xy_mm[:, 1], xy_mm[:, 0])
    return np.column_stack([zernike(n, m, rho, theta) for n, m, _ in terms])


class ZernikeFitter:
    """Fits many wafer maps on the same site pattern with one matrix multiply."""

    def __init__(self, xy_mm, r_max_mm: float = 147.0, terms=TERMS):
        self.terms = terms
        self.names = [t[2] for t in terms]
        self.B = basis(xy_mm, r_max_mm, terms)
        self.pinv = np.linalg.pinv(self.B)

    def fit(self, maps):
        """maps: (n_wafers, n_sites) or (n_sites,). Returns coefficients."""
        maps = np.atleast_2d(maps)
        return maps @ self.pinv.T

    def reconstruct(self, coeffs, xy_mm=None, r_max_mm: float = 147.0):
        B = self.B if xy_mm is None else basis(xy_mm, r_max_mm, self.terms)
        return np.atleast_2d(coeffs) @ B.T


def normalise_map(thickness_nm):
    """Percent deviation from each wafer's own mean: (t / mean - 1) * 100."""
    t = np.atleast_2d(thickness_nm)
    return (t / t.mean(axis=1, keepdims=True) - 1.0) * 100.0
