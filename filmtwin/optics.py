"""Thin-film optics: material models and transfer-matrix reflectance.

A spectral reflectometer shines broadband light on the wafer at normal
incidence and records reflectance R(lambda). Light reflected from the top of
the film interferes with light reflected from the film/substrate interface,
so R oscillates with wavelength. The spacing of those oscillations encodes
the optical thickness n*d, and their shape encodes n.

This module computes R(lambda) for any stack of layers on a substrate using
the characteristic (transfer) matrix method at normal incidence.

Silicon optical constants below are APPROXIMATE values read from standard
literature tables (e.g. Green, 2008). They are good enough for a simulator
where the "measurement" and the "model" use the same table. For real data,
replace them with a validated table.
"""
from __future__ import annotations

import numpy as np

# wavelength (nm), n, k for crystalline Si (approximate)
_SI_TABLE = np.array([
    [400, 5.57, 0.387],
    [450, 4.67, 0.130],
    [500, 4.29, 0.070],
    [550, 4.08, 0.041],
    [600, 3.94, 0.025],
    [633, 3.88, 0.018],
    [700, 3.79, 0.012],
    [800, 3.69, 0.006],
    [900, 3.63, 0.003],
    [1000, 3.57, 0.001],
])


def silicon_nk(wl_nm):
    """Complex refractive index N = n + ik of crystalline silicon."""
    wl = np.asarray(wl_nm, dtype=float)
    n = np.interp(wl, _SI_TABLE[:, 0], _SI_TABLE[:, 1])
    k = np.interp(wl, _SI_TABLE[:, 0], _SI_TABLE[:, 2])
    return n + 1j * k


# Cauchy model for transparent dielectrics: n(lambda) = A + B / lambda^2 (lambda in um)
CAUCHY_B_SIO2 = 0.0035  # um^2, typical for SiO2-like films


def cauchy_n(wl_nm, A, B=CAUCHY_B_SIO2):
    wl_um = np.asarray(wl_nm, dtype=float) / 1000.0
    return A + B / wl_um ** 2


def cauchy_A_from_n633(n633, B=CAUCHY_B_SIO2):
    """Cauchy A that gives refractive index n633 at 633 nm (the usual reporting wavelength)."""
    return n633 - B / 0.633 ** 2


def stack_reflectance(wl_nm, layers, n_substrate, n_ambient=1.0):
    """Normal-incidence reflectance of a multilayer stack.

    wl_nm       (W,) wavelengths
    layers      list of (N, d_nm) from TOP to BOTTOM. N is complex index,
                scalar or shape (W,). d_nm may be a scalar or an array of
                shape (G, 1) to evaluate G thicknesses at once (grid search).
    n_substrate complex index of the substrate, shape (W,)

    Returns R with shape (W,) or (G, W).
    """
    wl = np.asarray(wl_nm, dtype=float)
    m11 = np.ones(1, dtype=complex)
    m12 = np.zeros(1, dtype=complex)
    m21 = np.zeros(1, dtype=complex)
    m22 = np.ones(1, dtype=complex)
    for N, d in layers:
        N = np.asarray(N, dtype=complex)
        delta = 2.0 * np.pi * N * np.asarray(d, dtype=float) / wl
        c, s = np.cos(delta), np.sin(delta)
        a11, a12, a21, a22 = c, 1j * s / N, 1j * N * s, c
        m11, m12, m21, m22 = (m11 * a11 + m12 * a21, m11 * a12 + m12 * a22,
                              m21 * a11 + m22 * a21, m21 * a12 + m22 * a22)
    Ns = np.asarray(n_substrate, dtype=complex)
    B = m11 + m12 * Ns
    C = m21 + m22 * Ns
    r = (n_ambient * B - C) / (n_ambient * B + C)
    return np.abs(r) ** 2
