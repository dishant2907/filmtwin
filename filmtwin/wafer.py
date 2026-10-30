"""Wafer geometry: measurement site patterns and dense grids.

All positions are in millimetres, with the wafer centre at (0, 0).
The notch is assumed at the bottom of the wafer (theta = -90 degrees),
which is the usual convention for wafer map plots.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Wafer:
    """A round wafer with an edge exclusion zone."""

    diameter_mm: float = 300.0
    edge_exclusion_mm: float = 3.0

    @property
    def radius_mm(self) -> float:
        return self.diameter_mm / 2.0

    @property
    def usable_radius_mm(self) -> float:
        return self.radius_mm - self.edge_exclusion_mm


def sites_49pt(wafer: Wafer = Wafer()) -> np.ndarray:
    """A common 49-site polar pattern: 1 centre + rings of 8, 16 and 24 sites.

    Rings sit at 1/3, 2/3 and 3/3 of the usable radius. Returns an array of
    shape (49, 2) with x, y in mm.
    """
    r_max = wafer.usable_radius_mm
    pts = [(0.0, 0.0)]
    for n_sites, frac in ((8, 1 / 3), (16, 2 / 3), (24, 1.0)):
        r = r_max * frac
        # Small angular offset per ring so sites do not all line up radially.
        offset = np.pi / n_sites
        for k in range(n_sites):
            th = offset + 2 * np.pi * k / n_sites
            pts.append((r * np.cos(th), r * np.sin(th)))
    return np.array(pts)


def sites_line_scan(wafer: Wafer = Wafer(), n: int = 49, angle_deg: float = 0.0) -> np.ndarray:
    """A diameter line scan through the centre (useful for radial profiles)."""
    r_max = wafer.usable_radius_mm
    s = np.linspace(-r_max, r_max, n)
    a = np.deg2rad(angle_deg)
    return np.column_stack([s * np.cos(a), s * np.sin(a)])


def dense_grid(wafer: Wafer = Wafer(), n: int = 121):
    """Square grid masked to the usable wafer area (for contour plots).

    Returns (X, Y, mask) where mask is True inside the usable radius.
    """
    r_max = wafer.usable_radius_mm
    v = np.linspace(-r_max, r_max, n)
    X, Y = np.meshgrid(v, v)
    mask = X ** 2 + Y ** 2 <= r_max ** 2
    return X, Y, mask


def to_polar(xy: np.ndarray, wafer: Wafer = Wafer()):
    """Convert x, y (mm) to normalised radius rho in [0, 1] and angle theta (rad)."""
    x, y = xy[..., 0], xy[..., 1]
    rho = np.sqrt(x ** 2 + y ** 2) / wafer.usable_radius_mm
    theta = np.arctan2(y, x)
    return rho, theta
