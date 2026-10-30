"""Film uniformity metrics as they are usually reported in a fab."""
from __future__ import annotations

import numpy as np


def uniformity(values) -> dict:
    """Mean, 1-sigma non-uniformity and range non-uniformity.

    nu_1sigma_pct  = std / mean * 100          (most common in thin films)
    nu_range_pct   = (max - min) / (2 * mean) * 100
    Different companies quote different definitions, so always state which
    one you used in a report.
    """
    v = np.asarray(values, dtype=float)
    mean = float(np.mean(v))
    std = float(np.std(v, ddof=1))
    return {
        "mean": mean,
        "std": std,
        "min": float(np.min(v)),
        "max": float(np.max(v)),
        "nu_1sigma_pct": std / mean * 100.0,
        "nu_range_pct": (float(np.max(v)) - float(np.min(v))) / (2.0 * mean) * 100.0,
    }


def radial_profile(xy_mm, values, n_bins: int = 4, r_max_mm: float = 147.0):
    """Average value per radial ring (centre to edge)."""
    r = np.sqrt(xy_mm[:, 0] ** 2 + xy_mm[:, 1] ** 2)
    edges = np.linspace(0, r_max_mm * 1.0001, n_bins + 1)
    centers, means = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (r >= a) & (r < b)
        if m.any():
            centers.append(float(np.mean(r[m])))
            means.append(float(np.mean(np.asarray(values)[m])))
    return np.array(centers), np.array(means)
