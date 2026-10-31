"""Turn raw wafer data into a compact feature vector per wafer.

A process engineer does not stare at 49 numbers per wafer. They summarise
each wafer by a few physically meaningful numbers: how thick, how tilted,
how bowl-shaped, what's happening at the centre and edge, how the film
quality looks, and what the tool sensors said. This module builds that
summary for every wafer; Module 3's diagnosis works on it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .wafer import Wafer, sites_49pt
from .zernike import ZernikeFitter, normalise_map

# Site index groups for the 49-site pattern (centre, then rings of 8, 16, 24)
CENTRE = np.array([0])
RING1 = np.arange(1, 9)
RING2 = np.arange(9, 25)
RING3 = np.arange(25, 49)

# Radial information on a 49-site map has only 3 degrees of freedom (4 radii,
# minus the mean). defocus + centre_delta + edge_delta already use all three,
# so the Zernike 'spherical' term would be an exact linear combination of them
# and make the healthy covariance matrix singular. It is left out on purpose.
ZERNIKE_FEATURES = ["tilt_x", "tilt_y", "defocus", "astig_0", "astig_45"]
MAP_FEATURES = ["mean_nm"] + ZERNIKE_FEATURES + ["centre_delta_pct", "edge_delta_pct"]
FILM_FEATURES = ["ri_mean", "werr_edge_delta"]
SENSOR_FEATURES = ["s_rf_reflected_W", "s_vpp_V", "s_heater_edge_pct", "s_heater_center_pct",
                   "s_throttle_valve_pct"]
FEATURES = MAP_FEATURES + FILM_FEATURES + SENSOR_FEATURES

_FITTER = None


def _fitter():
    global _FITTER
    if _FITTER is None:
        _FITTER = ZernikeFitter(sites_49pt(Wafer()))
    return _FITTER


def features_from_arrays(thk, ri, werr, sensors: dict = None) -> pd.DataFrame:
    """thk, ri, werr: (n_wafers, 49) arrays in the standard site order."""
    thk = np.atleast_2d(thk)
    ri = np.atleast_2d(ri)
    werr = np.atleast_2d(werr)
    mean = thk.mean(axis=1)
    coeffs = _fitter().fit(normalise_map(thk))
    names = _fitter().names
    out = {"mean_nm": mean}
    for f in ZERNIKE_FEATURES:
        out[f] = coeffs[:, names.index(f)]
    out["centre_delta_pct"] = (thk[:, CENTRE].mean(axis=1) - thk[:, RING1].mean(axis=1)) / mean * 100
    out["edge_delta_pct"] = (thk[:, RING3].mean(axis=1) - thk[:, RING2].mean(axis=1)) / mean * 100
    out["ri_mean"] = ri.mean(axis=1)
    inner = np.concatenate([CENTRE, RING1])
    out["werr_edge_delta"] = werr[:, RING3].mean(axis=1) - werr[:, inner].mean(axis=1)
    if sensors:
        for k in SENSOR_FEATURES:
            out[k] = np.atleast_1d(sensors[k])
    return pd.DataFrame(out)


def wafer_features(wafers: pd.DataFrame, sites: pd.DataFrame) -> pd.DataFrame:
    """Feature table for a simulated (or real) run: one row per wafer.

    Keeps context columns (wafers_since_clean etc.) and any ground-truth
    `gt_*` columns so methods can be scored.
    """
    s = sites.sort_values(["wafer_index", "site"])
    n_sites = s["site"].nunique()
    if n_sites != 49:
        raise ValueError("feature extraction expects the standard 49-site pattern")
    thk = s["thickness_nm"].to_numpy().reshape(-1, 49)
    ri = s["refractive_index"].to_numpy().reshape(-1, 49)
    werr = s["werr"].to_numpy().reshape(-1, 49)
    w = wafers.sort_values("wafer_index").reset_index(drop=True)
    sensors = {k: w[k].to_numpy() for k in SENSOR_FEATURES}
    feats = features_from_arrays(thk, ri, werr, sensors)
    keep = ["wafer_index", "wafer_id", "lot_id", "wafers_since_clean", "pm_event",
            "m_nu_1sigma_pct"] + [c for c in w.columns if c.startswith("gt_")]
    if "scenario" in w.columns:
        keep = ["scenario"] + keep
    return pd.concat([w[keep], feats], axis=1)
