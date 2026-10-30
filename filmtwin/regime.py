"""Reaction-limited vs transport-limited regime analysis (the Grove curve).

This is one of the most common interview topics for deposition roles:
"Your rate dropped 5 %. Is it temperature or gas delivery?" The answer
depends on which regime the process runs in, and this module shows it.
"""
from __future__ import annotations

import dataclasses

import numpy as np

from .physics import Condition, Hardware, Recipe, deposit


def thermal_cvd_hardware() -> Hardware:
    """A hypothetical thermal (non-plasma) CVD chamber with a strong activation energy.

    Used only to contrast with PECVD. With Ea ~150 kJ/mol the process is
    reaction-limited at low temperature and transport-limited at high
    temperature, giving the textbook Grove curve.
    """
    return Hardware(name="THERMAL", Ea_kJ_mol=150.0, k_s_ref=1.3, T_ref_K=673.15)


def arrhenius_sweep(hw: Hardware, recipe: Recipe, temps_C) -> dict:
    """Centre-point deposition rate and Damkohler number vs temperature."""
    centre = np.array([[0.0, 0.0]])
    rates, das = [], []
    for T in temps_C:
        r = deposit(centre, dataclasses.replace(recipe, temperature_C=float(T)), hw,
                    Condition(), noise=False)
        rates.append(float(r.rate_nm_min[0]))
        das.append(float(r.damkohler[0]))
    return {"T_C": np.asarray(temps_C, float), "rate_nm_min": np.array(rates), "Da": np.array(das)}


def local_sensitivities(hw: Hardware, recipe: Recipe, rel_step: float = 0.02) -> dict:
    """Percent change in mean rate for a small change in each recipe knob.

    Returns {knob: % rate change per +2 % (or +2 C for temperature)}. This is
    the quick 'which knob matters here' table an engineer builds before a DOE.
    """
    from .wafer import sites_49pt

    xy = sites_49pt()
    base = deposit(xy, recipe, hw, Condition(), noise=False).rate_nm_min.mean()
    out = {}
    for knob in ("temperature_C", "pressure_torr", "rf_power_W", "precursor_sccm",
                 "oxidizer_sccm", "spacing_mm"):
        v = getattr(recipe, knob)
        step = 2.0 if knob == "temperature_C" else v * rel_step
        r = deposit(xy, dataclasses.replace(recipe, **{knob: v + step}), hw, Condition(),
                    noise=False).rate_nm_min.mean()
        out[knob] = (r / base - 1.0) * 100.0
    return out
