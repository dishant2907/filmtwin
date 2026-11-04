"""Module 4a: Design of Experiments on the chamber twin.

The standard process-development workflow, end to end:

1. SCREEN    A two-level fractional factorial (2^(6-2), resolution IV, 16 runs
             plus centre points) tests six recipe knobs at once and finds the
             few that matter. Centre points give pure error and a curvature check.
2. OPTIMISE  A face-centred central composite design (CCD) on the vital few
             fits a full quadratic response surface for each response.
3. CHOOSE    A desirability function trades the responses off against each
             other (rate vs uniformity vs film quality) and finds the best
             recipe inside the experimental region: the process window.
4. CONFIRM   Run the chosen recipe and check that the results fall inside the
             model's prediction intervals. If they don't, the model is wrong.

Every run is a real simulated wafer with run-to-run noise, in randomised order.
"""
from __future__ import annotations

import dataclasses
import itertools

import numpy as np
import pandas as pd

from .metrics import uniformity
from .physics import Condition, Hardware, Recipe, deposit
from .stats import ols
from .wafer import sites_49pt

# Factor ranges (low, high). Chosen around the default recipe, inside a
# plausible PECVD oxide window.
FACTORS = {
    "temperature_C": (350.0, 450.0),
    "pressure_torr": (2.0, 3.0),
    "rf_power_W": (500.0, 700.0),
    "precursor_sccm": (120.0, 180.0),
    "spacing_mm": (8.0, 12.0),
    "edge_temp_offset_C": (0.0, 15.0),
}
RESPONSES = ["rate_nm_min", "nu_pct", "profile_pct", "ri", "werr_mean", "werr_edge_delta"]
# profile_pct: signed edge-minus-centre thickness (% of mean). Unlike 1-sigma NU,
# which is a standard deviation and folds into a V at its minimum, the signed
# profile changes smoothly through zero, so a quadratic model fits it well.


# ---------------------------------------------------------------------------
# Coding
# ---------------------------------------------------------------------------
def to_real(coded: pd.DataFrame, factors=FACTORS) -> pd.DataFrame:
    out = pd.DataFrame(index=coded.index)
    for f, (lo, hi) in factors.items():
        if f in coded:
            out[f] = (lo + hi) / 2 + coded[f] * (hi - lo) / 2
    return out


def to_coded(real: pd.DataFrame, factors=FACTORS) -> pd.DataFrame:
    out = pd.DataFrame(index=real.index)
    for f, (lo, hi) in factors.items():
        if f in real:
            out[f] = (real[f] - (lo + hi) / 2) / ((hi - lo) / 2)
    return out


# ---------------------------------------------------------------------------
# Designs
# ---------------------------------------------------------------------------
def fractional_factorial_2_6_2(n_center: int = 4, seed: int = 0) -> pd.DataFrame:
    """2^(6-2) resolution IV design: E = ABC, F = BCD.

    Main effects are not aliased with two-factor interactions (resolution IV);
    two-factor interactions are aliased with each other.
    """
    names = list(FACTORS)
    base = np.array(list(itertools.product([-1, 1], repeat=4)))[:, ::-1]  # standard order
    A, B, C, D = base.T
    E = A * B * C
    F = B * C * D
    X = np.column_stack([A, B, C, D, E, F]).astype(float)
    X = np.vstack([X, np.zeros((n_center, 6))])
    d = pd.DataFrame(X, columns=names)
    d.insert(0, "std_order", np.arange(1, len(d) + 1))
    d["center_point"] = (d[names] == 0).all(axis=1)
    rng = np.random.default_rng(seed)
    d["run_order"] = rng.permutation(len(d)) + 1
    return d.sort_values("run_order").reset_index(drop=True)


def ccd_face_centered(factors: list, n_center: int = 6, seed: int = 0) -> pd.DataFrame:
    """Face-centred central composite design (alpha = 1) for a quadratic model."""
    k = len(factors)
    cube = np.array(list(itertools.product([-1, 1], repeat=k)), dtype=float)
    axial = np.vstack([np.eye(k), -np.eye(k)])
    center = np.zeros((n_center, k))
    X = np.vstack([cube, axial, center])
    d = pd.DataFrame(X, columns=factors)
    d.insert(0, "std_order", np.arange(1, len(d) + 1))
    d["point_type"] = ["cube"] * len(cube) + ["axial"] * len(axial) + ["center"] * n_center
    rng = np.random.default_rng(seed)
    d["run_order"] = rng.permutation(len(d)) + 1
    return d.sort_values("run_order").reset_index(drop=True)


# ---------------------------------------------------------------------------
# Running the experiment on the twin
# ---------------------------------------------------------------------------
def run_on_twin(design: pd.DataFrame, base_recipe: Recipe = None, hw: Hardware = None,
                cond: Condition = None, seed: int = 0, factors=FACTORS) -> pd.DataFrame:
    """Process one wafer per design row (in run order) and record the responses."""
    base_recipe = base_recipe or Recipe()
    hw = hw or Hardware()
    cond = cond or Condition()
    rng = np.random.default_rng(seed)
    xy = sites_49pt()
    names = [f for f in factors if f in design.columns]
    real = to_real(design[names], factors)
    rows = []
    for idx in design.index:
        recipe = dataclasses.replace(base_recipe, **{f: float(real.loc[idx, f]) for f in names})
        r = deposit(xy, recipe, hw, cond, rng)
        u = uniformity(r.thickness_nm)
        outer = r.werr[25:].mean()
        inner = r.werr[:9].mean()
        t = r.thickness_nm
        profile = (t[25:].mean() - t[:9].mean()) / t.mean() * 100.0
        rows.append({
            "rate_nm_min": u["mean"] / recipe.dep_time_s * 60.0,
            "nu_pct": u["nu_1sigma_pct"],
            "profile_pct": float(profile),
            "ri": float(np.mean(r.refractive_index)),
            "werr_mean": float(np.mean(r.werr)),
            "werr_edge_delta": float(outer - inner),
        })
    out = pd.concat([design.reset_index(drop=True), real.add_suffix("_real").reset_index(drop=True),
                     pd.DataFrame(rows)], axis=1)
    return out


# ---------------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------------
def model_matrix(coded: pd.DataFrame, factors: list, kind: str = "linear", curvature: bool = False):
    """Columns: intercept, mains, [2FI], [squares]. Returns (X, term names)."""
    cols, names = [np.ones(len(coded))], ["intercept"]
    for f in factors:
        cols.append(coded[f].to_numpy(float))
        names.append(f)
    if kind in ("interaction", "quadratic"):
        for a, b in itertools.combinations(factors, 2):
            cols.append(coded[a].to_numpy(float) * coded[b].to_numpy(float))
            names.append(f"{a}*{b}")
    if kind == "quadratic":
        for f in factors:
            cols.append(coded[f].to_numpy(float) ** 2)
            names.append(f"{f}^2")
    if curvature:
        cp = (coded[factors] == 0).all(axis=1).to_numpy(float)
        cols.append(cp)
        names.append("curvature(center)")
    return np.column_stack(cols), names


def analyse_screening(results: pd.DataFrame, response: str, factors=None):
    """Main-effects model with a centre-point curvature check.

    Effect = 2 x coefficient (the change in response going from low to high).
    """
    factors = factors or list(FACTORS)
    X, names = model_matrix(results, factors, "linear", curvature=True)
    fit = ols(X, results[response], names)
    tab = fit.table()
    tab["effect"] = np.where(tab.term.isin(factors), 2 * tab.estimate, np.nan)
    return fit, tab


def fit_response_surface(results: pd.DataFrame, response: str, factors: list):
    X, names = model_matrix(results, factors, "quadratic")
    return ols(X, results[response], names)


# ---------------------------------------------------------------------------
# Desirability (Derringer-Suich)
# ---------------------------------------------------------------------------
def d_smaller(y, target, upper):
    return np.clip((upper - y) / (upper - target), 0, 1)


def d_larger(y, lower, target):
    return np.clip((y - lower) / (target - lower), 0, 1)


def d_target(y, lower, target, upper):
    return np.where(y <= target, np.clip((y - lower) / (target - lower), 0, 1),
                    np.clip((upper - y) / (upper - target), 0, 1))


DEFAULT_GOALS = {
    # response: (kind, params)
    "profile_pct": ("target", (-3.0, 0.0, 3.0)),  # flat centre-to-edge profile
    "rate_nm_min": ("larger", (250.0, 340.0)),  # throughput
    "ri": ("target", (1.470, 1.478, 1.486)),  # stoichiometry window
    "werr_mean": ("smaller", (1.3, 1.9)),  # film density
}


def desirability(preds: dict, goals=DEFAULT_GOALS):
    """Overall desirability = geometric mean of the individual desirabilities."""
    ds = []
    for resp, (kind, p) in goals.items():
        y = preds[resp]
        if kind == "smaller":
            ds.append(d_smaller(y, *p))
        elif kind == "larger":
            ds.append(d_larger(y, *p))
        else:
            ds.append(d_target(y, *p))
    ds = np.array(ds)
    return np.prod(ds, axis=0) ** (1.0 / len(ds)), ds


def optimise(models: dict, factors: list, goals=DEFAULT_GOALS, n_grid: int = 41):
    """Grid search over the coded experimental region for maximum desirability."""
    g = np.linspace(-1, 1, n_grid)
    mesh = np.array(list(itertools.product(g, repeat=len(factors))))
    coded = pd.DataFrame(mesh, columns=factors)
    X, _ = model_matrix(coded, factors, "quadratic")
    preds = {r: X @ m.coef for r, m in models.items()}
    D, _ = desirability(preds, goals)
    best = int(np.argmax(D))
    return coded.iloc[best], {r: float(p[best]) for r, p in preds.items()}, float(D[best])
