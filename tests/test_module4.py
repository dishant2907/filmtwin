"""Module 4 tests: statistics, DOE, SPC/OCAP and run-to-run control.

    python3 tests/test_module4.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from filmtwin import doe  # noqa: E402
from filmtwin.chamber import ChamberSimulator, run_scenario  # noqa: E402
from filmtwin.r2r import EWMAController, estimate_post_pm_boost, run_closed_loop, thickness_error_stats  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402
from filmtwin.spc import (ContextLimits, Limits, capability, count_events, ewma_chart,  # noqa: E402
                          western_electric)
from filmtwin.stats import f_sf, ols, t_ppf, t_sf2  # noqa: E402


# ------------------------------------------------------------------ stats
def test_t_and_f_against_known_values():
    assert abs(t_ppf(0.975, 10) - 2.228) < 1e-3  # textbook table value
    assert abs(t_sf2(2.228, 10) - 0.05) < 1e-3
    assert abs(f_sf(4.10, 2, 10) - 0.05) < 2e-3  # F(0.95; 2, 10) = 4.10


def test_ols_recovers_coefficients():
    rng = np.random.default_rng(0)
    X = np.column_stack([np.ones(200), rng.normal(size=200)])
    y = 3.0 + 2.0 * X[:, 1] + rng.normal(0, 0.1, 200)
    fit = ols(X, y, ["b0", "b1"])
    assert abs(fit.coef[0] - 3) < 0.05 and abs(fit.coef[1] - 2) < 0.05
    assert fit.r2 > 0.99 and fit.p[1] < 1e-10


# ------------------------------------------------------------------ DOE
def test_fractional_factorial_is_orthogonal_resolution_iv():
    d = doe.fractional_factorial_2_6_2(n_center=0)
    X = d[list(doe.FACTORS)].to_numpy()
    assert X.shape == (16, 6)
    assert np.allclose(X.T @ X, 16 * np.eye(6)), "main effects must be orthogonal"
    # Resolution IV: no main effect is aliased with any two-factor interaction
    for j in range(6):
        for a in range(6):
            for b in range(a + 1, 6):
                assert abs(X[:, j] @ (X[:, a] * X[:, b])) < 1e-9


def test_ccd_structure():
    d = doe.ccd_face_centered(["a", "b", "c"], n_center=6)
    assert len(d) == 8 + 6 + 6
    assert (d.point_type == "axial").sum() == 6


def test_screening_finds_known_physics():
    scr = doe.run_on_twin(doe.fractional_factorial_2_6_2(seed=1), seed=2)
    _, rate = doe.analyse_screening(scr, "rate_nm_min")
    _, nu = doe.analyse_screening(scr, "nu_pct")
    _, ri = doe.analyse_screening(scr, "ri")
    eff = lambda tab, f: float(tab.loc[tab.term == f, "effect"].iloc[0])
    assert eff(rate, "precursor_sccm") > 0, "more precursor -> faster deposition"
    assert eff(rate, "spacing_mm") < 0, "wider gap -> thicker boundary layer -> slower"
    assert eff(nu, "edge_temp_offset_C") < 0, "edge trim flattens the edge-thin profile"
    top_ri = ri[ri.term.isin(doe.FACTORS)].assign(a=lambda t: t.t.abs()).nlargest(1, "a").term.iloc[0]
    assert top_ri == "precursor_sccm", "flow ratio controls stoichiometry"


def test_profile_models_better_than_nu():
    """The signed profile is smooth through zero; 1-sigma NU folds into a V."""
    f = {"edge_temp_offset_C": (0.0, 30.0), "pressure_torr": (2.0, 3.0), "precursor_sccm": (120.0, 180.0)}
    rsm = doe.run_on_twin(doe.ccd_face_centered(list(f), seed=3), seed=4, factors=f)
    r2_profile = doe.fit_response_surface(rsm, "profile_pct", list(f)).r2
    r2_nu = doe.fit_response_surface(rsm, "nu_pct", list(f)).r2
    assert r2_profile > r2_nu


def test_desirability_bounds():
    D, _ = doe.desirability({"profile_pct": np.array([0.0, 5.0]), "rate_nm_min": np.array([340.0, 340.0]),
                             "ri": np.array([1.478, 1.478]), "werr_mean": np.array([1.3, 1.3])})
    assert np.isclose(D[0], 1.0) and D[1] == 0.0


# ------------------------------------------------------------------ SPC
def test_limits_and_capability():
    rng = np.random.default_rng(1)
    x = rng.normal(500, 2, 5000)
    lim = Limits.from_phase1(x)
    assert abs(lim.sigma - 2) < 0.1
    cap = capability(x, 485, 515)
    assert abs(cap["Cp"] - 2.5) < 0.15 and abs(cap["Cpk"] - 2.5) < 0.15


def test_we_rules_fire_on_shift():
    rng = np.random.default_rng(2)
    lim = Limits(0.0, 1.0, 1.128)
    x = np.concatenate([rng.normal(0, 1, 100), rng.normal(1.5, 1, 100)])
    we = western_electric(x, lim)
    assert we["any"][100:].sum() > 5 * max(we["any"][:100].sum(), 1)


def test_ewma_detects_small_shift():
    rng = np.random.default_rng(3)
    lim = Limits(0.0, 1.0, 1.128)
    x = np.concatenate([rng.normal(0, 1, 200), rng.normal(1.0, 1, 200)])
    ew = ewma_chart(x, lim, lam=0.2)
    assert ew.signal[:200].sum() < 5 and ew.signal[200:].mean() > 0.5


def test_count_events():
    assert count_events([0, 1, 1, 0, 1, 0, 0, 1, 1, 1]) == 3


# ------------------------------------------------------------------ R2R
def test_r2r_holds_thickness_through_mfc_drift():
    sc = SCENARIOS["mfc_drift"]
    open_w, _ = run_closed_loop(ChamberSimulator(faults=list(sc.faults), seed=7), 600, None)
    closed_w, _ = run_closed_loop(ChamberSimulator(faults=list(sc.faults), seed=7), 600, EWMAController(lam=0.6))
    o = thickness_error_stats(open_w, 500, 15)
    c = thickness_error_stats(closed_w, 500, 15)
    assert o["pct_out_of_spec"] > 20 and c["pct_out_of_spec"] < 1
    assert closed_w.r2r_time_s.iloc[-1] > 106, "controller should stretch time to compensate"


def test_feedforward_fixes_first_wafer():
    _, hist, _ = run_scenario(SCENARIOS["baseline"], seed=100)
    boost = estimate_post_pm_boost(hist)
    assert 0.005 < boost[0] < 0.03
    fb, _ = run_closed_loop(ChamberSimulator(seed=5), 600, EWMAController(lam=0.6))
    ff, _ = run_closed_loop(ChamberSimulator(seed=5), 600, EWMAController(lam=0.6, post_pm_boost=boost))
    assert thickness_error_stats(ff, 500, 15)["max_abs_error_nm"] < thickness_error_stats(fb, 500, 15)["max_abs_error_nm"]


def test_controller_adjustment_alarm():
    _, hist, _ = run_scenario(SCENARIOS["baseline"], seed=100)
    boost = estimate_post_pm_boost(hist)
    mk = lambda seed, faults=(): run_closed_loop(ChamberSimulator(faults=list(faults), seed=seed), 600,
                                                 EWMAController(lam=0.6, post_pm_boost=boost))[0]
    p1 = pd.concat([mk(6), mk(8)], ignore_index=True)
    ctx = ContextLimits.from_phase1(p1.r2r_time_s, p1.wafers_since_clean)
    healthy = mk(9)
    a_h, _ = ctx.adjustment_alarm(healthy.r2r_time_s, healthy.wafers_since_clean, 1.0)
    faulty = mk(7, SCENARIOS["mfc_drift"].faults)
    a_f, _ = ctx.adjustment_alarm(faulty.r2r_time_s, faulty.wafers_since_clean, 1.0)
    assert a_h.sum() == 0, "no false alarms on a healthy chamber"
    assert a_f[400:].all(), "sustained alarm once the drift is large"


if __name__ == "__main__":
    tests = [(n, f) for n, f in sorted(globals().items()) if n.startswith("test_") and callable(f)]
    failed = 0
    for name, fn in tests:
        try:
            fn()
            print(f"PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {name}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    sys.exit(1 if failed else 0)
