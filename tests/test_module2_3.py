"""Module 2 (metrology) and Module 3 (diagnosis) tests.

    python3 tests/test_module2_3.py
"""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from filmtwin.chamber import run_scenario  # noqa: E402
from filmtwin.diagnosis import Diagnoser, chi2_quantile, nnls, score  # noqa: E402
from filmtwin.features import wafer_features  # noqa: E402
from filmtwin.metrology import Reflectometer, gauge_study  # noqa: E402
from filmtwin.optics import cauchy_n, silicon_nk, stack_reflectance  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402
from filmtwin.wafer import sites_49pt  # noqa: E402
from filmtwin.zernike import ZernikeFitter, basis, normalise_map  # noqa: E402

# ---------------------------------------------------------------- Module 2
def test_bare_silicon_reflectance():
    """Bare Si at 633 nm reflects ~35 % (textbook value)."""
    wl = np.array([633.0])
    R = stack_reflectance(wl, [], silicon_nk(wl))
    assert 0.33 < R[0] < 0.37


def test_quarter_wave_minimum():
    """Quarter-wave film (n*d = lambda/4) on a non-dispersive substrate: the minimum sits
    exactly at lambda with R = ((n0*ns - n1^2) / (n0*ns + n1^2))^2. On real Si the minimum
    shifts a little because Si's index changes with wavelength."""
    wl = np.linspace(400, 1000, 6001)
    n1, ns = 1.46, 3.9
    d = 600.0 / (4 * n1)
    R = stack_reflectance(wl, [(n1 + 0j, d)], np.full(wl.shape, ns + 0j))
    assert abs(wl[np.argmin(R)] - 600.0) < 0.2
    assert np.isclose(R.min(), ((ns - n1 ** 2) / (ns + n1 ** 2)) ** 2, rtol=1e-6)


def test_fit_recovers_truth():
    m = Reflectometer(seed=3)
    f = m.fit(m.measure(500.0, 1.478))
    assert f.ok
    assert abs(f.thickness_nm - 500.0) < 1.0
    assert abs(f.n633 - 1.478) < 0.003


def test_thin_film_falls_back_to_fixed_index():
    m = Reflectometer(seed=4)
    f = m.fit(m.measure(20.0, 1.478))
    assert "index fixed" in f.message
    assert abs(f.thickness_nm - 20.0) < 2.0


def test_bad_start_is_caught_by_gof():
    """Skipping the global search can land in the wrong valley; GOF must flag it."""
    m = Reflectometer(seed=5)
    f = m.fit(m.measure(500.0, 1.478), initial_thickness=380, fit_index=False)
    assert not f.ok


def test_gauge_capability():
    g = gauge_study(Reflectometer(seed=6), 500.0, 1.478, repeats=15, tolerance_nm=30.0)
    assert g["repeatability_nm"] < 0.5
    assert g["p_to_t"] < 0.1, "gauge should be capable for a +/-3 % spec"


def test_wrong_stack_model_flagged():
    m = Reflectometer(seed=7)
    sin = (cauchy_n(m.wl, 1.95, B=0.018).astype(complex), 80.0)
    f = m.fit(m.measure(500.0, 1.478, under_layers=[sin]))
    assert not f.ok, "a single-layer model on a two-layer stack should fail the GOF check"


# ---------------------------------------------------------------- Module 3
def test_zernike_recovers_known_shape():
    xy = sites_49pt()
    zf = ZernikeFitter(xy)
    true = np.zeros(len(zf.names))
    true[zf.names.index("tilt_x")] = 0.8
    true[zf.names.index("defocus")] = -0.5
    maps = (true @ basis(xy).T)[None, :]
    est = zf.fit(maps)[0]
    assert np.allclose(est, true, atol=1e-8)


def test_normalise_map_is_percent():
    m = normalise_map(np.array([[99.0, 101.0]]))
    assert np.allclose(m, [[-1.0, 1.0]])


def test_nnls_matches_known_solution():
    A = np.array([[1.0, 0.0], [0.0, 1.0], [1.0, 1.0]])
    x = nnls(A, np.array([2.0, -1.0, 1.0]))
    assert x[1] == 0.0 and x[0] > 0


def test_chi2_quantile_approximation():
    assert abs(chi2_quantile(0.95, 10) - 18.307) < 0.1
    assert abs(chi2_quantile(0.999, 15) - 37.697) < 0.5


_DIAG = None


def _diagnoser():
    global _DIAG
    if _DIAG is None:
        train = pd.concat([wafer_features(*run_scenario(SCENARIOS["baseline"], seed=s)[1:])
                           for s in (100, 101)], ignore_index=True)
        _DIAG = Diagnoser.train(train)
    return _DIAG


def test_no_false_alarms_on_fresh_healthy_run():
    f = wafer_features(*run_scenario(SCENARIOS["baseline"], seed=999)[1:])
    r = _diagnoser().run(f)
    assert r["alarm"].mean() < 0.01


def test_each_fault_detected_and_diagnosed():
    for i, name in enumerate(["edge_tc_drift", "showerhead_clog", "rf_match_degradation", "leveling_tilt"]):
        sc = SCENARIOS[name]
        f = wafer_features(*run_scenario(sc, seed=700 + i)[1:])
        r = _diagnoser().run(f)
        s = score(r, f.gt_label, sc.faults[0].start_wafer)
        assert s["detected"], name
        assert s["diagnosis_accuracy"] > 0.9, (name, s)


def test_mixed_faults_both_found():
    sc = SCENARIOS["mixed_clog_and_match"]
    f = wafer_features(*run_scenario(sc, seed=801)[1:])
    r = _diagnoser().run(f)
    assert (r.tail(20).diagnosis == "rf_match_degradation+showerhead_clog").mean() > 0.8


def test_severity_estimate_close_to_truth():
    sc = SCENARIOS["rf_match_degradation"]
    f = wafer_features(*run_scenario(sc, seed=802)[1:])
    r = _diagnoser().run(f)
    est = r.tail(30)["est_rf_match_degradation"].mean()
    assert abs(est - 0.10) < 0.01


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
