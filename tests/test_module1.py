"""Module 1 tests. Run with either:

    python3 tests/test_module1.py
    python3 -m pytest tests          (if you have pytest)

Each test checks a piece of physics or simulator behaviour you should be
able to explain in an interview.
"""
import dataclasses
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from filmtwin.chamber import ChamberSimulator, run_scenario  # noqa: E402
from filmtwin.faults import FaultEvent  # noqa: E402
from filmtwin.metrics import radial_profile, uniformity  # noqa: E402
from filmtwin.physics import Condition, Hardware, Recipe, deposit  # noqa: E402
from filmtwin.regime import arrhenius_sweep, thermal_cvd_hardware  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402
from filmtwin.wafer import sites_49pt  # noqa: E402

XY = sites_49pt()
CENTRE = np.array([[0.0, 0.0]])


def mean_thk(recipe=Recipe(), hw=Hardware(), cond=Condition()):
    return deposit(XY, recipe, hw, cond, noise=False).thickness_nm.mean()


def test_49_sites():
    assert XY.shape == (49, 2)
    r = np.hypot(XY[:, 0], XY[:, 1])
    assert np.isclose(r.max(), 147.0)


def test_baseline_in_realistic_window():
    r = deposit(XY, Recipe(), noise=False)
    u = uniformity(r.thickness_nm)
    assert 150 < r.rate_nm_min.mean() < 600, "PECVD oxide rate should be a few hundred nm/min"
    assert 0.5 < u["nu_1sigma_pct"] < 2.5, "healthy chamber should be ~1-2 % 1-sigma"
    assert 0.2 < r.damkohler.mean() < 2.0, "default recipe sits in the mixed regime"


def test_grove_limits():
    """Reaction-limited at low T (steep Arrhenius), transport-limited at high T (flat)."""
    hw = thermal_cvd_hardware()
    sweep = arrhenius_sweep(hw, Recipe(), [300, 325, 600, 625])
    r = sweep["rate_nm_min"]
    low_T_gain = r[1] / r[0]
    high_T_gain = r[3] / r[2]
    assert sweep["Da"][0] < 0.1 and sweep["Da"][-1] > 10
    assert low_T_gain > 2.0, "reaction-limited rate should jump with temperature"
    assert high_T_gain < 1.1, "transport-limited rate should barely move"


def test_transport_limited_is_spacing_sensitive():
    """In the transport-limited regime the gap (boundary layer) matters more."""
    rec = Recipe(temperature_C=600)
    hw_t = thermal_cvd_hardware()
    wide = dataclasses.replace(rec, spacing_mm=12.0)
    drop_transport = 1 - mean_thk(wide, hw_t) / mean_thk(rec, hw_t)
    rec_low = Recipe(temperature_C=300)
    wide_low = dataclasses.replace(rec_low, spacing_mm=12.0)
    drop_reaction = 1 - mean_thk(wide_low, hw_t) / mean_thk(rec_low, hw_t)
    assert drop_transport > 5 * drop_reaction


def test_knob_directions():
    base = mean_thk()
    assert mean_thk(dataclasses.replace(Recipe(), precursor_sccm=200)) > base
    assert mean_thk(dataclasses.replace(Recipe(), oxidizer_sccm=4000)) < base  # dilution
    assert mean_thk(dataclasses.replace(Recipe(), rf_power_W=800)) > base
    assert mean_thk(dataclasses.replace(Recipe(), pressure_torr=3.5)) > base
    assert mean_thk(dataclasses.replace(Recipe(), dep_time_s=206)) > 1.99 * base


def test_fault_signatures():
    base = deposit(XY, Recipe(), noise=False).thickness_nm
    _, b = radial_profile(XY, base)

    clog = deposit(XY, Recipe(), cond=Condition(center_gas_loss=0.12), noise=False).thickness_nm
    _, c = radial_profile(XY, clog)
    assert (c[0] - b[0]) < -0.05 * b[0], "clog should thin the centre"
    assert abs(c[-1] - b[-1]) < 0.01 * b[-1], "clog should leave the edge alone"

    ring = deposit(XY, Recipe(), cond=Condition(extra_edge_plasma=0.10), noise=False).thickness_nm
    _, e = radial_profile(XY, ring)
    assert e[-1] > b[-1] * 1.03, "edge ring erosion should thicken the edge"

    tilt = deposit(XY, Recipe(), cond=Condition(extra_tilt=0.015, extra_tilt_angle_deg=0), noise=False).thickness_nm
    right = XY[:, 0] > 100
    left = XY[:, 0] < -100
    assert tilt[right].mean() > tilt[left].mean()

    rf = deposit(XY, Recipe(), cond=Condition(rf_delivered_factor=0.9), noise=False)
    mfc = deposit(XY, Recipe(), cond=Condition(precursor_flow_error=-0.05), noise=False)
    assert rf.thickness_nm.mean() < base.mean() and mfc.thickness_nm.mean() < base.mean()
    # Same thickness symptom, different refractive index: the key to telling them apart.
    assert mfc.refractive_index.mean() < rf.refractive_index.mean() - 0.001


def test_edge_tc_drift_shows_in_heater_power():
    _, w, _ = run_scenario(SCENARIOS["edge_tc_drift"], seed=3)
    early = w.iloc[:100]["s_heater_edge_pct"].mean()
    late = w.iloc[-50:]["s_heater_edge_pct"].mean()
    assert late < early - 8
    # ...while the centre zone stays put
    assert abs(w.iloc[-50:]["s_heater_center_pct"].mean() - w.iloc[:100]["s_heater_center_pct"].mean()) < 0.5


def test_rf_match_shows_in_reflected_power():
    _, w, _ = run_scenario(SCENARIOS["rf_match_degradation"], seed=4)
    assert w.iloc[-50:]["s_rf_reflected_W"].mean() > 30
    assert w.iloc[:150]["s_rf_reflected_W"].mean() < 5


def test_pm_resets_wall_aging():
    _, w, _ = run_scenario(SCENARIOS["baseline"], seed=5)
    before = w[(w.wafers_since_clean > 250)]["s_vpp_V"].mean()
    after = w[(w.wafers_since_clean.between(5, 30))]["s_vpp_V"].mean()
    assert before > after
    first = w[w.wafers_since_clean == 0]["m_thickness_mean_nm"].mean()
    settled = w[w.wafers_since_clean.between(5, 30)]["m_thickness_mean_nm"].mean()
    assert first < settled, "first wafer after a clean should run low"


def test_reproducible_with_seed():
    a = ChamberSimulator(seed=7).run(30)[0]["m_thickness_mean_nm"].values
    b = ChamberSimulator(seed=7).run(30)[0]["m_thickness_mean_nm"].values
    assert np.allclose(a, b)


def test_ground_truth_labels():
    sim = ChamberSimulator(faults=[FaultEvent("mfc_drift", start_wafer=10, severity=-0.03)], seed=1)
    w, s = sim.run(20)
    assert (w.loc[w.wafer_index < 10, "gt_label"] == "none").all()
    assert (w.loc[w.wafer_index >= 10, "gt_label"] == "mfc_drift").all()
    assert len(s) == 20 * 49


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
