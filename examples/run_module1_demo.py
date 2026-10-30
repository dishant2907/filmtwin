"""Module 1 demo: run the chamber twin and generate figures + datasets.

Run from the project root:

    python3 examples/run_module1_demo.py            # all scenarios
    python3 examples/run_module1_demo.py --quick    # skip the CSV export

Outputs:
    outputs/figures/*.png   maps, fault gallery, regime plot, run charts
    data/*.csv              wafer-level and site-level tables per scenario
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from filmtwin.chamber import run_scenario  # noqa: E402
from filmtwin.faults import Aging, build_condition  # noqa: E402
from filmtwin.physics import Hardware, Recipe  # noqa: E402
from filmtwin.plots import plot_baseline, plot_fault_gallery, plot_regime, plot_trend  # noqa: E402
from filmtwin.regime import arrhenius_sweep, local_sensitivities, thermal_cvd_hardware  # noqa: E402
from filmtwin.scenarios import SCENARIOS  # noqa: E402

ROOT = os.path.join(os.path.dirname(__file__), "..")
FIG = os.path.join(ROOT, "outputs", "figures")
DATA = os.path.join(ROOT, "data")


def main(quick: bool = False):
    os.makedirs(FIG, exist_ok=True)
    os.makedirs(DATA, exist_ok=True)
    hw, recipe = Hardware(), Recipe()

    # 1. Baseline maps
    plot_baseline(recipe, hw, os.path.join(FIG, "01_baseline_maps.png"))
    print("[1/4] Baseline maps saved")

    # 2. Fault gallery: each single fault at full severity, aging switched off
    no_aging = Aging(wall_loss=0.0, first_wafer_dip=0.0)
    conditions = {}
    for name, sc in SCENARIOS.items():
        if len(sc.faults) != 1:
            continue
        cond, _ = build_condition(sc.faults, wafer_index=10_000, wafers_since_clean=1000,
                                  last_pm_wafer=0, aging=no_aging)
        conditions[name] = cond
    plot_fault_gallery(recipe, hw, conditions, os.path.join(FIG, "02_fault_gallery.png"))
    print("[2/4] Fault gallery saved")

    # 3. Regime analysis
    temps = np.arange(250, 651, 25)
    pecvd = arrhenius_sweep(hw, recipe, temps)
    thermal = arrhenius_sweep(thermal_cvd_hardware(), recipe, temps)
    plot_regime(pecvd, thermal, os.path.join(FIG, "03_regime_grove_curve.png"))
    print("[3/4] Regime plot saved. Rate sensitivity at the default recipe:")
    for knob, pct in local_sensitivities(hw, recipe).items():
        step = "+2 C" if knob == "temperature_C" else "+2 %"
        print(f"      {knob:<16} {step}: {pct:+.2f} % rate")

    # 4. Scenarios: run charts + datasets
    print("[4/4] Running scenarios (600 wafers each)")
    for i, (name, sc) in enumerate(SCENARIOS.items()):
        _, wafers, sites = run_scenario(sc, seed=42 + i)
        plot_trend(wafers, sc, os.path.join(FIG, f"10_trend_{name}.png"))
        if not quick:
            wafers.to_csv(os.path.join(DATA, f"{name}_wafers.csv"), index=False)
            sites.to_csv(os.path.join(DATA, f"{name}_sites.csv"), index=False)
        late = wafers.tail(50)
        print(f"      {name:<22} late mean {late['m_thickness_mean_nm'].mean():6.1f} nm | "
              f"NU {late['m_nu_1sigma_pct'].mean():.2f} % | label: {wafers['gt_label'].iloc[-1]}")
    print(f"\nDone. Figures in {os.path.abspath(FIG)}")
    if not quick:
        print(f"Datasets in {os.path.abspath(DATA)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip CSV export")
    main(quick=ap.parse_args().quick)
