"""Module 4c: Run-to-run (R2R) control of deposition time.

Deposition rate drifts: the walls coat up between cleans, parts wear, flows
drift. Without control, thickness drifts with it. An EWMA run-to-run
controller fixes that the way fabs do:

    model         thickness = rate * time
    after a wafer is measured:
                  rate_est <- lam * (measured / time) + (1 - lam) * rate_est
    next wafer:   time = target / rate_est   (clamped to safe limits)

Real-world details included, because they decide whether R2R works:
  - only sampled wafers are measured (slots 1 and 13 of each lot by default),
  - measurements arrive late (metrology queue delay, in wafers),
  - step size and absolute time are clamped,
  - optional reset after a chamber clean, when the rate jumps back up.

The key lesson: a controller that holds thickness on target also HIDES the
drift that caused the problem. The controller's own output (the deposition
time it chose) has to be monitored like any other sensor.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np
import pandas as pd

from .chamber import ChamberSimulator


@dataclass
class EWMAController:
    target_nm: float = 500.0
    lam: float = 0.4
    nominal_time_s: float = 103.0
    initial_rate_nm_s: float = None
    max_step_frac: float = 0.10  # max change in time per adjustment
    time_limits_s: tuple = (80.0, 130.0)
    reset_on_pm: bool = False  # True only helps if the clean actually fixes the drift
    # Feedforward: known extra time for the first wafers after a clean (seasoning
    # effect). Feedback can't fix a wafer before it is measured; feedforward can.
    post_pm_boost: tuple = ()  # e.g. (0.02, 0.005) = +2 % on wafer 0, +0.5 % on wafer 1

    def __post_init__(self):
        self.rate = self.initial_rate_nm_s or self.target_nm / self.nominal_time_s
        self.post_pm_rate = self.rate
        self.time = self.nominal_time_s

    def update(self, measured_nm: float, time_used_s: float, wafers_since_clean: int = None):
        observed = measured_nm / time_used_s
        # A just-cleaned chamber runs slow for a wafer or two. Convert that
        # measurement to its settled-chamber equivalent before learning from
        # it, or the controller over-corrects the following wafers.
        if wafers_since_clean is not None and wafers_since_clean < len(self.post_pm_boost):
            observed *= 1.0 + self.post_pm_boost[wafers_since_clean]
        self.rate = self.lam * observed + (1 - self.lam) * self.rate

    def on_pm(self):
        if self.reset_on_pm:
            self.rate = self.post_pm_rate

    def next_time(self, wafers_since_clean: int = None) -> float:
        desired = self.target_nm / self.rate
        if wafers_since_clean is not None and wafers_since_clean < len(self.post_pm_boost):
            boost = 1.0 + self.post_pm_boost[wafers_since_clean]
            return float(np.clip(desired * boost, *self.time_limits_s))  # does not move the base time
        lo = self.time * (1 - self.max_step_frac)
        hi = self.time * (1 + self.max_step_frac)
        t = float(np.clip(desired, lo, hi))
        t = float(np.clip(t, *self.time_limits_s))
        self.time = t
        return t


def run_closed_loop(sim: ChamberSimulator, n_wafers: int, controller: EWMAController = None,
                    metrology_lag: int = 2, sampled_only: bool = True):
    """Run the chamber wafer by wafer with (or without) R2R control.

    controller=None runs open loop at the recipe's fixed deposition time.
    Returns (wafers_df, sites_df) with extra columns: r2r_time_s, r2r_rate_est.
    """
    recipe = sim.recipe
    pending = []  # (available_at_wafer, measured_nm, time_used)
    wafer_rows, site_rows = [], []
    for i in range(n_wafers):
        if controller is not None:
            if i % sim.pm_interval_wafers == 0 and i > 0:
                controller.on_pm()
            for item in [p for p in pending if p[0] <= i]:
                controller.update(item[1], item[2], item[3])
                pending.remove(item)
            t = controller.next_time(i % sim.pm_interval_wafers)
            rec = dataclasses.replace(recipe, dep_time_s=t)
        else:
            rec = recipe
        row, sites = sim.process_wafer(i, rec)
        row["r2r_time_s"] = rec.dep_time_s
        row["r2r_rate_est"] = controller.rate if controller is not None else np.nan
        if controller is not None and (row["sampled_for_metrology"] or not sampled_only):
            pending.append((i + 1 + metrology_lag, row["m_thickness_mean_nm"], rec.dep_time_s,
                            row["wafers_since_clean"]))
        wafer_rows.append(row)
        site_rows.extend(sites)
    return pd.DataFrame(wafer_rows), pd.DataFrame(site_rows)


def estimate_post_pm_boost(wafers: pd.DataFrame, n: int = 2, settle=(5, 30)):
    """Learn the feedforward boost from historical (open-loop) data: how much
    thinner are the first wafers after a clean than settled wafers?"""
    settled = wafers[wafers.wafers_since_clean.between(*settle)]["m_thickness_mean_nm"].mean()
    out = []
    for k in range(n):
        first = wafers[wafers.wafers_since_clean == k]["m_thickness_mean_nm"].mean()
        out.append(float(settled / first - 1.0))
    return tuple(out)


def thickness_error_stats(wafers: pd.DataFrame, target: float, spec_half_width: float):
    e = wafers["m_thickness_mean_nm"].to_numpy() - target
    return {
        "rmse_nm": float(np.sqrt(np.mean(e ** 2))),
        "mean_error_nm": float(e.mean()),
        "max_abs_error_nm": float(np.max(np.abs(e))),
        "pct_out_of_spec": float(np.mean(np.abs(e) > spec_half_width) * 100),
    }
