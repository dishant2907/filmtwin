"""Chamber simulator: runs wafers through a chamber that ages and develops faults.

Produces two tables, the same shape a fab data system would give you:

  wafers  one row per wafer: context (lot, slot, time), recipe, tool sensor
          summaries (FDC-style), metrology summary, and ground-truth fault
          labels (columns starting with `gt_`). In a real fab you never get the
          `gt_` columns; here they let you score your root-cause methods.
  sites   one row per measured site: x, y, thickness, refractive index, wet
          etch rate ratio, local temperature.
"""
from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field
from datetime import datetime, timedelta

import numpy as np
import pandas as pd

from .faults import FAULT_KINDS, Aging, FaultEvent, build_condition
from .metrics import uniformity
from .physics import Condition, Hardware, Recipe, deposit, gas_state
from .wafer import Wafer, sites_49pt

WAFERS_PER_LOT = 25


@dataclass
class SensorModel:
    """Turns the true chamber state into what the tool's sensors report."""

    tv_ref_pct: float = 38.0  # throttle valve opening at the reference flow/pressure
    vpp_ref_V: float = 290.0
    heater_center_ref_pct: float = 52.0
    heater_edge_ref_pct: float = 61.0

    def read(self, recipe: Recipe, cond: Condition, hw: Hardware, gas: dict, rng) -> dict:
        ref_flow = 3150.0
        delivered = recipe.rf_power_W * cond.rf_delivered_factor
        power_lost = recipe.rf_power_W - delivered
        wall_coating = 1.0 - cond.radical_factor  # rough proxy for wall deposits
        return {
            "s_pressure_torr": recipe.pressure_torr + rng.normal(0, 0.004),
            # Valve opens wider for more gas, closes for higher pressure setpoints.
            "s_throttle_valve_pct": self.tv_ref_pct
            * (gas["total_actual"] / ref_flow)
            * (2.5 / recipe.pressure_torr)
            + rng.normal(0, 0.25),
            "s_rf_forward_W": recipe.rf_power_W + rng.normal(0, 1.0),
            # Most of the lost power shows up as reflected power.
            "s_rf_reflected_W": 2.0 + 0.8 * power_lost + abs(rng.normal(0, 0.6)),
            # Vpp scales with delivered power and shifts as the walls coat up.
            "s_vpp_V": self.vpp_ref_V
            * np.sqrt(max(delivered, 1.0) / 600.0)
            * (1.0 + 2.0 * wall_coating)
            + rng.normal(0, 1.5),
            "s_heater_center_pct": self.heater_center_ref_pct
            + 0.10 * (recipe.temperature_C - 400.0)
            + rng.normal(0, 0.3),
            # Controller trusts a thermocouple that reads high, so it cuts power.
            "s_heater_edge_pct": self.heater_edge_ref_pct
            + 0.10 * (recipe.temperature_C - 400.0)
            + 0.12 * recipe.edge_temp_offset_C
            + 0.9 * cond.edge_temp_error_C
            + rng.normal(0, 0.3),
            # The MFC reports its setpoint even when its calibration has drifted.
            "s_precursor_flow_sccm": recipe.precursor_sccm + rng.normal(0, 0.3),
            "s_wall_temp_C": 65.0 + rng.normal(0, 0.4),
        }


@dataclass
class ChamberSimulator:
    """A single chamber processing wafers over time."""

    hardware: Hardware = field(default_factory=Hardware)
    recipe: Recipe = field(default_factory=Recipe)
    faults: list = field(default_factory=list)
    aging: Aging = field(default_factory=Aging)
    sensors: SensorModel = field(default_factory=SensorModel)
    wafer: Wafer = field(default_factory=Wafer)
    pm_interval_wafers: int = 300  # chamber wet clean every N wafers
    metrology_slots: tuple = (1, 13)  # slots measured per lot (sampling plan)
    seed: int = 0
    start_time: datetime = datetime(2026, 1, 5, 6, 0, 0)
    cycle_time_s: float = 210.0  # deposition plus transfer and pump/purge

    def __post_init__(self):
        self.rng = np.random.default_rng(self.seed)
        self.sites = sites_49pt(self.wafer)

    # -- single wafer --------------------------------------------------------
    def condition_at(self, wafer_index: int):
        """Chamber condition and ground truth at a given wafer count."""
        wafers_since_clean = wafer_index % self.pm_interval_wafers
        last_pm = wafer_index - wafers_since_clean
        return build_condition(self.faults, wafer_index, wafers_since_clean, last_pm, self.aging)

    def process_wafer(self, i: int, recipe: Recipe):
        """Process one wafer at wafer count i. Returns (wafer_row, site_rows).

        Call wafers in order (0, 1, 2, ...). This is what closed-loop control
        (Module 4) uses: decide the recipe, run one wafer, look at the result,
        decide the next recipe.
        """
        wafers_since_clean = i % self.pm_interval_wafers
        if wafers_since_clean == 0 or not hasattr(self, "_rf_seconds"):
            self._rf_seconds = 0.0
        cond, truth = self.condition_at(i)

        res = deposit(self.sites, recipe, self.hardware, cond, self.rng,
                      usable_radius_mm=self.wafer.usable_radius_mm)
        gas = gas_state(recipe, self.hardware, cond, float(np.mean(res.temperature_C) + 273.15))
        sensors = self.sensors.read(recipe, cond, self.hardware, gas, self.rng)
        self._rf_seconds += recipe.dep_time_s

        lot = i // WAFERS_PER_LOT + 1
        slot = i % WAFERS_PER_LOT + 1
        wafer_id = f"L{lot:04d}-W{slot:02d}"
        u = uniformity(res.thickness_nm)
        active = [k for k in FAULT_KINDS if abs(truth[k]) > 0]

        row = {
            "wafer_index": i,
            "wafer_id": wafer_id,
            "lot_id": f"L{lot:04d}",
            "slot": slot,
            "chamber": self.hardware.name,
            "timestamp": self.start_time + timedelta(seconds=i * self.cycle_time_s),
            "wafers_since_clean": wafers_since_clean,
            "rf_hours_since_clean": self._rf_seconds / 3600.0,
            "pm_event": wafers_since_clean == 0 and i > 0,
            "sampled_for_metrology": slot in self.metrology_slots,
        }
        row.update({f"r_{k}": v for k, v in dataclasses.asdict(recipe).items()})
        row.update(sensors)
        row.update({
            "m_thickness_mean_nm": u["mean"],
            "m_nu_1sigma_pct": u["nu_1sigma_pct"],
            "m_nu_range_pct": u["nu_range_pct"],
            "m_ri_mean": float(np.mean(res.refractive_index)),
            "m_werr_mean": float(np.mean(res.werr)),
            "m_rate_nm_min": float(np.mean(res.rate_nm_min)),
        })
        row.update({f"gt_{k}": truth[k] for k in FAULT_KINDS})
        row["gt_label"] = "+".join(active) if active else "none"

        site_rows = [{
            "wafer_index": i,
            "wafer_id": wafer_id,
            "site": s + 1,
            "x_mm": x,
            "y_mm": y,
            "thickness_nm": res.thickness_nm[s],
            "refractive_index": res.refractive_index[s],
            "werr": res.werr[s],
            "temperature_C": res.temperature_C[s],
        } for s, (x, y) in enumerate(self.sites)]
        return row, site_rows

    def run(self, n_wafers: int, recipe_changes: dict = None):
        """Process n_wafers and return (wafers_df, sites_df).

        recipe_changes: optional {wafer_index: Recipe} to switch recipes mid-run.
        """
        recipe_changes = recipe_changes or {}
        recipe = self.recipe
        wafer_rows, site_rows = [], []
        for i in range(n_wafers):
            if i in recipe_changes:
                recipe = recipe_changes[i]
            row, sites = self.process_wafer(i, recipe)
            wafer_rows.append(row)
            site_rows.extend(sites)
        return pd.DataFrame(wafer_rows), pd.DataFrame(site_rows)

    # -- ideal maps for plotting --------------------------------------------
    def ideal_map(self, xy: np.ndarray, wafer_index: int = 50, recipe: Recipe = None):
        """Noise-free film at arbitrary points for a given point in the run."""
        cond, _ = self.condition_at(wafer_index)
        return deposit(xy, recipe or self.recipe, self.hardware, cond, noise=False,
                       usable_radius_mm=self.wafer.usable_radius_mm)


def run_scenario(scenario, n_wafers: int = None, seed: int = 0, hardware: Hardware = None,
                 recipe: Recipe = None):
    """Convenience wrapper: build a simulator from a Scenario and run it."""
    sim = ChamberSimulator(
        hardware=hardware or Hardware(),
        recipe=recipe or Recipe(),
        faults=list(scenario.faults),
        pm_interval_wafers=scenario.pm_interval_wafers,
        seed=seed,
    )
    wafers, sites = sim.run(n_wafers or scenario.n_wafers)
    wafers.insert(0, "scenario", scenario.name)
    sites.insert(0, "scenario", scenario.name)
    return sim, wafers, sites
