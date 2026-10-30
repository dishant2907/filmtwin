"""Chamber faults and normal aging.

Each fault changes one part of the chamber `Condition`. Some faults are
visible in tool sensors, some only show up in metrology. That split is
deliberate: it mirrors real root-cause work, where the sensor data and the
wafer map have to be read together.

Fault kinds (severity units in brackets):

  edge_tc_drift          Edge-zone thermocouple reads high, so the heater
                         controller cuts power and the real wafer edge runs
                         colder. [deg C of real edge cooling]
                         Sensor clue: edge heater power % drifts down.
  showerhead_clog        Centre showerhead holes partly blocked by deposits,
                         starving the wafer centre of precursor.
                         [fraction of centre gas lost]
                         Sensor clue: almost none. Map clue: centre goes thin.
  mfc_drift              Precursor mass-flow controller calibration drifts.
                         The MFC still *reports* its setpoint.
                         [fractional flow error, e.g. -0.04]
                         Sensor clue: essentially none. Map clue: whole-wafer
                         mean shift plus a refractive index shift.
  rf_match_degradation   Matching network degrades, part of the forward power
                         is reflected or lost. [fraction of power lost]
                         Sensor clue: reflected power climbs, Vpp drops.
  edge_ring_erosion      Process kit edge ring wears, plasma concentrates at
                         the wafer edge. [extra edge plasma amplitude]
                         Sensor clue: weak. Map clue: thick edge ring.
  leveling_tilt          Heater or lift-pin leveling is off (often right
                         after a maintenance event). [fractional tilt per
                         radius] Sensor clue: none. Map clue: one side thick.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .physics import Condition

FAULT_KINDS = (
    "edge_tc_drift",
    "showerhead_clog",
    "mfc_drift",
    "rf_match_degradation",
    "edge_ring_erosion",
    "leveling_tilt",
)


@dataclass
class FaultEvent:
    """A fault that starts at a given wafer count and ramps to full severity."""

    kind: str
    start_wafer: int
    severity: float
    ramp_wafers: int = 0  # 0 means a step change
    angle_deg: float = 0.0  # only used by leveling_tilt
    cleared_by_pm: bool = False  # True if a chamber clean fixes it

    def __post_init__(self):
        if self.kind not in FAULT_KINDS:
            raise ValueError(f"Unknown fault kind '{self.kind}'. Choose from {FAULT_KINDS}")

    def progress(self, wafer_index: int, last_pm_wafer: int) -> float:
        """0 before the fault starts, rising to 1 at full severity."""
        if wafer_index < self.start_wafer:
            return 0.0
        if self.cleared_by_pm and last_pm_wafer > self.start_wafer:
            return 0.0
        if self.ramp_wafers <= 0:
            return 1.0
        return float(min(1.0, (wafer_index - self.start_wafer + 1) / self.ramp_wafers))


@dataclass
class Aging:
    """Normal, expected chamber aging between cleans (not a fault)."""

    wall_loss: float = 0.015  # radical supply lost as walls coat up
    wall_tau_wafers: float = 150.0
    first_wafer_dip: float = 0.02  # first wafer after a clean runs low
    first_wafer_tau: float = 1.5

    def radical_factor(self, wafers_since_clean: int) -> float:
        n = wafers_since_clean
        wall = 1.0 - self.wall_loss * (1.0 - np.exp(-n / self.wall_tau_wafers))
        seasoning = 1.0 - self.first_wafer_dip * np.exp(-n / self.first_wafer_tau)
        return float(wall * seasoning)


def build_condition(
    faults,
    wafer_index: int,
    wafers_since_clean: int,
    last_pm_wafer: int,
    aging: Aging,
):
    """Combine aging and every active fault into one Condition.

    Returns (condition, ground_truth) where ground_truth maps each fault kind
    to its current effective severity (0 when inactive).
    """
    cond = Condition(radical_factor=aging.radical_factor(wafers_since_clean))
    truth = {k: 0.0 for k in FAULT_KINDS}

    for f in faults:
        p = f.progress(wafer_index, last_pm_wafer)
        if p <= 0:
            continue
        s = f.severity * p
        truth[f.kind] += s
        if f.kind == "edge_tc_drift":
            cond.edge_temp_error_C -= s
        elif f.kind == "showerhead_clog":
            cond.center_gas_loss += s
        elif f.kind == "mfc_drift":
            cond.precursor_flow_error += s
        elif f.kind == "rf_match_degradation":
            cond.rf_delivered_factor *= (1.0 - s)
        elif f.kind == "edge_ring_erosion":
            cond.extra_edge_plasma += s
        elif f.kind == "leveling_tilt":
            cond.extra_tilt += s
            cond.extra_tilt_angle_deg = f.angle_deg

    return cond, truth
