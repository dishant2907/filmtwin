"""Scenario library: ready-made fault stories with known root causes.

Each scenario is a run of wafers (default 600 = 24 lots of 25) with chamber
cleans every 300 wafers and one or more faults. The descriptions are what a
process engineer would write in an escalation summary once the root cause is
known. Use them to check your own diagnosis, not as hints before you try.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .faults import FaultEvent


@dataclass
class Scenario:
    name: str
    description: str
    map_signature: str
    sensor_clue: str
    typical_fix: str
    faults: list = field(default_factory=list)
    n_wafers: int = 600
    pm_interval_wafers: int = 300


SCENARIOS = {
    "baseline": Scenario(
        name="baseline",
        description="Healthy chamber. Only normal aging: slow rate loss as the walls coat up, and a low first wafer after each clean.",
        map_signature="Stable centre-slightly-thick bowl, ~1.2-1.5 % 1-sigma.",
        sensor_clue="Vpp rises slowly between cleans and resets at each clean.",
        typical_fix="None. Use this as your reference for SPC limits.",
    ),
    "edge_tc_drift": Scenario(
        name="edge_tc_drift",
        description="Edge-zone thermocouple drifts high. The controller cuts edge heater power, so the real wafer edge runs up to 15 C colder.",
        map_signature="Outer ring goes ~3 % thin (PECVD rate is only weakly temperature-dependent, so the thickness hit is modest), while edge wet etch rate ratio climbs. The film-quality change is the bigger risk.",
        sensor_clue="Edge heater power % trends down while the edge temperature reading looks perfect.",
        typical_fix="Replace or recalibrate the edge thermocouple; verify with a temperature-sensing wafer.",
        faults=[FaultEvent("edge_tc_drift", start_wafer=120, severity=15.0, ramp_wafers=300)],
    ),
    "showerhead_clog": Scenario(
        name="showerhead_clog",
        description="Deposits slowly block centre showerhead holes, starving the wafer centre.",
        map_signature="Centre goes thin. The normal centre-thick bowl first flattens (1-sigma NU actually IMPROVES for a while, so an NU-only SPC chart misses it) and then inverts into a centre-thin dish. Cleans do not fix it.",
        sensor_clue="Almost none. Throttle valve and pressure look normal.",
        typical_fix="Inspect and replace or clean the showerhead (faceplate).",
        faults=[FaultEvent("showerhead_clog", start_wafer=50, severity=0.12, ramp_wafers=500)],
    ),
    "mfc_drift": Scenario(
        name="mfc_drift",
        description="Precursor MFC calibration drifts low by up to 5 %. It still reports its setpoint.",
        map_signature="Whole-wafer mean thickness drops with no change in map shape. Refractive index drops too (film becomes more oxygen-rich).",
        sensor_clue="None on the flow reading. Only metrology (thickness + RI together) points at gas chemistry.",
        typical_fix="MFC rate-of-rise verification, recalibrate or replace the MFC.",
        faults=[FaultEvent("mfc_drift", start_wafer=100, severity=-0.05, ramp_wafers=400)],
    ),
    "rf_match_degradation": Scenario(
        name="rf_match_degradation",
        description="RF matching network degrades; up to 10 % of forward power is reflected or lost.",
        map_signature="Mean thickness drops, uniformity changes slightly. RI unchanged.",
        sensor_clue="Reflected power climbs and Vpp drops. Forward power reading stays at setpoint.",
        typical_fix="Inspect the match (capacitors, connections), retune or replace.",
        faults=[FaultEvent("rf_match_degradation", start_wafer=200, severity=0.10, ramp_wafers=250)],
    ),
    "edge_ring_erosion": Scenario(
        name="edge_ring_erosion",
        description="The process kit edge ring wears over its life, concentrating plasma at the wafer edge.",
        map_signature="The outer edge goes from slightly thin to thick. Like the showerhead case, NU first improves (the wear cancels the chamber's natural edge-thin shape) and then gets worse as a thick edge ring forms.",
        sensor_clue="Weak; nothing clear in the summary sensors.",
        typical_fix="Replace the edge ring; track ring RF-hours against the uniformity trend.",
        faults=[FaultEvent("edge_ring_erosion", start_wafer=0, severity=0.10, ramp_wafers=600)],
    ),
    "leveling_tilt": Scenario(
        name="leveling_tilt",
        description="After the clean at wafer 300, the heater was re-installed slightly out of level.",
        map_signature="Step change right after the maintenance event: one side of the wafer thick, the opposite side thin (planar tilt).",
        sensor_clue="None. The timing (right after PM) is the main clue.",
        typical_fix="Re-level the heater / check lift pins; add a post-PM qualification wafer with a uniformity check.",
        faults=[FaultEvent("leveling_tilt", start_wafer=300, severity=0.015, angle_deg=120.0)],
    ),
    "mixed_clog_and_match": Scenario(
        name="mixed_clog_and_match",
        description="Two problems at once: slow showerhead clogging plus a degrading RF match.",
        map_signature="Mean drops and the map shape changes. Separating the two requires using sensors and map shape together.",
        sensor_clue="Reflected power explains the mean drop; the shape change has no sensor clue.",
        typical_fix="Fix the match first (sensor-confirmed), then re-check the map shape for the remaining problem.",
        faults=[
            FaultEvent("showerhead_clog", start_wafer=60, severity=0.12, ramp_wafers=500),
            FaultEvent("rf_match_degradation", start_wafer=250, severity=0.08, ramp_wafers=250),
        ],
    ),
}


def get(name: str) -> Scenario:
    if name not in SCENARIOS:
        raise KeyError(f"Unknown scenario '{name}'. Available: {list(SCENARIOS)}")
    return SCENARIOS[name]
