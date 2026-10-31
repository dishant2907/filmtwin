# Incident report: chamber CH-A, scenario `leveling_tilt`

**Detection.** Multivariate alarm first raised at wafer 302 (lot L0013).

**Current diagnosis (wafer 599).** `leveling_tilt`, explaining 99% of the deviation.

**Estimated severity:**
- leveling_tilt: +0.0152 (tilt per radius)

**Evidence (largest deviations from healthy, in sigmas):**
- tilt_y: +23.9
- tilt_x: -13.4
- s_heater_center_pct: -1.8
- edge_delta_pct: -1.7

**Playbook for `leveling_tilt`.**
- Expected map pattern: Step change right after the maintenance event: one side of the wafer thick, the opposite side thin (planar tilt).
- Sensor check: None. The timing (right after PM) is the main clue.
- Action: Re-level the heater / check lift pins; add a post-PM qualification wafer with a uniformity check.

_Ground truth for scoring (not available in a real fab): leveling_tilt._