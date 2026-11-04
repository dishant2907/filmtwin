# Incident report: chamber CH-A, scenario `rf_match_degradation`

**Detection.** Multivariate alarm first raised at wafer 211 (lot L0009).

**Current diagnosis (wafer 599).** `rf_match_degradation`, explaining 100% of the deviation.

**Estimated severity:**
- rf_match_degradation: +0.0993 (fraction of RF power lost)

**Evidence (largest deviations from healthy, in sigmas):**
- s_rf_reflected_W: +128.8
- rate_nm_min: -17.6
- s_vpp_V: -10.0
- s_heater_edge_pct: +2.3

**Playbook for `rf_match_degradation`.**
- Expected map pattern: Mean thickness drops, uniformity changes slightly. RI unchanged.
- Sensor check: Reflected power climbs and Vpp drops. Forward power reading stays at setpoint.
- Action: Inspect the match (capacitors, connections), retune or replace.

_Ground truth for scoring (not available in a real fab): rf_match_degradation._