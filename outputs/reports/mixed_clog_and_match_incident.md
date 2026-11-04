# Incident report: chamber CH-A, scenario `mixed_clog_and_match`

**Detection.** Multivariate alarm first raised at wafer 89 (lot L0004).

**Current diagnosis (wafer 599).** `rf_match_degradation+showerhead_clog`, explaining 100% of the deviation.

**Estimated severity:**
- showerhead_clog: +0.121 (fraction of centre gas lost)
- rf_match_degradation: +0.0793 (fraction of RF power lost)

**Evidence (largest deviations from healthy, in sigmas):**
- s_rf_reflected_W: +103.2
- defocus: +66.9
- centre_delta_pct: -40.0
- rate_nm_min: -15.8

**Playbook for `rf_match_degradation`.**
- Expected map pattern: Mean thickness drops, uniformity changes slightly. RI unchanged.
- Sensor check: Reflected power climbs and Vpp drops. Forward power reading stays at setpoint.
- Action: Inspect the match (capacitors, connections), retune or replace.

**Playbook for `showerhead_clog`.**
- Expected map pattern: Centre goes thin. The normal centre-thick bowl first flattens (1-sigma NU actually IMPROVES for a while, so an NU-only SPC chart misses it) and then inverts into a centre-thin dish. Cleans do not fix it.
- Sensor check: Almost none. Throttle valve and pressure look normal.
- Action: Inspect and replace or clean the showerhead (faceplate).

_Ground truth for scoring (not available in a real fab): showerhead_clog+rf_match_degradation._