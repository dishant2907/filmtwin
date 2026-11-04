# Incident report: chamber CH-A, scenario `showerhead_clog`

**Detection.** Multivariate alarm first raised at wafer 80 (lot L0004).

**Current diagnosis (wafer 599).** `showerhead_clog`, explaining 100% of the deviation.

**Estimated severity:**
- showerhead_clog: +0.121 (fraction of centre gas lost)

**Evidence (largest deviations from healthy, in sigmas):**
- defocus: +67.6
- centre_delta_pct: -38.7
- rate_nm_min: -4.1
- s_vpp_V: -1.7

**Playbook for `showerhead_clog`.**
- Expected map pattern: Centre goes thin. The normal centre-thick bowl first flattens (1-sigma NU actually IMPROVES for a while, so an NU-only SPC chart misses it) and then inverts into a centre-thin dish. Cleans do not fix it.
- Sensor check: Almost none. Throttle valve and pressure look normal.
- Action: Inspect and replace or clean the showerhead (faceplate).

_Ground truth for scoring (not available in a real fab): showerhead_clog._