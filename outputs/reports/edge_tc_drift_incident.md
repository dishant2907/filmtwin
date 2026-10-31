# Incident report: chamber CH-A, scenario `edge_tc_drift`

**Detection.** Multivariate alarm first raised at wafer 146 (lot L0006).

**Current diagnosis (wafer 599).** `edge_tc_drift`, explaining 100% of the deviation.

**Estimated severity:**
- edge_tc_drift: +15.6 (C colder)

**Evidence (largest deviations from healthy, in sigmas):**
- s_heater_edge_pct: -47.3
- defocus: -18.3
- werr_edge_delta: +16.6
- edge_delta_pct: -12.5

**Playbook for `edge_tc_drift`.**
- Expected map pattern: Outer ring goes ~3 % thin (PECVD rate is only weakly temperature-dependent, so the thickness hit is modest), while edge wet etch rate ratio climbs. The film-quality change is the bigger risk.
- Sensor check: Edge heater power % trends down while the edge temperature reading looks perfect.
- Action: Replace or recalibrate the edge thermocouple; verify with a temperature-sensing wafer.

_Ground truth for scoring (not available in a real fab): edge_tc_drift._