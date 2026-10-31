# Incident report: chamber CH-A, scenario `mfc_drift`

**Detection.** Multivariate alarm first raised at wafer 183 (lot L0008).

**Current diagnosis (wafer 599).** `mfc_drift`, explaining 99% of the deviation.

**Estimated severity:**
- mfc_drift: -0.0504 (fractional flow error)

**Evidence (largest deviations from healthy, in sigmas):**
- mean_nm: -17.1
- ri_mean: -16.5
- astig_0: +1.5
- edge_delta_pct: +1.0

**Playbook for `mfc_drift`.**
- Expected map pattern: Whole-wafer mean thickness drops with no change in map shape. Refractive index drops too (film becomes more oxygen-rich).
- Sensor check: None on the flow reading. Only metrology (thickness + RI together) points at gas chemistry.
- Action: MFC rate-of-rise verification, recalibrate or replace the MFC.

_Ground truth for scoring (not available in a real fab): mfc_drift._