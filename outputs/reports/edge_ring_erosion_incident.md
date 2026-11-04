# Incident report: chamber CH-A, scenario `edge_ring_erosion`

**Detection.** Multivariate alarm first raised at wafer 40 (lot L0002).

**Current diagnosis (wafer 599).** `edge_ring_erosion`, explaining 100% of the deviation.

**Estimated severity:**
- edge_ring_erosion: +0.0955 (extra edge plasma)

**Evidence (largest deviations from healthy, in sigmas):**
- edge_delta_pct: +56.3
- defocus: +45.4
- rate_nm_min: +10.6
- werr_edge_delta: -7.6

**Playbook for `edge_ring_erosion`.**
- Expected map pattern: The outer edge goes from slightly thin to thick. Like the showerhead case, NU first improves (the wear cancels the chamber's natural edge-thin shape) and then gets worse as a thick edge ring forms.
- Sensor check: Weak; nothing clear in the summary sensors.
- Action: Replace the edge ring; track ring RF-hours against the uniformity trend.

_Ground truth for scoring (not available in a real fab): edge_ring_erosion._