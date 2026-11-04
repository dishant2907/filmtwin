# DOE report: PECVD oxide recipe optimisation on FilmTwin

## 1. Screening (2^(6-2) resolution IV, 16 runs + 4 centre points, randomised)

**rate_nm_min** (R² 0.979, model p 4.2e-09)

| term | effect | t | p |
|---|---|---|---|
| temperature_C | +28.41 | +4.48 | 0.000746 |
| pressure_torr | +77.59 | +12.25 | 3.86e-08 |
| rf_power_W | +44.15 | +6.97 | 1.5e-05 |
| precursor_sccm | +110.3 | +17.40 | 7.03e-10 |
| spacing_mm | -37.78 | -5.96 | 6.58e-05 |
| edge_temp_offset_C | +2.827 | +0.45 | 0.663 |
| curvature(center) |  | +0.68 | 0.509 |

**nu_pct** (R² 0.986, model p 3.8e-10)

| term | effect | t | p |
|---|---|---|---|
| temperature_C | +0.06422 | +2.14 | 0.0535 |
| pressure_torr | +0.2463 | +8.21 | 2.87e-06 |
| rf_power_W | +0.1346 | +4.49 | 0.000743 |
| precursor_sccm | -0.03639 | -1.21 | 0.248 |
| spacing_mm | -0.1337 | -4.46 | 0.000781 |
| edge_temp_offset_C | -0.8083 | -26.95 | 4.19e-12 |
| curvature(center) |  | -0.09 | 0.929 |

**ri** (R² 1.000, model p 2.8e-21)

| term | effect | t | p |
|---|---|---|---|
| temperature_C | -0.00199 | -30.48 | 9.77e-13 |
| pressure_torr | +9.916e-05 | +1.52 | 0.155 |
| rf_power_W | +5.444e-05 | +0.83 | 0.421 |
| precursor_sccm | +0.01603 | +245.55 | 1.4e-23 |
| spacing_mm | +4.742e-05 | +0.73 | 0.482 |
| edge_temp_offset_C | -0.0001921 | -2.94 | 0.0123 |
| curvature(center) |  | +1.28 | 0.224 |

**werr_mean** (R² 0.998, model p 4.1e-15)

| term | effect | t | p |
|---|---|---|---|
| temperature_C | -1.033 | -74.83 | 2.16e-17 |
| pressure_torr | -0.0008892 | -0.06 | 0.95 |
| rf_power_W | -0.1104 | -7.99 | 3.79e-06 |
| precursor_sccm | -0.001027 | -0.07 | 0.942 |
| spacing_mm | -0.000233 | -0.02 | 0.987 |
| edge_temp_offset_C | -0.0954 | -6.91 | 1.63e-05 |
| curvature(center) |  | -5.85 | 7.85e-05 |

Vital few chosen for the response surface: edge_temp_offset_C, pressure_torr, precursor_sccm.

## 2. Response surface (face-centred CCD, 20 runs)

| response | R² | adj R² | RMSE |
|---|---|---|---|
| rate_nm_min | 0.9998 | 0.9997 | 0.8597 |
| nu_pct | 0.9682 | 0.9396 | 0.09596 |
| profile_pct | 0.9962 | 0.9928 | 0.1317 |
| ri | 0.9999 | 0.9997 | 9.553e-05 |
| werr_mean | 0.9990 | 0.9982 | 0.002733 |
| werr_edge_delta | 0.9947 | 0.9900 | 0.01063 |

Uniformity is modelled through the signed edge-minus-centre profile, not 1-sigma NU. NU is a standard deviation and folds into a V at its minimum, which a quadratic cannot fit (compare the R² values).

## 3. Optimum (Derringer desirability)

| factor | value |
|---|---|
| edge_temp_offset_C | 26.2 |
| pressure_torr | 3 |
| precursor_sccm | 152 |
| overall desirability | 0.934 |

Pressure sits at the edge of the explored range. The next experiment should extend it, checking constraints this model does not include (particles from gas-phase nucleation at high pressure, film stress).

## 4. Confirmation (5 wafers)

Verdict compares the mean of the confirmation runs with the 95 % prediction interval for a mean of 5 new runs. The last column puts the miss in engineering terms: the difference as a share of the window that matters for that response.

| response | predicted | 95 % PI | observed mean | verdict | miss as % of window |
|---|---|---|---|---|---|
| profile_pct | -0.0159944 | ±0.204 | -0.17999 | confirmed | 5.5 % (±3) |
| rate_nm_min | 336.32 | ±1.33 | 335.554 | confirmed | 7.7 % (±10) |
| ri | 1.47829 | ±0.000148 | 1.47812 | outside PI | 2.2 % (±0.008) |
| werr_mean | 1.40277 | ±0.00422 | 1.40546 | confirmed | 2.7 % (±0.1) |
| nu_pct | 0.286133 | ±0.148 | 0.197478 | confirmed | 17.7 % (±0.5) |

**Reading this table like an engineer.** A response can fail the statistical test and still not matter. Refractive index is measured so precisely that its prediction interval is tiny, so a model error far too small to affect the film becomes statistically detectable. Also, with five responses each tested at 95 %, there is roughly a 23 % chance that at least one fails by chance alone (1 - 0.95^5). The decision rule: a statistical miss that is also a small fraction of the engineering window is noted, not acted on; a miss that is a large fraction of the window means the model is wrong in that region and needs more experiments.

**Result:** 1-sigma non-uniformity 1.27 % → 0.20 %, deposition rate 292 → 336 nm/min, refractive index held on target.

Raw tables (`screening_results.csv`, `ccd_results.csv`) can be re-analysed in JMP: Analyze → Fit Model, with the coded factor columns as effects.