# Module 1: what you should be able to explain

These are the questions an interviewer for a deposition or process engineering role is likely to ask. Every answer is something you can show with FilmTwin on screen.

## 1. "Your deposition rate dropped 5 %. Where do you look first?"

Start with the regime. In a reaction-limited process (`Da < 1`), temperature is the prime suspect: heater, thermocouple, wafer-to-heater contact. In a transport-limited process (`Da > 1`), look at gas delivery: MFCs, pressure control, spacing, showerhead.

PECVD sits in between, and plasma adds a third branch: RF delivery (match, generator, cables). Then split the symptom.

- **Did the map shape change, or only the mean?** A mean-only shift points at something global, such as flow, power or chamber wall condition. A shape change points at something local, such as one heater zone, the showerhead pattern or the edge ring.
- **What do the sensors say?** Reflected power up means RF. Edge heater power drifting means temperature control. If no sensor moved but refractive index did, suspect gas chemistry (MFC).

Run `rf_match_degradation` and `mfc_drift` side by side. The thickness maps are nearly identical, but RI and reflected power separate them.

## 2. "Explain the Grove model."

Two fluxes in series. Gas-phase transport brings precursor to the surface (`h_g`). The surface reaction consumes it (`k_s`). At steady state the effective rate constant is `k_s*h_g/(k_s+h_g)`, which works like two resistances in series, so the slower step controls the rate. This is the same mass-transfer-with-reaction math as a catalytic reactor in chemical engineering.

## 3. "Why is the wafer centre thicker than the edge?"

In this chamber the cause is precursor depletion as gas moves outward, a slightly centre-high plasma, and edge heat loss. Typical knobs to fix it:

- an edge heater trim (`edge_temp_offset_C`),
- more total flow (less depletion),
- a different spacing.

Every fix is a trade-off against rate and film properties. That trade-off is exactly what the Module 4 DOE will map.

## 4. "Uniformity got better. Is that good news?"

Not necessarily. In `showerhead_clog` the 1-sigma NU improves for about 200 wafers because the fault cancels the chamber's natural centre-thick shape, and then it gets worse. An SPC chart on NU alone would miss the start. Monitoring the *shape* (radial and Zernike terms, Module 3) catches it early.

## 5. "What happens after a chamber clean?"

The first wafer runs low because the chamber isn't seasoned yet: fresh walls absorb radicals. Then the walls coat up and the rate slowly declines until the next clean. That is why fabs run seasoning wafers and why SPC limits often depend on wafers-since-clean.

Post-maintenance problems such as `leveling_tilt` show up as step changes right after the clean. That's why a post-PM qualification wafer matters.

## 6. "How would you validate a model like this?"

Be honest. Constants are set to realistic windows, not calibrated to a tool. The next step is calibration against real data: published rate-vs-recipe trends first, then real wafers when available. Show you know the difference between a model that has the right *directions* and one that predicts *absolute values*.
