# OCAP execution log: MFC drift excursion

## The plan

```
OCAP: PECVD oxide mean thickness, chamber CH-A

Trigger: any Western Electric rule on the mean-thickness I chart, or an EWMA signal.

1. VERIFY THE MEASUREMENT. Re-measure the wafer. If the re-measurement is
   back inside the limits and the fit GOF is good, log a metrology false
   alarm and release the lot.
2. CHECK CONTEXT. Is this the first or second wafer after a chamber clean?
   If so, a low reading is the known seasoning effect: log it, continue.
3. CHECK TOOL SENSORS (FDC). Reflected RF power, edge and centre heater power,
   throttle valve. Any sensor outside its own limits points at hardware.
4. RUN THE DIAGNOSER (Module 3) on the last wafers. Record the probable
   cause, its estimated severity, and the share of the deviation explained.
5. DISPOSITION.
   - Out of spec: HOLD the lot, take the chamber down, escalate to the
     process owner. Restart only after a qualification wafer passes.
   - In spec, nothing confirms the chart (sensors normal, no multivariate
     anomaly): continue under watch, sample every wafer, re-run if it persists.
     Charts with several run rules do raise occasional false alarms.
   - In spec, cause identified: continue running, open a maintenance ticket
     for the diagnosed part, tighten sampling to every wafer.
   - In spec but an anomaly nobody can explain: HOLD and escalate.
```

## Trigger at wafer 116 (trend rule)

- Step 1, verify: re-measured 496.7 nm (first reading 496.6 nm). Confirms the signal.
- Step 2, context: 116 wafers since clean, not a seasoning effect.
- Step 3, sensors: all inside their limits.
- Step 4, diagnosis: none (explains 67% of the deviation).
- Step 5, disposition: CONTINUE UNDER WATCH: single-chart signal with no confirming evidence (likely a chart false alarm or a very early drift); sample every wafer for the next lot and re-run this OCAP if the signal persists.

**Disposition:** CONTINUE UNDER WATCH: single-chart signal with no confirming evidence (likely a chart false alarm or a very early drift); sample every wafer for the next lot and re-run this OCAP if the signal persists.

## Trigger at wafer 260 (we1)

- Step 1, verify: re-measured 489.0 nm (first reading 489.0 nm). Confirms the signal.
- Step 2, context: 260 wafers since clean, not a seasoning effect.
- Step 3, sensors: all inside their limits.
- Step 4, diagnosis: mfc_drift (explains 92% of the deviation).
- Step 5, disposition: CONTINUE with ticket: thickness 489.0 nm still in spec; maintenance ticket for mfc_drift; sample every wafer until fixed.

**Disposition:** CONTINUE with ticket: thickness 489.0 nm still in spec; maintenance ticket for mfc_drift; sample every wafer until fixed.

## Trigger at wafer 450 (we1)

- Step 1, verify: re-measured 480.1 nm (first reading 480.0 nm). Confirms the signal.
- Step 2, context: 150 wafers since clean, not a seasoning effect.
- Step 3, sensors: all inside their limits.
- Step 4, diagnosis: mfc_drift (explains 96% of the deviation).
- Step 5, disposition: HOLD LOT, chamber down (out of spec); escalate to process owner.

**Disposition:** HOLD LOT, chamber down (out of spec); escalate to process owner.
